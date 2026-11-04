"""Run reproducible morphology experiments; no changes to ai4alscn."""
import argparse
from dataclasses import replace, asdict
from pathlib import Path
import json, hashlib, csv, time, zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from kmc_model import Config, Model

ROOT = Path(__file__).resolve().parent
COLORS = ListedColormap(['#eef2f8', '#38a69b', '#d09026'])

def panel(ax, state, cfg, title):
    h = cfg.spacing_nm
    depth = 3*cfg.layer_cells*h
    # Electrode rows are displayed as strips outside the active 60 nm film.
    ax.imshow(state[1:-1], cmap=COLORS, vmin=0, vmax=2, interpolation='nearest',
              extent=[0, cfg.nx*h, depth, 0], aspect='equal')
    for y in [depth/3, 2*depth/3]: ax.axhline(y, color='#696c90', ls='--', lw=.8)
    ax.axhspan(-2*h, 0, color='#8896a5', clip_on=False)
    ax.axhspan(depth, depth+2*h, color='#5c6573', clip_on=False)
    ax.set(xlabel='x (nm)', ylabel='Depth from Ag (nm)', title=title, ylim=(depth+2*h, -2*h))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', default='results/baseline')
    ap.add_argument('--config', default='configs/baseline.json')
    ap.add_argument('--seeds', default='7,17,27')
    ap.add_argument('--controls', action='store_true')
    args = ap.parse_args()
    cfg = Config(**json.loads((ROOT/args.config).read_text(encoding='utf-8')))
    out = ROOT/args.output
    out.mkdir(parents=True, exist_ok=False)
    sources = [ROOT/'kmc_model.py', ROOT/'run.py', ROOT/'test_physics.py', ROOT/args.config]
    manifest = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    with zipfile.ZipFile(out/'source_snapshot.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for p in sources: archive.write(p, str(p.relative_to(ROOT)))
    variants = [('assist', cfg.polarization_C_m2)]
    if args.controls: variants += [('off', 0.0), ('oppose', -cfg.polarization_C_m2)]
    summaries, runs = [], []
    for name, p in variants:
        for seed in map(int, args.seeds.split(',')):
            c = replace(cfg, seed=seed, polarization_C_m2=p)
            start = time.perf_counter()
            model = Model(c)
            initial_phi = model.phi.reshape(model.shape).copy()
            summary = model.run()
            summary['potential_state'] = 'before_bridge_event' if model.connected else 'final_state'
            summary['wall_s'] = time.perf_counter()-start
            summary['variant'] = name
            prefix = out/f'{name}_seed_{seed}'
            prefix.with_suffix('.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
            np.savez_compressed(prefix.with_suffix('.npz'), snapshots=np.array(model.snapshots),
                                times_s=[r['time_s'] for r in model.trace],
                                potential_V=model.phi.reshape(model.shape),
                                initial_potential_V=initial_phi,
                                site_energy_eV=model.site.reshape(model.shape))
            summaries.append(summary); runs.append((model, summary))
            print(name, seed, summary['status'], summary['events'], f"{summary['time_s']:.5g}s", model.trace[-1], flush=True)
    fig, axs = plt.subplots(len(variants), len(args.seeds.split(',')), figsize=(14, 5*len(variants)), squeeze=False)
    for ax, (m, s) in zip(axs.ravel(), runs):
        panel(ax, m.snapshots[-1], m.c, f"{s['variant']} | seed {m.c.seed} | {s['status']}\nt={m.time:.3g} s, N(Ag)={m.trace[-1]['metal']}")
    fig.suptitle('Ag / AlN(20) / AlScN(20) / AlN(20) / Pt | 2D effective KMC\nGold: neutral Ag   Teal: Ag+   Dashed: layer interfaces | exploratory, uncalibrated', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .94))
    fig.savefig(out/'morphology.png', dpi=190); fig.savefig(out/'morphology.svg'); plt.close(fig)
    m, s = runs[0]
    ids = np.unique(np.linspace(0, len(m.snapshots)-1, min(6, len(m.snapshots))).astype(int))
    fig, axs = plt.subplots(1, len(ids), figsize=(3.3*len(ids), 5), squeeze=False)
    for ax, i in zip(axs.ravel(), ids):
        panel(ax, m.snapshots[i], m.c, f"t={m.trace[i]['time_s']:.3g} s\nAg={m.trace[i]['metal']}, Ag+={m.trace[i]['ions']}")
    fig.suptitle(f'Natural growth from virgin state | seed {m.c.seed} | no preloaded filament')
    fig.tight_layout(rect=(0, 0, 1, .92)); fig.savefig(out/'evolution.png', dpi=190); fig.savefig(out/'evolution.svg'); plt.close(fig)
    (out/'summary.json').write_text(json.dumps(summaries, indent=2), encoding='utf-8')
    with (out/'summary.csv').open('w', newline='', encoding='utf-8-sig') as f:
        keys = ['variant', 'seed', 'status', 'connected', 'time_s', 'events', 'metal', 'ions', 'ag_balance_error', 'wall_s']
        writer = csv.DictWriter(f, keys); writer.writeheader()
        for s in summaries:
            writer.writerow({k: (s['config']['seed'] if k == 'seed' else s['trace'][-1].get(k, s.get(k))) for k in keys})

if __name__ == '__main__': main()
