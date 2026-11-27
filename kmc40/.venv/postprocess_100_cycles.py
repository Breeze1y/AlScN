import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPLAY = ROOT / 'replay/cycles100'
CONTINUOUS = ROOT / 'replay/reset68/continuous'
REFERENCE = ROOT / 'results/cycle_reset_fit40/mid_clock_set'
PHASES = ('SET', 'RESET')
FIELDS = ('point', 'branch', 'voltage_V', 'current_A', 'endpoint_current_A')
POINTS = 202
CELL_NM = 0.5
COLORS = ['#ffffff', '#ced5dd', '#b36832', '#146a80']
SET_COLOR, RESET_COLOR = '#146a80', '#a34e2a'

sys.path.insert(0, str(ROOT / 'scripts'))
from form_hourglass40 import span

mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                     'font.size': 9, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                     'axes.labelcolor': '#222b33', 'text.color': '#222b33'})


def collect(cycles):
    data = []
    prior = np.load(REFERENCE / 'INIT_RESET/checkpoint/restart.npz')
    with prior:
        prior_metal = prior['metal_id']
        fixed_vn = prior['VN']
    with (CONTINUOUS / 'iv/iv_all.csv').open(encoding='utf-8-sig', newline='') as handle:
        merged = csv.DictReader(handle)
        assert set(('cycle', 'phase', *FIELDS)).issubset(merged.fieldnames)
        for cycle in range(1, cycles + 1):
            record = {}
            for phase in PHASES:
                folder = CONTINUOUS / f'cycle_{cycle:02d}' / phase
                with (folder / 'iv.csv').open(encoding='utf-8-sig', newline='') as individual:
                    source_rows = list(csv.DictReader(individual))
                assert len(source_rows) == POINTS, folder
                for point, raw in enumerate(source_rows):
                    row = next(merged)
                    assert (int(row['cycle']), row['phase'], int(row['point'])) == (cycle, phase, point), row
                    assert row['branch'] == ('outward' if point <= 100 else 'return'), row
                    assert all(row[field] == raw[field] for field in FIELDS), row
                summary = json.loads((folder / 'summary.json').read_text(encoding='utf-8-sig'))
                with np.load(folder / 'states.npz') as state, np.load(folder / 'checkpoint/restart.npz') as checkpoint:
                    voltage = state['voltage_V']
                    current = state['current_A']
                    for field in FIELDS[2:]:
                        np.testing.assert_array_equal(state[field], [float(row[field]) for row in source_rows])
                    assert np.isfinite(current).all() and np.isfinite(state['endpoint_current_A']).all(), folder
                    assert voltage[0] == voltage[-1] == 0 and voltage[100] == (1.5 if phase == 'SET' else -1.5), folder
                    metal = state['metal_id']
                    assert metal.shape == (POINTS, 80, 100), folder
                    np.testing.assert_array_equal(state['VN'], fixed_vn)
                    np.testing.assert_array_equal(metal[0], prior_metal)
                    np.testing.assert_array_equal(metal[-1], checkpoint['metal_id'])
                    prior_metal = metal[-1].copy()
                    np.testing.assert_allclose(state['Ag_total'] + state['returned'], state['injected'], rtol=0, atol=1e-12)
                    connected = state['metal_connected']
                    assert bool(span(prior_metal > 0).any()) == connected[-1] == summary['final_metal_connected'], folder
                    resistance = float(1. / state['read_G_S'][-1])
                    assert np.isclose(resistance, summary['final_read_ohm'], rtol=1e-12), folder
                    changes = np.flatnonzero((~connected[:-1] & connected[1:]) if phase == 'SET' else
                                             (connected[:-1] & ~connected[1:])) + 1
                    record[phase] = dict(voltage=voltage, current=current, connected=bool(connected[-1]),
                                         resistance=resistance, ag=float(state['Ag_total'][-1]),
                                         switch_v=float(voltage[changes[0]]) if len(changes) else None,
                                         end=prior_metal > 0, delayed=summary.get('delayed_return_actions', 0),
                                         elapsed=float(summary.get('elapsed_seconds', 0)))
            data.append(record)
        if cycles == 100:
            assert next(merged, None) is None, 'Unexpected rows after cycle 100'
    return data, fixed_vn


def draw_iv(out, data):
    cycle = np.arange(1, len(data) + 1)
    fig, (ax, rs) = plt.subplots(1, 2, figsize=(9.4, 4.15), width_ratios=[1.4, 1], layout='constrained')
    for record in data:
        for phase in PHASES:
            trace = record[phase]
            current = np.abs(trace['current'])
            ax.plot(trace['voltage'], np.ma.masked_equal(current, 0), color='#8896a3', alpha=0.20, lw=0.7)
    for phase, color in zip(PHASES, (SET_COLOR, RESET_COLOR)):
        trace = data[-1][phase]
        for branch, linestyle, label in [(slice(0, 101), '-', f'Cycle {len(data)} {phase}'), (slice(100, None), '--', None)]:
            amps = np.ma.masked_equal(np.abs(trace['current'][branch]), 0)
            ax.plot(trace['voltage'][branch], amps, color=color, ls=linestyle, lw=1.5, label=label)
        resistances = [record[phase]['resistance'] for record in data]
        rs.plot(cycle, resistances, color=color, lw=1.4, label=f'{phase} endpoint')
    for axis in (ax, rs):
        axis.set_yscale('log')
        axis.grid(axis='y', which='major', color='#e5e9ed', lw=0.7)
        axis.legend(loc='upper center', bbox_to_anchor=(0.5, -0.19), ncols=2, fontsize=8)
    ax.set(xlabel='Programmed voltage (V)', ylabel='Dwell-averaged |current| (A)',
           title=f'I–V overlay across {len(data)} cycles')
    rs.set(xlabel='Cycle', ylabel='Endpoint read resistance (Ω)', title='Resistance after SET and RESET')
    fig.savefig(out / 'iv_100cycles.png', dpi=300)
    fig.savefig(out / 'iv_100cycles.pdf')
    fig.savefig(out / 'iv_100cycles.svg')
    plt.close(fig)


def morphology_page(data, vn, groups, title):
    fig = plt.figure(figsize=(10.0, 1.35 + 2.8 * len(groups)))
    grid = fig.add_gridspec(2*len(groups), 5, left=0.104, right=0.97, top=0.89, bottom=0.145,
                            wspace=0.13, hspace=0.33 if len(groups) == 1 else 0.26)
    fig.suptitle(title, y=0.97, fontweight='bold', fontsize=13)
    cmap = ListedColormap(COLORS)
    for group_number, group in enumerate(groups):
        for row, phase in enumerate(PHASES):
            for col, cycle in enumerate(group):
                ax = fig.add_subplot(grid[2 * group_number + row, col])
                layer = data[cycle-1][phase]['end']
                pixels = np.where(span(layer), 3, np.where(layer, 2, vn.astype(np.uint8)))
                ax.imshow(pixels, cmap=cmap, vmin=0, vmax=3, interpolation='nearest',
                          extent=(0, layer.shape[1]*CELL_NM, layer.shape[0]*CELL_NM, 0))
                ax.set(xticks=[], yticks=[])
                for spine in ax.spines.values():
                    spine.set_visible(True)
                    spine.set_linewidth(0.6)
                    spine.set_color('#a9b4be')
                if not row:
                    ax.set_title(f'Cycle {cycle}', fontsize=9, pad=4)
                if not col:
                    ax.set_ylabel(phase, rotation=0, rotation_mode='anchor', ha='right', va='center', labelpad=7, fontweight='bold')
    handles = [Patch(facecolor=color, edgecolor='#9aa3ad', label=label)
               for color, label in zip(COLORS, ('Host', 'N vacancy', 'Other Ag', 'Spanning Ag'))]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.053), ncols=4, fontsize=9)
    fig.text(0.5, 0.033, 'Field of view 50 × 40 nm per tile. Top electrode at upper edge. Grid spacing 0.5 nm.',
             ha='center', va='bottom', fontsize=8)
    return fig


def draw_morphology(out, data, vn):
    pdf, preview = out / 'morphology_100cycles.pdf', out / 'morphology_overview.png'
    milestones = sorted(set(np.rint(np.linspace(1, len(data), 5)).astype(int)))
    with PdfPages(pdf) as book:
        figure = morphology_page(data, vn, [milestones], 'Ag morphology at selected cycle endpoints')
        book.savefig(figure, dpi=300)
        figure.savefig(preview, dpi=300)
        plt.close(figure)
        for start in range(1, len(data)+1, 10):
            end = min(start+9, len(data))
            groups = [range(group, min(group+5, end+1)) for group in range(start, end+1, 5)]
            figure = morphology_page(data, vn, groups, f'Ag morphology at cycle endpoints {start}–{end}')
            book.savefig(figure, dpi=300)
            plt.close(figure)


def describe(data):
    def triplet(values):
        return [float(np.min(values)), float(np.median(values)), float(np.max(values))] if len(values) else []
    for phase in PHASES:
        items = [record[phase] for record in data]
        bridges = [n for n, item in enumerate(items, 1) if item['connected']]
        print(phase, dict(bridged=len(bridges), exception_cycles=(sorted(set(range(1, len(data)+1))-set(bridges)) if phase == 'SET' else bridges),
                          R_ohm=triplet([a['resistance'] for a in items]),
                          Ag_packets=triplet([a['ag'] for a in items]),
                          first_spanning_change_V=triplet([a['switch_v'] for a in items if a['switch_v'] is not None]),
                          waited_cycles={n: a['delayed'] for n, a in enumerate(items, 1) if a['delayed']},
                          replay_stage_seconds=sum(a['elapsed'] for a in items[5:])))
        if len(items) == 100:
            for endpoints in (items[:10], items[-10:]):
                print(phase, '10-cycle medians', {name: float(np.median([a[name] for a in endpoints]))
                                                     for name in ('ag', 'resistance')})
    print('R_RESET/R_SET min / median / max', triplet([r['RESET']['resistance']/r['SET']['resistance'] for r in data]))
    for n in sorted({1, min(5, len(data)), min(16, len(data)), round(len(data)/2), len(data)}):
        print('Cycle', n, {phase: {k: data[n-1][phase][k] for k in ('ag', 'resistance', 'connected', 'switch_v')} for phase in PHASES})
    print(f'Verified {len(data)} cycles and {len(data)*2*POINTS} IV rows against individual IV and NPZ outputs')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cycles', type=int, default=100)
    parser.add_argument('--output', type=Path, default=REPLAY / 'report')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    results, vacancies = collect(args.cycles)
    describe(results)
    if args.cycles == 100:
        from plot_continuous_iv import draw, read_curves
        from reset68_morphology_atlas import main as draw_atlas
        draw(read_curves(), args.output)
        draw_atlas(args.output)
    else:
        draw_iv(args.output, results)
        draw_morphology(args.output, results, vacancies)
