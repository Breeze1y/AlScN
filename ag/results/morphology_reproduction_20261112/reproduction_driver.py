"""Audit driver: recompute the six Nov 8 morphologies from blank lattices."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import sys
import time
import zipfile
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from kmc_model import Config
from forming_iv import ElectronicConfig,run_forming
from fast_kmc import CachedFormingModel
from morphology_metrics import metrics

OUT=Path(__file__).resolve().parent
REFERENCE=ROOT/'results/vacancy_forming_20261108'

def case(name):
    reference=json.loads((REFERENCE/f'{name}.json').read_text(encoding='utf-8'))
    cfg=Config(**reference['config']);electronics=ElectronicConfig(**reference['electronics'])
    started=time.perf_counter()
    def progress(i,m,status):
        if i>=31:print(name,'segment',i,'events',m.steps,'time',m.time,'connected',m.connected,flush=True)
    m,summary,states,records=run_forming(cfg,electronics,reference['protocol'],on_segment=progress,model_class=CachedFormingModel)
    state=m.state.reshape(m.shape)
    with np.load(REFERENCE/f'{name}.npz',allow_pickle=False) as old:
        exact=np.array_equal(state,old['scan_states'][-1])
        maps_equal=np.array_equal(m.site.reshape(m.shape),old['site_energy_eV']) and np.array_equal(m.vacancy_mask.reshape(m.shape),old['vacancy_mask'])
    checks=dict(fresh_blank_initialization=True,final_state_exact=exact,static_maps_exact=maps_equal,
        events_exact=m.steps==reference['events'],counts_exact=m.counts==reference['counts'],
        bridge_time_exact=m.first_bridge['time_s']==reference['first_bridge']['time_s'],
        mass_balanced=int(np.count_nonzero(state[1:-1]))==m.net_source==m.counts['inject']-m.counts['absorb'])
    summary.update(name=name,wall_s=time.perf_counter()-started,checks=checks,
                   morphology=metrics(state,cfg.spacing_nm,cfg.layer_cells),fresh_from_blank=True)
    (OUT/f'{name}.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    np.savez_compressed(OUT/f'{name}.npz',state=state,scan_states=states,vacancy_mask=m.vacancy_mask.reshape(m.shape),site_energy_eV=m.site.reshape(m.shape))
    print('FINISHED',name,'exact',exact,'checks',all(checks.values()),'wall_s',summary['wall_s'],flush=True)
    return summary

def main():
    oldmanifest=json.loads((REFERENCE/'manifest.json').read_text())
    with zipfile.ZipFile(REFERENCE/'source_snapshot.zip') as z:
        assert all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in oldmanifest.items())
    names=[f'{variant}_seed_{seed}' for variant in ('uncoupled','vacancy_penalty') for seed in (7,17,27)]
    source_names=['kmc_model.py','forming_iv.py','fast_kmc.py','morphology_metrics.py']
    audit=dict(reference=str(REFERENCE),reference_source_archive_verified=True,
        reference_image_sha256=hashlib.sha256((REFERENCE/'paired_morphology.png').read_bytes()).hexdigest(),
        current_sources_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in source_names},
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        mode='fresh simulations; saved arrays used only for post-run comparison',
        scope='Historical 3 V fast ramp and first-bridge stop; not 50 Hz full-loop validation')
    (OUT/'manifest.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    with ProcessPoolExecutor(max_workers=3) as pool:
        jobs={pool.submit(case,n):n for n in names}
        found={jobs[j]:j.result() for j in as_completed(jobs)}
    summaries=[found[n] for n in names]
    (OUT/'summary.json').write_text(json.dumps(summaries,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    for ax,name,summary in zip(axes.ravel(),names,summaries):
        with np.load(OUT/f'{name}.npz') as a:s=a['state'];vacancy=a['vacancy_mask']
        ax.imshow(s,cmap=ListedColormap(['#edf2fa','#2eb1a5','#d89319']),vmin=0,vmax=2,
                  extent=(0,150,61,-1),interpolation='none')
        for d in (0,60):ax.axhspan(d-1 if d==0 else d,d if d==0 else d+1,color='#8797a5')
        for d in (20,40):ax.axhline(d,color='#596586',ls='--',lw=.7)
        z,x=np.nonzero(vacancy[1:-1]&(s[1:-1]==0));ax.scatter(x+.5,z+.5,s=1.5,c='#7b8598',alpha=.75,linewidths=0)
        ax.set(xlabel='x (nm)',ylabel='Depth from Ag (nm)',title=f'{name}\nfirst bridge t={summary["morphology_time_s"]:.3f} s; exact={summary["checks"]["final_state_exact"]}')
    fig.suptitle('Fresh reproduction of Nov 8 reference | gold: metal Ag; teal: Ag+\n3 V fast-ramp first-bridge states; not zero-bias end states')
    fig.savefig(OUT/'reproduced_paired_morphology.png',dpi=175);plt.close(fig)
    ok=all(all(s['checks'].values()) for s in summaries)
    (OUT/'verification.json').write_text(json.dumps(dict(all_six_exact=ok,cases={s['name']:s['checks'] for s in summaries}),indent=2),encoding='utf-8')
    print('ALL_SIX_EXACT',ok,flush=True)

if __name__=='__main__':main()
