"""Controlled cohesion assays and fresh-forming validation; original kernel unchanged."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse,copy,hashlib,json,sys,time
from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from kmc_model import Config
from forming_iv import ElectronicConfig
from dynamic_kmc import DynamicCachedModel
from morphology_metrics import metrics
from run_kmc import hashes,write_json
BASE=json.loads((OUT/'baseline_settings.json').read_text(encoding='utf-8'))

def assay(bond,seed):
    folder=OUT/'assays'/f'J{bond:.3f}_s{seed}';folder.mkdir(parents=True,exist_ok=True)
    source=ROOT/BASE['restore_reference']
    with np.load(source,allow_pickle=False) as ar:
        state=ar['state'].copy();meta=json.loads(str(ar['metadata']))
    settings=copy.deepcopy(BASE);settings['model'].update(bond_eV=bond,seed=seed,voltage=0.)
    m=DynamicCachedModel(Config(**settings['model']),ElectronicConfig(**settings['electronics']))
    m.state=state.ravel().copy();m.steps=meta['forming_events'];m.counts=meta['counts'].copy();m.net_source=meta['net_source']
    m.time=0.;m.command_voltage=1.35;m.first_bridge=dict(time_s=0.,events=m.steps,command_voltage_V=1.35)
    m.update_structure();m.topology_log=[];m.iv=[];m.trace=[];m.snapshots=[];m.record('assay_initial');m.save()
    initial_steps=m.steps;initial_counts=m.counts.copy();started=time.perf_counter()
    m.advance(.1,m.steps+1000000);m.record('assay_end')
    on=connected=0.
    for a,b in zip(m.iv,m.iv[1:]):
        dt=b['time_s']-a['time_s'];on+=dt*(abs(a['current_A'])>=.0005);connected+=dt*a['connected']
    summary=dict(kind='shared_geometry_positive_bias_assay_not_fresh_forming',bond_eV=bond,seed=seed,
        command_voltage_V=1.35,duration_s=m.time,connected_time_fraction=connected/m.time,
        on_time_fraction=on/m.time,ruptures=sum(not v['connected'] for v in m.topology_log),
        final_connected=m.connected,events=m.steps-initial_steps,counts={k:m.counts[k]-initial_counts[k] for k in m.counts},
        ag_balance_error=int(np.count_nonzero(m.state[m.interior]))-m.net_source,
        reaction_balance_error=int(np.sum(m.state[m.interior]==2))-(m.counts['reduce']-m.counts['oxidize']),
        initial_state_sha256=hashlib.sha256(state.tobytes()).hexdigest(),initial_source=str(source),
        initial_rng='New independent event RNG per seed; not an exact continuation of the old trajectory',
        endpoint=m.iv[-1],morphology=metrics(m.state.reshape(m.shape),1,20),wall_s=time.perf_counter()-started,
        config=asdict(m.c),physics_hashes=hashes())
    assert summary['ag_balance_error']==summary['reaction_balance_error']==0
    write_json(folder/'result.json',summary)
    np.savez_compressed(folder/'states.npz',initial=state,state=m.state.reshape(m.shape))
    print('ASSAY',bond,seed,'on_fraction',summary['on_time_fraction'],'ruptures',summary['ruptures'],flush=True)
    return summary

def scout(bond,nu,seed):
    sys.path.insert(0,str(ROOT/'results/opening_calibration_20261113'))
    import calibration_driver as runner
    runner.BASE=copy.deepcopy(BASE);runner.BASE['model']['bond_eV']=bond
    runner.OUT=OUT/f'bond_{bond:.3f}';runner.OUT.mkdir(exist_ok=True)
    result=runner.trial(nu,seed)
    folder=Path(result['folder']);audit=json.loads((folder/'manifest.json').read_text())
    audit.update(experiment_driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        experiment='Stronger effective Ag-Ag cohesion; original redox and migration energy formulas retained',bond_eV=bond)
    write_json(folder/'manifest.json',audit)
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['assay','scout'])
    parser.add_argument('--bonds',type=float,nargs='+',default=[.06,.08,.10]);parser.add_argument('--nus',type=float,nargs='+',default=[2.42e12])
    parser.add_argument('--seeds',type=int,nargs='+',default=[7,17,27]);parser.add_argument('--workers',type=int,default=6)
    args=parser.parse_args()
    jobs=[(b,s) for b in args.bonds for s in args.seeds] if args.mode=='assay' else [(b,n,s) for b in args.bonds for n in args.nus for s in args.seeds]
    with ProcessPoolExecutor(max_workers=min(args.workers,len(jobs))) as pool:
        futures=[pool.submit(assay if args.mode=='assay' else scout,*j) for j in jobs]
        summaries=[f.result() for f in as_completed(futures)]
    write_json(OUT/(args.mode+'_round_summary.json'),summaries)
    print('ROUND_COMPLETE',args.mode,flush=True)

if __name__=='__main__':main()
