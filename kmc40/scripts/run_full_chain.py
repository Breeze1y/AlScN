"""Reproduce the named forming morphology, then run a continuous cycle locally.

FORMING retains the user's archived top-source dynamics. Its electrical
observations are frozen endpoint reads of the common cycle network. Later
phases use the original continuous circuit/electronic/event-clock driver.
"""
from pathlib import Path
from dataclasses import asdict
import argparse, copy, csv, hashlib, json, shutil, time
import numpy as np
import source_cycle_vn40 as source
import cycle_unified40 as cycle
import fit_unified_reset40 as network_table
import form_hourglass40 as geom
from cycle_activity40 import ActivityParams, ActivityState
from mixed_network40 import NetworkParams, prepare, solve
from set_growth40 import GrowthSequence, checkpoint_at_event as growth_checkpoint
from rate_reset40 import ErosionSequence, RateParams, mixed_field_solver, checkpoint_at_event as erosion_checkpoint

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / 'inputs'
cycle.OLD = INPUTS / 'reference_parameters'


def js(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def dump(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    with np.load(path) as z:
        return {k:z[k].copy() for k in z.files}


def save_arrays(folder, arrays):
    folder.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(folder/'states.npz', **arrays)
    keys = [k for k,v in arrays.items() if np.asarray(v).ndim == 1 and len(v) == len(arrays['voltage_V'])]
    with (folder/'iv.csv').open('w', encoding='utf-8-sig', newline='') as f:
        w=csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for i in range(len(arrays['voltage_V'])):
            w.writerow({k:arrays[k][i] for k in keys})


def attach_electronics(folder, state):
    folder=Path(folder); m=js(folder/'restart.json');m['electronic_state']=state
    dump(folder/'restart.json',m);dump(folder/'electronic_state.json',state)
    a=load(folder/'restart.npz');a.update(activity=state['q'],D=state['D'],H=state['H'],elapsed_units=state['elapsed_units'])
    np.savez_compressed(folder/'restart.npz',**a)


def common_parameters():
    pos=js(INPUTS/'cycle_reference/SET_summary.json')['clock']
    neg=js(INPUTS/'cycle_reference/RESET_summary.json')['clock']
    return pos,neg,NetworkParams(**js(INPUTS/'reference_parameters/network_params.json'))


def measured(phase):
    with (INPUTS/'cycle_reference/cycle074.csv').open(encoding='utf-8-sig') as f:
        rows=[r for r in csv.DictReader(f) if r['phase']==phase]
    return np.array([float(r['voltage_V']) for r in rows]),np.array([float(r['current_A']) for r in rows])


def run_forming(out, pos, network):
    folder=out/'FORMING'; summary=folder/'summary.json'
    if summary.exists():
        return folder/'checkpoint'
    reference=INPUTS/'forming_reference'; settings=js(reference/'summary.json')
    model=source.Pilot(settings['seed'],settings['eta'],settings['migration_saddle_reduction_eV'],settings['pt_mode'],settings['VN_fraction'])
    model.form_scale=settings['forming_supply_prefactor']
    print('FORMING: exact top-source replay starts', flush=True)
    model.scan('forming',0,model.form_scale)
    model.save(folder/'legacy')
    arrays=load(folder/'legacy/states.npz'); original=load(reference/'states.npz')
    pick=np.flatnonzero(original['phase'].astype(str)=='forming')
    compared=[]
    for k in ['metal_id','voltage_V','current_A','Ag_total','injected','returned','connected']:
        np.testing.assert_array_equal(arrays[k],original[k][pick]);compared.append(k)
    np.testing.assert_array_equal(arrays['VN'],original['VN'])
    state=ActivityState(params=ActivityParams(**pos['activity_params']));initial=state.to_dict()
    circuit=cycle.CircuitParams(**pos['circuit']); records=[]
    for i,v in enumerate(arrays['voltage_V']):
        g=solve(prepare(arrays['metal_id'][i]>0,arrays['VN']),network,activity=state.q)['conductance_S']
        r=cycle.read_circuit(g,state.q,float(v),circuit)
        state.advance(r['material_voltage_V'],1.,current_A=max(r['material_current_A'],0.))
        records.append(dict(current_A=r['current_A'],endpoint_current_A=r['current_A'],
            read_G_S=r['read_G_S'],terminal_voltage_V=r['terminal_voltage_V'],
            material_voltage_V=r['material_voltage_V'],activity=state.q,D=state.D,H=state.H,
            event_index=sum(e['point']<=i for e in model.events)))
    arrays['legacy_current_A']=arrays['current_A'].copy(); arrays['legacy_read_G_S']=arrays['read_G_S'].copy()
    arrays.update({k:np.array([r[k] for r in records]) for k in records[0]})
    arrays['metal_connected']=arrays['connected'].copy()
    arrays['branch']=np.where(np.arange(len(records))<=100,'outward','return')
    save_arrays(folder,arrays)
    checkpoint=folder/'checkpoint';checkpoint.mkdir(parents=True,exist_ok=True)
    for name in ['restart.npz','restart.json']:
        shutil.copy2(folder/'legacy'/name,checkpoint/name)
    attach_electronics(checkpoint,state.to_dict())
    dump(summary,dict(schema='full_chain40/forming/v1',phase='FORMING',
        source_reference='inputs/forming_reference',source_sha256=sha(reference/'states.npz'),
        original_fields_exact=compared+['VN'],initial_electronic_state=initial,final_electronic_state=state.to_dict(),
        clock=pos,network_params=asdict(network),final_Ag=int(arrays['Ag_total'][-1]),
        final_metal_connected=bool(arrays['metal_connected'][-1]),final_read_G_S=float(arrays['read_G_S'][-1]),
        final_read_ohm=1/float(arrays['read_G_S'][-1]),
        current_semantics='Common-network frozen endpoint current at each archived forming morphology; legacy_current_A retains original readout. No within-dwell event times are asserted for the original forming kernel.',
        forming_dynamics='Original Pilot forming protocol, including its return-sweep oxidation rule; dynamics were not recalibrated',
        following_phase_reason='Metal disconnected but still electrically conductive; a real INIT_RESET creates the subsequent SET high-resistance initial state',
        units='0.5 nm coarse Ag packets; effective sampling units, not calibrated seconds'))
    print('FORMING complete',dict(Ag=int(arrays['Ag_total'][-1]),R=1/arrays['read_G_S'][-1]),flush=True)
    return checkpoint


def new_sequence(kind, checkpoint, network):
    if kind=='growth':
        return GrowthSequence.from_checkpoint(checkpoint)
    m=js(checkpoint/'restart.json');a=load(checkpoint/'restart.npz')
    provenance=dict(path=str(checkpoint.relative_to(ROOT)) if checkpoint.is_relative_to(ROOT) else str(checkpoint),
                    sha256=sha(checkpoint/'restart.npz'),injected=m['injected'],returned=m['returned'],electronic_state=m['electronic_state'])
    return ErosionSequence(a['metal_id'],a['VN'],provenance,RateParams(),m.get('seed',314),
        mixed_field_solver(network,1.),rng_state=m['rng_state'])


def run_phase(out,name,kind,checkpoint,params,network,substeps=16,phase_factory=None):
    folder=out/name
    if (folder/'summary.json').exists():return folder/'checkpoint'
    cat=out/'catalogues'/name; initial_meta=js(checkpoint/'restart.json')
    initial=ActivityState.from_dict(initial_meta['electronic_state'])
    v,reference=measured('SET' if kind=='growth' else 'RESET')
    if (cat/'sequence.npz').exists():
        sequence=GrowthSequence.restore(cat) if kind=='growth' else ErosionSequence.restore(cat)
    else:
        print(name,'generate 400 spatial events',flush=True)
        sequence=new_sequence(kind,checkpoint,network).extend(400);sequence.save(cat)
    while True:
        phase=(phase_factory or cycle.Phase)(cat)
        rr,end,info=cycle.simulate(phase,v,initial,params['rate_scale'],params['voltage_scale_V'],
                                  cycle.CircuitParams(**params['circuit']),substeps=substeps)
        if not info['sequence_exhausted']:break
        print(name,'extending catalogue',len(sequence.events),flush=True)
        if len(sequence.events)>3000:raise RuntimeError('Unexpected event budget exhaustion')
        sequence.extend(300).save(cat)
    saved=cycle.save_phase(folder,phase,rr,initial,end,info,params,reference)
    a=load(folder/'states.npz');a['metal_connected']=np.array([bool(geom.span(m>0).any()) for m in a['metal_id']])
    a['Ag_bottom10nm']=(a['metal_id'][:,-20:]>0).sum((1,2))
    a['Ag_bottom_interface']=(a['metal_id'][:,-1]>0).sum(1)
    a['branch']=np.where(np.arange(len(v))<=100,'outward','return')
    save_arrays(folder,a)
    cp=folder/'checkpoint'; fn=growth_checkpoint if kind=='growth' else erosion_checkpoint
    kwargs={'observation_exposure':info['event_exposure']}
    if kind=='growth':kwargs['electronic_state']=end.to_dict()
    fn(cat,info['final_event'],cp,**kwargs);attach_electronics(cp,end.to_dict())
    boundary=load(checkpoint/'restart.npz')
    np.testing.assert_array_equal(a['metal_id'][0],boundary['metal_id'])
    np.testing.assert_array_equal(a['VN'],boundary['VN'])
    assert int(a['injected'][0])==initial_meta['injected'] and int(a['returned'][0])==initial_meta['returned']
    saved.update(schema='full_chain40/phase/v1',phase=name,
        source_checkpoint=str(checkpoint.relative_to(ROOT)) if checkpoint.is_relative_to(ROOT) else str(checkpoint),
        source_checkpoint_sha256=sha(checkpoint/'restart.npz'),network_params=asdict(network),
        final_Ag=int(a['Ag_total'][-1]),final_metal_connected=bool(a['metal_connected'][-1]),
        final_read_G_S=float(a['read_G_S'][-1]),final_read_ohm=1/float(a['read_G_S'][-1]),
        current_semantics='Within-dwell time averaged current; morphology and endpoint_current_A are sampled endpoint states',
        measured_role='C74 mature cycle reference; new forming-history trajectory, not the original calibrated geometry',
        protocol_substeps=substeps,parameter_role='Frozen cycle_unified40 selected parameters',
        inherited_rng_state=initial_meta['rng_state'])
    dump(folder/'summary.json',saved)
    print(name,'complete',dict(events=info['final_event'],Ag=int(a['Ag_total'][-1]),R=saved['final_read_ohm'],RMSE=saved['metrics']['all']),flush=True)
    return cp


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path)
    ap.add_argument('--forming-only',action='store_true');ap.add_argument('--substeps',type=int,default=16)
    ap.add_argument('--refine-from',type=Path,help='Reuse exact forming and matching spatial catalogues, refine only numerical q interpolation')
    args=ap.parse_args()
    if args.out is None:
        selected=ROOT/'results/selected_run.json'
        args.out=ROOT/js(selected)['run_directory'] if selected.exists() else ROOT/'results/full_chain'
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    if args.refine_from:
        grid=np.unique(np.r_[np.linspace(-10.,0.,21),np.arange(-5.75,0.,.5)])
    elif (out/'parameters.json').exists():
        grid=np.array(js(out/'parameters.json').get('logq_grid',np.linspace(-10.,0.,21)),float)
    elif (ROOT/'results/selected_run.json').exists():
        selected_parameters=ROOT/js(ROOT/'results/selected_run.json')['run_directory']/'parameters.json'
        grid=np.array(js(selected_parameters)['logq_grid'],float)
    else:grid=np.linspace(-10.,0.,21)
    network_table.LOGQ=grid;cycle.LOGQ=grid
    pos,neg,network=common_parameters()
    dump(out/'parameters.json',dict(schema='full_chain40/parameters/v1',positive=pos,negative=neg,network=asdict(network),
        source_model='cycle_unified40 joint_simple',forming_model='vacancy_assisted40 assisted',
        phases=['FORMING','INIT_RESET','SET','RESET'],initial_electronic_state='D=H=0, q=1 at pristine forming; continuous thereafter',
        logq_grid=grid.tolist(),q_grid_reason='Midpoint refinement across the observed log10(q) range -6 to 0; physics and clocks unchanged' if len(grid)>21 else 'Original 21-node network lookup',
        units='effective sampling units; Ag packets are coarse volumes'))
    if args.refine_from and not (out/'FORMING').exists():
        shutil.copytree(args.refine_from/'FORMING',out/'FORMING')
    cp=run_forming(out,pos,network)
    if args.forming_only:return
    for name,kind,p in [('INIT_RESET','erosion',neg),('SET','growth',pos),('RESET','erosion',neg)]:
        cat=out/'catalogues'/name
        if args.refine_from and not cat.exists():
            oldcat=args.refine_from/'catalogues'/name
            oldmeta=js(oldcat/'restart.json');newmeta=js(cp/'restart.json')
            oldstate=load(oldcat/'sequence.npz');newstate=load(cp/'restart.npz')
            if (np.array_equal(oldstate['metal_id'][0],newstate['metal_id']) and
                np.array_equal(oldstate['VN'],newstate['VN']) and
                oldmeta['initial_rng_state']==newmeta['rng_state']):
                shutil.copytree(oldcat,cat)
                meta=js(cat/'restart.json')
                meta['initial_source']=dict(meta['initial_source'],path=str(cp),
                    injected=newmeta['injected'],returned=newmeta['returned'],
                    electronic_state=newmeta['electronic_state'])
                if kind=='growth':meta['electronic_state']=newmeta['electronic_state']
                dump(cat/'restart.json',meta)
                print(name,'reusing exact spatial catalogue; recomputing q table with',len(grid),'nodes',flush=True)
        cp=run_phase(out,name,kind,cp,p,network,args.substeps)
    summaries={name:js(out/name/'summary.json') for name in ['FORMING','INIT_RESET','SET','RESET']}
    dump(out/'chain_summary.json',dict(schema='full_chain40/chain/v1',complete=True,phases=summaries,
        state_chain='Pristine -> reproduced forming -> electrical initialization RESET -> SET -> RESET',
        no_manual_geometry_or_q_resets=True))
    print('CHAIN COMPLETE',out,flush=True)


if __name__=='__main__':main()
