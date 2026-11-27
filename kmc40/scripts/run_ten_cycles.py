"""Extend the saved local-retreat experiment to a continuous ten-cycle chain.

The first cycle is reused exactly. Later spatial catalogues are newly generated
from the preceding physical checkpoint and RNG. Frozen circuit and kinetics
are evaluated on the same experimental voltage protocol, without per-cycle fit.
"""
from pathlib import Path
from dataclasses import asdict
import argparse,csv,hashlib,json,math,shutil,time
import numpy as np
import cycle_unified40 as cycle
import run_full_chain as full
import run_local_retreat as local
import multi_cycle_retreat as multi
import local_fraction_network as fraction_net
from mixed_network40 import prepare,solve
from cycle_activity40 import ActivityState
from set_growth40 import GrowthSequence,checkpoint_at_event as growth_checkpoint
from rate_reset40 import ErosionSequence,checkpoint_at_event as erosion_checkpoint
from form_hourglass40 import span

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/full_chain_refined'
FIRST_RESET=ROOT/'results/local_retreat40/local_tau8_p12'
OUT=ROOT/'results/ten_cycles_local_retreat40'
load,js,dump,sha=full.load,full.js,full.dump,full.sha

class LazyNetwork:
    """The old log-q interpolation with its grid nodes evaluated only on demand."""
    def __init__(self,ids,vn,params,grid,folder):
        self.ids=ids;self.vn=vn;self.params=params;self.grid=np.asarray(grid,float)
        self.folder=Path(folder);self.cache={};self.topology={};self.calls=0;self.keys=[]
        self.identity=sha_arrays(vn,self.grid,extra=json.dumps(asdict(params),sort_keys=True))
        for m in ids:self.keys.append(hashlib.sha256((m>0).tobytes()).hexdigest())
        path=self.folder/'network_activity_lazy.npz'
        if path.exists():
            z=load(path)
            if str(z['identity'])!=self.identity:raise ValueError('Network cache parameters changed')
            for k,row in zip(z['geometry_sha256'].astype(str),z['logG']):self.cache[k]=row.copy()
    def g(self,index,q):
        key=self.keys[index]
        row=self.cache.setdefault(key,np.full(len(self.grid),np.nan))
        x=math.log10(max(q,1e-30));x=float(np.clip(x,self.grid[0],self.grid[-1]))
        hi=int(np.searchsorted(self.grid,x,side='left'))
        nodes=[hi] if hi==0 or self.grid[hi]==x else [hi-1,hi]
        if np.any(~np.isfinite(row[nodes])):
            if key not in self.topology:self.topology[key]=prepare(self.ids[index]>0,self.vn)
            for j in nodes:
                if not np.isfinite(row[j]):
                    g=solve(self.topology[key],self.params,activity=10.**self.grid[j])['conductance_S']
                    row[j]=math.log(max(g,1e-30));self.calls+=1
        return math.exp(float(np.interp(x,self.grid,row)))
    def save(self):
        self.folder.mkdir(parents=True,exist_ok=True)
        rows=np.array([self.cache.get(k,np.full(len(self.grid),np.nan)) for k in self.keys])
        np.savez_compressed(self.folder/'network_activity_lazy.npz',logG=rows,log_activity=self.grid,
            geometry_sha256=np.array(self.keys),identity=self.identity,params=json.dumps(asdict(self.params)))
        return dict(grid_nodes=len(self.grid),evaluated_unique_cells=sum(np.isfinite(r).sum().item() for r in self.cache.values()),
            geometries=len(self.keys),unique_geometries=len(set(self.keys)),new_solves=self.calls,
            method='Same 33-node log-q interpolation as the frozen chain; only queried bracketing grid entries are computed.')

def sha_arrays(*arrays,extra=''):
    h=hashlib.sha256()
    for a in arrays:h.update(np.asarray(a).tobytes())
    h.update(extra.encode());return h.hexdigest()

class ContinuousPhase:
    def __init__(self,folder,network,grid):
        self.folder=Path(folder);z=load(self.folder/'sequence.npz')
        self.ids=z['metal_id'];self.vn=z['VN'];self.clock=z['cumulative_hazard_clock'];self.params=network
        self.injected=z['injected'];self.returned=z['returned']
        self.lookup=LazyNetwork(self.ids,self.vn,network,grid,self.folder)
    def g(self,index,q):return self.lookup.g(index,q)

class ResetCatalogue:
    def __init__(self,folder,checkpoint,negative,network,grid):
        self.folder=Path(folder);z=load(self.folder/'sequence.npz')
        self.ids=z['metal_id'];self.vn=z['VN'];self.clock=z['cumulative_hazard_clock']
        self.injected=int(z['initial_injected']);self.returned=int(z['initial_returned'])
        with (self.folder/'events.csv').open(encoding='utf-8',newline='') as f:self.events=list(csv.DictReader(f))
        routes=load(self.folder/'transport.npz');self.routes=[routes['nodes'][a:b] for a,b in zip(routes['offsets'][:-1],routes['offsets'][1:])]
        v,i=full.measured('RESET');self.base=dict(voltage_V=v,measured_A=i,branch=np.where(np.arange(len(v))<=100,'outward','return'))
        self.summary=dict(clock=negative,initial_electronic_state=js(checkpoint/'restart.json')['electronic_state'])
        self.network=network;self.topology=fraction_net.prepare(self.ids[0]>0,self.vn)
        self.lookup=LazyNetwork(self.ids,self.vn,network,grid,self.folder)
    def binary_g(self,index,q):return self.lookup.g(index,q)

def check_initial(sequence,checkpoint):
    a=load(checkpoint/'restart.npz');m=js(checkpoint/'restart.json')
    np.testing.assert_array_equal(sequence.states[0],a['metal_id']);np.testing.assert_array_equal(sequence.vn,a['VN'])
    assert sequence.initial_rng_state==m['rng_state']
    assert np.count_nonzero(a['metal_id'])==m['injected']-m['returned']
    if 'Ag_units' in a:assert np.all((a['Ag_units']==0)|(a['Ag_units']==int(a['subunits'])))

def prepare_sequence(folder,kind,checkpoint,network):
    if (folder/'sequence.npz').exists():
        sequence=GrowthSequence.restore(folder) if kind=='growth' else ErosionSequence.restore(folder)
    else:
        sequence=full.new_sequence(kind,checkpoint,network)
        sequence.extend(96 if kind=='growth' else 320);sequence.save(folder)
    check_initial(sequence,checkpoint);return sequence

def source_info(checkpoint,folder,number,name):
    return dict(schema='ten_cycles40/phase/v1',cycle_number=number,phase=name,
        source_checkpoint=str(checkpoint.relative_to(ROOT)),source_checkpoint_sha256=sha(checkpoint/'restart.npz'),
        source_catalogue=str(folder.relative_to(ROOT)),inherited_rng_state=js(checkpoint/'restart.json')['rng_state'],
        measured_role='Same C74 reference protocol for all cycles; only cycle 1 inherited calibration, later cycles continuous frozen-parameter evolution',
        parameter_role='Frozen from full_chain_refined and local_retreat40/local_tau8_p12; no cycle-specific fit')

def run_set(root,number,checkpoint,positive,network,grid,steps):
    target=root/f'cycle_{number:02d}'/'SET';folder=target.parent/'catalogues/SET'
    if (target/'summary.json').exists() and (target/'checkpoint/restart.npz').exists():return target/'checkpoint'
    start=time.perf_counter();sequence=prepare_sequence(folder,'growth',checkpoint,network)
    initial=ActivityState.from_dict(js(checkpoint/'restart.json')['electronic_state']);v,reference=full.measured('SET')
    while True:
        phase=ContinuousPhase(folder,network,grid)
        records,end,info=cycle.simulate(phase,v,initial,positive['rate_scale'],positive['voltage_scale_V'],
            cycle.CircuitParams(**positive['circuit']),substeps=steps)
        cache_info=phase.lookup.save()
        if not info['sequence_exhausted']:break
        print(f'cycle {number} SET extends catalogue from {len(sequence.events)}',flush=True)
        sequence.extend(128).save(folder)
    s=cycle.save_phase(target,phase,records,initial,end,info,positive,reference)
    a=load(target/'states.npz');a['metal_connected']=np.array([bool(span(m>0).any()) for m in a['metal_id']])
    a['branch']=np.where(np.arange(len(v))<=100,'outward','return');a['global_elapsed_units']=initial.elapsed_units+np.arange(1,len(v)+1)
    full.save_arrays(target,a)
    cp=target/'checkpoint';growth_checkpoint(folder,info['final_event'],cp,info['event_exposure'],end.to_dict());full.attach_electronics(cp,end.to_dict())
    s.update(source_info(checkpoint,folder,number,'SET'),network_params=asdict(network),integration_substeps=steps,
        final_Ag=int(a['Ag_total'][-1]),final_read_ohm=float(1/a['read_G_S'][-1]),final_metal_connected=bool(a['metal_connected'][-1]),
        lazy_network=cache_info,elapsed_seconds=time.perf_counter()-start)
    dump(target/'summary.json',s);print(f'cycle {number} SET complete Ag={s["final_Ag"]} R={s["final_read_ohm"]:.6g}',flush=True)
    return cp

def run_reset(root,number,checkpoint,negative,network,grid,retreat,steps):
    target=root/f'cycle_{number:02d}'/'RESET';folder=target.parent/'catalogues/RESET'
    if (target/'summary.json').exists() and (target/'checkpoint/restart.npz').exists():return target/'checkpoint'
    start=time.perf_counter();sequence=prepare_sequence(folder,'erosion',checkpoint,network)
    while True:
        cat=ResetCatalogue(folder,checkpoint,negative,network,grid)
        try:a,s,actions,snaps=multi.simulate(cat,retreat,steps)
        except RuntimeError as exc:
            cat.lookup.save()
            if str(exc)!='Spatial catalogue exhausted':raise
            print(f'cycle {number} RESET extends catalogue from {len(sequence.events)}',flush=True)
            sequence.extend(160).save(folder);continue
        cache_info=cat.lookup.save();break
    assert not s['pending_microsteps'],'Final HRS still has incomplete volume retreats'
    n=retreat.subunits;assert np.all((a['Ag_units'][-1]==0)|(a['Ag_units'][-1]==n))
    a['global_elapsed_units']=s['initial_electronic_state']['elapsed_units']+a['elapsed_units']
    full.save_arrays(target,a);dump(target/'snapshots.json',snaps[1]);np.savez_compressed(target/'snapshots.npz',Ag_units=snaps[0],VN=cat.vn,subunits=n)
    np.savez_compressed(target/'actual_transport.npz',offsets=np.r_[0,np.cumsum([len(r) for r in cat.actual_routes])],
        nodes=np.concatenate(cat.actual_routes) if cat.actual_routes else np.array([],np.int16))
    if actions:
        with (target/'events.csv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=actions[0]);w.writeheader();w.writerows(actions)
    cp=target/'checkpoint';erosion_checkpoint(folder,s['final_event'],cp,s['final_event_exposure']);full.attach_electronics(cp,s['final_electronic_state'])
    saved=load(cp/'restart.npz');np.testing.assert_array_equal(saved['metal_id'],a['metal_id'][-1])
    saved.update(Ag_units=a['Ag_units'][-1],subunits=n);np.savez_compressed(cp/'restart.npz',**saved)
    dump(cp/'local_retreat.json',dict(parameters=asdict(retreat),trigger=s['trigger'],pending_microsteps=[],
        elapsed_units=s['elapsed_units'],observation_exposure=s['final_event_exposure'],final_event=s['final_event']))
    s.update(source_info(checkpoint,folder,number,'RESET'),source_chain=str(root.relative_to(ROOT)),
        network_params=asdict(network),final_Ag=float(a['Ag_total'][-1]),final_read_ohm=float(1/a['read_G_S'][-1]),
        final_metal_connected=bool(a['metal_connected'][-1]),lazy_network=cache_info,elapsed_seconds=time.perf_counter()-start)
    dump(target/'summary.json',s);print(f'cycle {number} RESET complete Ag={s["final_Ag"]} R={s["final_read_ohm"]:.6g} partial_sites={s["fractional_sites"]}',flush=True)
    return cp

def inherit_first(root):
    copies=[(BASE/'FORMING',root/'FORMING'),(BASE/'INIT_RESET',root/'INIT_RESET'),
        (BASE/'SET',root/'cycle_01/SET'),(FIRST_RESET,root/'cycle_01/RESET'),
        (BASE/'catalogues/SET',root/'cycle_01/catalogues/SET'),(BASE/'catalogues/RESET',root/'cycle_01/catalogues/RESET')]
    for source,dest in copies:
        if not dest.exists():shutil.copytree(source,dest)
    dump(root/'inherited_sources.json',{str(dest.relative_to(root)):dict(source=str(source.relative_to(ROOT)),states_sha256=sha(source/('sequence.npz' if 'catalogues' in str(source) else 'states.npz'))) for source,dest in copies})

def cycle_row(root,number):
    sf=root/f'cycle_{number:02d}'/'SET';rf=sf.parent/'RESET';a,b=load(sf/'states.npz'),load(rf/'states.npz')
    ss,rs=js(sf/'summary.json'),js(rf/'summary.json');limit=ss['clock']['circuit']['positive_limit_A']
    above=abs(a['current_A'][:101])>=.9*limit;on=np.flatnonzero(above[:-1]&above[1:]);j=int(on[0]) if len(on) else None
    connected=b['metal_connected'];loss=np.flatnonzero(connected[:-1]&~connected[1:]);br=int(loss[0]+1) if len(loss) else None
    metric=local.metrics(b,b['measured_A'])
    return dict(cycle=number,set_open_V=float(a['voltage_V'][j]) if j is not None else None,
        set_open_previous_V=float(a['voltage_V'][max(0,j-1)]) if j is not None else None,
        LRS_ohm=float(1/a['read_G_S'][-1]),HRS_ohm=float(1/b['read_G_S'][-1]),
        resistance_ratio=float(a['read_G_S'][-1]/b['read_G_S'][-1]),
        SET_final_Ag=float(a['Ag_total'][-1]),RESET_final_Ag=float(b['Ag_total'][-1]),
        SET_final_q=float(a['activity'][-1]),RESET_final_q=float(b['activity'][-1]),
        SET_metal_connected=bool(a['metal_connected'][-1]),RESET_metal_connected=bool(connected[-1]),
        RESET_break_V=float(b['voltage_V'][br]) if br is not None else None,
        RESET_early_drop_dec=metric['first_break_max_drop_dec'],RESET_outward_drop_dec=metric['outward_max_drop_dec'],
        RESET_reference_log_RMSE_dec=metric['whole_log_RMSE_dec'],
        SET_events=int(ss['final_event']),RESET_events=int(rs['final_event']),fractional_sites=rs['fractional_sites'])

def update_summary(root,total,complete=False):
    rows=[cycle_row(root,k) for k in range(1,total+1) if (root/f'cycle_{k:02d}/RESET/checkpoint/restart.npz').exists() and (root/f'cycle_{k:02d}/RESET/summary.json').exists()]
    dump(root/'chain_summary.json',dict(schema='ten_cycles40/chain/v1',complete=complete,requested_cycles=total,completed_cycles=len(rows),cycles=rows,
        sequence='FORMING -> INIT_RESET -> (SET -> RESET) x 10',
        protocol='Each SET 0 to +1.5 to 0 V; each RESET 0 to -1.5 to 0 V; 202 original C74 voltage samples per phase',
        parameter_fit=False,continuity='Physical geometry, VN, cumulative Ag balance, electronic state and RNG inherited between all stages',
        set_open_definition='First of two consecutive outward samples at >=90% of positive current compliance; previous sample brackets onset'))
    if rows:
        with (root/'cycle_summary.csv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    return rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cycles',type=int,default=10);ap.add_argument('--out',type=Path,default=OUT);ap.add_argument('--substeps',type=int,default=16);args=ap.parse_args()
    root=args.out.resolve();assert 1<=args.cycles<=100;root.mkdir(parents=True,exist_ok=True)
    positive,negative,network=full.common_parameters();grid=np.array(js(BASE/'parameters.json')['logq_grid'])
    retreat=local.Params(**js(FIRST_RESET/'summary.json')['parameters'])
    dump(root/'parameters.json',dict(positive=positive,negative=negative,network=asdict(network),local_retreat=asdict(retreat),logq_grid=grid.tolist(),substeps=args.substeps,
        first_cycle='Exact archived SET plus selected local-retreat RESET',units='0.5nm coarse Ag volumes; effective dwell units, not seconds'))
    inherit_first(root);checkpoint=root/'cycle_01/RESET/checkpoint';update_summary(root,args.cycles)
    for number in range(2,args.cycles+1):
        dump(root/'progress.json',dict(status='running',cycle=number,phase='SET',requested_cycles=args.cycles))
        checkpoint=run_set(root,number,checkpoint,positive,network,grid,args.substeps)
        dump(root/'progress.json',dict(status='running',cycle=number,phase='RESET',requested_cycles=args.cycles))
        checkpoint=run_reset(root,number,checkpoint,negative,network,grid,retreat,args.substeps)
        update_summary(root,args.cycles)
    rows=update_summary(root,args.cycles,complete=True)
    dump(root/'progress.json',dict(status='complete',completed_cycles=len(rows)))
    print(json.dumps(dict(complete=True,cycles=len(rows),out=str(root))),flush=True)

if __name__=='__main__':main()
