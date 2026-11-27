"""Local, finite-duration Ag volume loss near the first prospective rupture.

The saved spatial KMC catalogue supplies event order and hazard thresholds.
Current, material voltage, electronic activity and event exposure are coupled
again. No measured-current interpolation or output smoothing is applied.
"""
from pathlib import Path
from dataclasses import dataclass,asdict
import argparse,csv,heapq,json,math,time,hashlib
import numpy as np
from cycle_activity40 import ActivityState
from cycle_unified40 import CircuitParams,read_circuit
from mixed_network40 import NetworkParams
from form_hourglass40 import span
import local_fraction_network as net

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/full_chain_refined'
OUT=ROOT/'results/local_retreat40'
def load(p):
    with np.load(p) as z:return {k:z[k].copy() for k in z.files}
def js(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

@dataclass(frozen=True)
class Params:
    enabled:bool=True
    subunits:int=16
    retreat_duration:float=4.
    local_window:float=12.
    radius_nm:float=1.5
    contact_exponent:float=4.

class Catalogue:
    def __init__(self):
        self.folder=BASE/'catalogues/RESET';z=load(self.folder/'sequence.npz')
        self.ids=z['metal_id'];self.vn=z['VN'];self.clock=z['cumulative_hazard_clock']
        self.injected=int(z['initial_injected']);self.returned=int(z['initial_returned'])
        t=load(self.folder/'network_activity_table.npz');self.logq=t['log_activity'];self.logg=np.log(np.maximum(t['G_S'],1e-30))
        with (self.folder/'events.csv').open(encoding='utf-8',newline='') as f:self.events=list(csv.DictReader(f))
        self.base=load(BASE/'RESET/states.npz');self.summary=js(BASE/'RESET/summary.json')
        self.network=NetworkParams(**self.summary['network_params'])
        self.topology=net.prepare(self.ids[0]>0,self.vn)
        routes=load(self.folder/'transport.npz')
        self.routes=[routes['nodes'][a:b] for a,b in zip(routes['offsets'][:-1],routes['offsets'][1:])]
    def binary_g(self,index,q):return math.exp(float(np.interp(math.log10(max(q,1e-30)),self.logq,self.logg[index])))

def metrics(a,ref):
    v=a['voltage_V'];x=abs(v);nz=x>0;y=np.log10(np.maximum(abs(a['current_A']),1e-30));r=np.log10(np.maximum(abs(ref),1e-30))
    outward=(np.arange(len(v))<=100)&nz;ix=np.flatnonzero(outward);drops=-np.diff(y[ix]);loc=(x[ix[:-1]]>=.15-1e-12)&(x[ix[1:]]<=.36+1e-12)
    windows=[]
    for lo,hi in [(.15,.3),(.3,.6),(.6,.9),(.9,1.2),(1.2,1.5)]:
        m=outward&(x>=lo-1e-12)&(x<=hi+1e-12)
        sm=float(np.polyfit(x[m],y[m],1)[0]);sr=float(np.polyfit(x[m],r[m],1)[0]);windows.append(dict(absV_min=lo,absV_max=hi,model=sm,measured=sr))
    br=np.flatnonzero(~a['metal_connected']);k=int(br[0]) if len(br) else None
    mid=outward&(x>=.3-1e-12)&(x<=.9+1e-12)
    local=float(max(0,drops[loc].max()))
    whole=float(np.sqrt(np.mean((y[nz]-r[nz])**2)))
    return dict(whole_log_RMSE_dec=whole,outward_log_RMSE_dec=float(np.sqrt(np.mean((y[outward]-r[outward])**2))),
        first_break_max_drop_dec=local,outward_max_drop_dec=float(max(0,drops.max())),
        mid_log_RMSE_dec=float(np.sqrt(np.mean((y[mid]-r[mid])**2))),windows=windows,
        first_disconnected_sample=k,first_disconnected_voltage_V=float(v[k]) if k is not None else None,
        final_read_ohm=float(1/a['read_G_S'][-1]),final_Ag=float(a['Ag_total'][-1]),final_q=float(a['activity'][-1]),
        ranking_score=local+0.5*whole)

def simulate(cat,p,substeps=16):
    cfg=cat.summary['clock'];circ=CircuitParams(**cfg['circuit']);state=ActivityState.from_dict(cat.summary['initial_electronic_state'])
    units=((cat.ids[0]>0)*p.subunits).astype(np.uint8);initial_units=int(units.sum())
    index=0;exposure=0.;now=0.;queue=[];serial=0;trigger=None;records=[];states=[];actions=[];snapshots=[];snapshot_meta=[];exact_solves=0
    connected=bool(span(units>0).any());actual_break_recorded=False
    def read(v):
        nonlocal exact_solves
        if not np.array_equal(units,(cat.ids[index]>0)*p.subunits):
            g=net.solve(cat.topology,units.astype(float)/p.subunits,cat.network,activity=state.q,contact_exponent=p.contact_exponent)['conductance_S'];exact_solves+=1
        else:g=cat.binary_g(index,state.q)
        return read_circuit(g,state.q,float(v),circ)
    def snap(label,point,v):
        snapshots.append(units.copy());snapshot_meta.append(dict(label=label,point=point,time=now,voltage_V=float(v),q=state.q,event_index=index,current_A=read(v)['current_A']))
    def remove(y,x,amount,point,v,kind,event_id,step):
        nonlocal connected,actual_break_recorded
        before=connected;previous=int(units[y,x]);assert previous>=amount
        route=cat.routes[event_id-1]
        assert int(route[0])==y*units.shape[1]+x
        # Released Ag uses an existing dielectric return route. In these runs
        # no path beyond its source passes through any residual Ag volume.
        assert not np.any(units.ravel()[route[1:]])
        trial=units.copy();trial[y,x]-=amount;after=bool(span(trial>0).any())
        if before and not after and not actual_break_recorded:snap('metal_break_before',point,v)
        units[y,x]-=amount;connected=after
        actions.append(dict(time=now,point=point,voltage_V=float(v),event_index=event_id,y=y,x=x,kind=kind,
            unit_step=step,removed_units=amount,remaining_units=int(units[y,x]),Ag_total=float(units.sum()/p.subunits),
            returned_effective=float(cat.returned+(initial_units-int(units.sum()))/p.subunits),q=state.q,
            connected_before=before,connected_after=after))
        if before and not after and not actual_break_recorded:
            actual_break_recorded=True;snap('metal_break_after',point,v)
    for point,v in enumerate(cat.base['voltage_V']):
        remaining=1.;mean_i=mean_terminal=mean_bulk=0.;before=index;micro_before=len(actions);count=0
        while remaining>1e-10:
            rd=read(v);speed=cfg['rate_scale']*math.expm1(min(abs(rd['material_voltage_V'])/cfg['voltage_scale_V'],60.))
            if index+1>=len(cat.clock) and speed>0:raise RuntimeError('Spatial catalogue exhausted')
            to_event=(cat.clock[index+1]-exposure)/speed if speed>0 and index+1<len(cat.clock) else math.inf
            to_micro=max(0.,queue[0][0]-now) if queue else math.inf
            cap=min(remaining,1/substeps);dt=max(0.,min(cap,to_event,to_micro))
            event_hit=to_event<=dt+1e-12;micro_hit=to_micro<=dt+1e-12
            mean_i+=rd['current_A']*dt;mean_terminal+=rd['terminal_voltage_V']*dt;mean_bulk+=rd['material_voltage_V']*dt
            state.advance(rd['material_voltage_V'],dt,current_A=max(rd['material_current_A'],0.))
            exposure+=speed*dt;remaining-=dt;now+=dt;count+=1
            if micro_hit:
                while queue and queue[0][0]<=now+1e-10:
                    due,order,y,x,ev,step=heapq.heappop(queue)
                    # Keep binary/partial read mode valid while recording a last microstep.
                    remove(y,x,1,point,v,'retreat_microstep',ev,step)
            if event_hit:
                index+=1;exposure=float(cat.clock[index]);e=cat.events[index-1];y,x=int(e['y']),int(e['x'])
                assert units[y,x]==p.subunits,'A catalogue site may be initiated only once'
                if p.enabled and trigger is None and connected:
                    trial=units.copy();trial[y,x]=0
                    if not span(trial>0).any():
                        trigger=dict(event_index=index,point=point,time=now,voltage_V=float(v),y=y,x=x,until=now+p.local_window)
                        snap('retreat_onset_before',point,v)
                local=p.enabled and trigger is not None and now<=trigger['until']+1e-12 and math.hypot(y-trigger['y'],x-trigger['x'])*.5<=p.radius_nm+1e-12
                if local:
                    for step in range(2,p.subunits+1):
                        serial+=1;heapq.heappush(queue,(now+p.retreat_duration*(step-1)/(p.subunits-1),serial,y,x,index,step))
                    remove(y,x,1,point,v,'retreat_start',index,1)
                    if trigger['event_index']==index:snap('retreat_onset_after',point,v)
                else:remove(y,x,p.subunits,point,v,'whole_cell',index,1)
            if count>20000:raise RuntimeError('Integration did not advance')
        last=read(v);total=float(units.sum()/p.subunits);returned=cat.returned+(initial_units-int(units.sum()))/p.subunits
        assert abs(total-cat.injected+returned)<1e-12
        records.append(dict(point=point,voltage_V=float(v),current_A=mean_i,endpoint_current_A=last['current_A'],
            terminal_voltage_V=mean_terminal,material_voltage_V=mean_bulk,read_G_S=last['read_G_S'],activity=state.q,D=state.D,H=state.H,
            event_index=index,events=index-before,microevents=len(actions)-micro_before,event_exposure=exposure,
            Ag_total=total,Ag_occupied_cells=int(np.count_nonzero(units)),injected=cat.injected,returned=returned,
            metal_connected=connected,pending_retreats=len({item[4] for item in queue}),elapsed_units=now))
        states.append(units.copy())
    arrays={k:np.array([r[k] for r in records]) for k in records[0]}
    arrays.update(Ag_units=np.asarray(states),VN=cat.vn,subunits=np.array(p.subunits),initial_metal_id=cat.ids[0],
        metal_id=np.where(np.asarray(states)>0,cat.ids[0],0),measured_A=cat.base['measured_A'],branch=cat.base['branch'])
    summary=dict(schema='local_retreat40/run/v1',parameters=asdict(p),clock=cfg,network=asdict(cat.network),integration_substeps=substeps,
        source_chain='results/full_chain_refined',source_checkpoint='results/full_chain_refined/SET/checkpoint',
        initial_electronic_state=cat.summary['initial_electronic_state'],final_electronic_state=state.to_dict(),
        trigger=trigger,final_event=index,final_event_exposure=exposure,elapsed_units=now,pending_microsteps=[list(item) for item in queue],
        microsteps=len(actions),fractional_sites=len({r['event_index'] for r in actions if r['kind']=='retreat_start'}),exact_fractional_solves=exact_solves,
        metrics=metrics(arrays,cat.base['measured_A']),
        semantics='Catalogue entries initiate dissolution; Ag effective volume returns only on actual microsteps. Fixed spatial q=1 event order; dynamic q and event timing recomputed. Every volume removal checks its inherited dielectric return route; excluding its source, no route crosses remaining Ag in these runs.',
        time_units='effective sampling units; no conversion to seconds',Ag_units='0.5 nm coarse-volume packets, subdivided into effective volumes, not atoms')
    return arrays,summary,actions,(np.asarray(snapshots),snapshot_meta)

def run(name,p,steps=16,save=True,force=False):
    target=OUT/name
    if save and not force and (target/'summary.json').exists():return js(target/'summary.json')
    cat=Catalogue();start=time.perf_counter();a,s,actions,snaps=simulate(cat,p,steps);s['elapsed_seconds']=time.perf_counter()-start
    if save:
        target.mkdir(parents=True,exist_ok=True);np.savez_compressed(target/'states.npz',**a)
        s['source_sha256']={str(f.relative_to(ROOT)):sha(f) for f in [BASE/'RESET/states.npz',BASE/'RESET/summary.json',BASE/'SET/checkpoint/restart.npz',cat.folder/'sequence.npz',cat.folder/'network_activity_table.npz']}
        s['code_sha256']={str(f.relative_to(ROOT)):sha(f) for f in [Path(__file__),ROOT/'scripts/local_fraction_network.py']}
        dump(target/'summary.json',s);dump(target/'snapshots.json',snaps[1]);np.savez_compressed(target/'snapshots.npz',Ag_units=snaps[0],VN=cat.vn,subunits=p.subunits)
        for path,rows in [(target/'events.csv',actions),(target/'iv.csv',[{k:np.asarray(v[i]).item() for k,v in a.items() if np.asarray(v).ndim==1 and len(v)==len(a['point'])} for i in range(len(a['point']))])]:
            with path.open('w',encoding='utf-8-sig',newline='') as f:
                w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    print(name,json.dumps(dict(metrics=s['metrics'],sites=s['fractional_sites'],seconds=s['elapsed_seconds'])),flush=True)
    return s

def save_final_checkpoint(name):
    from rate_reset40 import checkpoint_at_event
    from run_full_chain import attach_electronics
    target=OUT/name;s=js(target/'summary.json');a=load(target/'states.npz')
    assert not s['pending_microsteps'],'Pending partial cells require the scheduler to resume'
    n=s['parameters']['subunits'];assert np.all((a['Ag_units'][-1]==0)|(a['Ag_units'][-1]==n))
    checkpoint_at_event(BASE/'catalogues/RESET',s['final_event'],target/'checkpoint',s['final_event_exposure'])
    attach_electronics(target/'checkpoint',s['final_electronic_state'])
    saved=load(target/'checkpoint/restart.npz')
    np.testing.assert_array_equal(saved['metal_id'],a['metal_id'][-1])
    saved.update(Ag_units=a['Ag_units'][-1],subunits=n)
    np.savez_compressed(target/'checkpoint/restart.npz',**saved)
    dump(target/'checkpoint/local_retreat.json',dict(schema='local_retreat40/checkpoint/v1',parameters=s['parameters'],
        trigger=s['trigger'],pending_microsteps=s['pending_microsteps'],elapsed_units=s['elapsed_units'],
        observation_exposure=s['final_event_exposure'],final_event=s['final_event'],
        continuation='Completed RESET at zero bias. All local retreats finished. Binary checkpoint can initialize the next SET; electronic activity and RNG are preserved. A new RESET starts a new rupture trigger.'))
    return target/'checkpoint'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--name',default='local_tau4_p4');ap.add_argument('--tau',type=float,default=4.);ap.add_argument('--window',type=float,default=12.);ap.add_argument('--radius',type=float,default=1.5);ap.add_argument('--power',type=float,default=4.);ap.add_argument('--subunits',type=int,default=16);ap.add_argument('--substeps',type=int,default=16);ap.add_argument('--disabled',action='store_true');ap.add_argument('--force',action='store_true');ap.add_argument('--checkpoint',action='store_true');a=ap.parse_args()
    assert a.tau>0 and a.window>0 and 2<=a.subunits<=255 and a.radius>0 and a.power>0
    run(a.name,Params(not a.disabled,a.subunits,a.tau,a.window,a.radius,a.power),a.substeps,force=a.force)
    if a.checkpoint:print(save_final_checkpoint(a.name))
if __name__=='__main__':main()
