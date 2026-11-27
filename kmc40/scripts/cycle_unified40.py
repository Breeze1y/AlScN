"""Continuous circuit, electronic state and voltage-driven event clocks.

Spatial event sequences use declared reference-field probabilities. Only the
actual current-limited internal voltage drives their clock and electronic
state. Every dwell is integrated across events and electronic-state updates.
"""
from pathlib import Path
from dataclasses import dataclass,asdict,replace
import argparse,csv,json,math,time
import numpy as np
from scipy.optimize import least_squares
from cycle_activity40 import ActivityState,ActivityParams
from mixed_network40 import NetworkParams,prepare,solve
from fit_unified_reset40 import table,LOGQ,scores

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/cycle_unified40'
OLD=ROOT/'results/reset_unified40'

@dataclass(frozen=True)
class CircuitParams:
    network_scale:float=9.987870449406756
    series_contact_ohm:float=80.95973988295985
    leakage_S:float=8.944110346098897e-7
    positive_limit_A:float=.001
    negative_limit_A:float=.01
    # Optional effective Ag packets per dwell at positive current compliance.
    positive_exchange_per_dwell:float=0.

def read_circuit(gbase,q,voltage,params):
    material=params.network_scale*gbase
    series=material/(1+params.series_contact_ohm*material)
    total=series+params.leakage_S
    cap=params.positive_limit_A if voltage>=0 else params.negative_limit_A
    terminal=math.copysign(min(abs(voltage),cap/max(total,1e-30)),voltage)
    bulk=terminal/(1+params.series_contact_ohm*material)
    return dict(current_A=terminal*total,material_current_A=bulk*material,
                terminal_voltage_V=terminal,material_voltage_V=bulk,read_G_S=total)

class Phase:
    def __init__(self,folder,max_states=None):
        self.folder=Path(folder)
        with np.load(self.folder/'sequence.npz') as z:
            stop=min(len(z['metal_id']),max_states or len(z['metal_id']))
            self.ids=z['metal_id'][:stop].copy();self.vn=z['VN'].copy()
            self.clock=z['cumulative_hazard_clock'][:stop].copy()
            self.injected=z['injected'][:stop].copy() if 'injected' in z else np.repeat(int(z['initial_injected']),stop)
            self.returned=z['returned'][:stop].copy() if 'returned' in z else int(z['initial_returned'])+np.arange(stop)
        np.testing.assert_array_equal((self.ids>0).sum((1,2)),self.injected-self.returned)
        if not np.all(np.diff(self.clock)>0):raise ValueError('Event clocks must strictly increase')
        params=NetworkParams(**json.loads((OLD/'network_params.json').read_text()))
        self.params=params
        self.gt=table(self.ids>0,self.vn,params,self.folder/'network_activity_table.npz')
        self.logg=np.log(np.maximum(self.gt,1e-30))

    def g(self,index,q):
        return math.exp(float(np.interp(math.log10(max(q,1e-30)),LOGQ,self.logg[index])))

def simulate(phase,v,state,rate_scale,voltage_scale,circuit=CircuitParams(),substeps=8,exact=False):
    """One uninterrupted polarity branch. Input state is copied, never reset."""
    activity=state.copy();records=[];index=0;exposure=0.;exhausted=False;topologies={}
    def read(index,voltage):
        if exact:
            if index not in topologies:topologies[index]=prepare(phase.ids[index]>0,phase.vn)
            gg=solve(topologies[index],phase.params,activity=activity.q)['conductance_S']
        else:gg=phase.g(index,activity.q)
        return read_circuit(gg,activity.q,float(voltage),circuit)
    for point,voltage in enumerate(v):
        remaining=1.;mean_i=mean_terminal=mean_bulk=0.;count=0
        before=index
        while remaining>1e-10:
            r=read(index,voltage)
            u=abs(r['material_voltage_V'])
            speed=rate_scale*math.expm1(min(u/voltage_scale,60.))
            if voltage > 0 and circuit.positive_exchange_per_dwell:
                speed += (circuit.positive_exchange_per_dwell *
                          max(r['material_current_A'], 0.) / circuit.positive_limit_A)
            dt=min(remaining,1/substeps)
            if voltage>0:
                recovery=activity.recovery_hazard(r['material_voltage_V'],max(r['material_current_A'],0.))
                if recovery>0:dt=min(dt,.04/recovery)
            to_event=(phase.clock[index+1]-exposure)/speed if speed>0 and index+1<len(phase.clock) else math.inf
            hit=to_event<=dt+1e-12
            dt=max(0.,min(dt,to_event))
            if index+1==len(phase.clock) and speed>0:exhausted=True
            mean_i+=r['current_A']*dt;mean_terminal+=r['terminal_voltage_V']*dt;mean_bulk+=r['material_voltage_V']*dt
            activity.advance(r['material_voltage_V'],dt,current_A=max(r['material_current_A'],0.))
            exposure+=speed*dt;remaining-=dt;count+=1
            if hit:
                index+=1;exposure=float(phase.clock[index])
            if count>10000:raise RuntimeError('Event integration failed to advance')
        last=read(index,voltage)
        records.append(dict(point=point,voltage_V=float(voltage),current_A=mean_i,
            endpoint_current_A=last['current_A'],terminal_voltage_V=mean_terminal,
            material_voltage_V=mean_bulk,read_G_S=last['read_G_S'],activity=activity.q,
            D=activity.D,H=activity.H,event_index=index,events=index-before,
            event_exposure=exposure,Ag_total=int((phase.ids[index]>0).sum()),
            injected=int(phase.injected[index]),returned=int(phase.returned[index])))
    return records,activity,dict(final_event=index,event_exposure=exposure,
                               sequence_exhausted=exhausted,available_events=len(phase.clock)-1)

def experimental(which):
    file=OUT/'experimental'/('next_cycle_set.csv' if which=='reSET' else 'cycle074.csv')
    with file.open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    if which!='reSET':rows=[r for r in rows if r['phase']==which]
    return np.array([float(r['voltage_V']) for r in rows]),np.array([float(r['current_A']) for r in rows])

def save_phase(dest,phase,records,initial,final,runinfo,parameters,measured=None):
    if runinfo['sequence_exhausted']:raise ValueError('Extend the event sequence before saving a fitted result')
    dest=Path(dest)
    voltage=np.array([r['voltage_V'] for r in records])
    if dest.name.casefold()=='reset' and np.any(voltage>0):
        raise ValueError('RESET requires negative voltage; use SET_next for the next SET on Windows')
    if dest.name.casefold() in ('set','set_next') and np.any(voltage<0):
        raise ValueError('SET output requires positive voltage')
    dest.mkdir(parents=True,exist_ok=True)
    with (dest/'iv.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=records[0]);writer.writeheader();writer.writerows(records)
    arrays={k:np.array([r[k] for r in records]) for k in records[0]}
    arrays.update(metal_id=phase.ids[arrays['event_index']],VN=phase.vn)
    if measured is not None:arrays['measured_A']=measured
    np.savez_compressed(dest/'states.npz',**arrays)
    info=dict(source=str(phase.folder),initial_electronic_state=initial.to_dict(),
              final_electronic_state=final.to_dict(),clock=parameters,**runinfo)
    if measured is not None:info['metrics']=scores(arrays['voltage_V'],measured,arrays['current_A'])
    (dest/'summary.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
    return info

def fit_set():
    phase=Phase(OUT/'growth_initial',701);v,measured=experimental('SET')
    original=ActivityState.from_retention(OLD/'retention/HRS')
    # Unit attempt clock near the first true bridge gives a reproducible start.
    bridge=402.;options=[];start=time.monotonic()
    for vs in [.025,.04,.06]:
        scale=bridge/max(np.expm1(v[:23]/vs).sum(),1e-30)
        p0=[math.log(scale),math.log(vs),math.log(.5)]
        bounds=([p0[0]-5,math.log(.012),math.log(.08)],
                [p0[0]+5,math.log(.12),math.log(4.)])
        def run(p,steps=4):
            st=original.copy();st.params=replace(st.params,recovery_rate=math.exp(p[2]),recovery_voltage_scale_V=.3)
            rr,end,meta=simulate(phase,v,st,math.exp(p[0]),math.exp(p[1]),substeps=steps)
            return rr,st,end,meta
        def residual(p):
            rr,st,end,meta=run(p);pred=np.array([r['current_A'] for r in rr]);ok=v!=0
            error=np.log10(np.maximum(abs(pred[ok]),1e-30)/abs(measured[ok]))
            # Same-circle low-bias return points directly constrain the LRS.
            tail=np.array([r['current_A'] for r in rr])[198:201]
            return np.r_[error,2*np.log10(np.maximum(abs(tail),1e-30)/abs(measured[198:201])),
                         5. if meta['sequence_exhausted'] else 0.]
        opt=least_squares(residual,p0,bounds=bounds,max_nfev=100,diff_step=.002)
        rr,st,end,meta=run(opt.x,8)
        metric=scores(v,measured,np.array([r['current_A'] for r in rr]))
        options.append((metric['all'],opt.x,rr,st,end,meta,metric))
        pp=opt.x
        candidate=dict(rate_scale=math.exp(pp[0]),voltage_scale_V=math.exp(pp[1]),
                       recovery_rate=math.exp(pp[2]),recovery_voltage_scale_V=.3,circuit=asdict(CircuitParams()))
        if not meta['sequence_exhausted']:
            save_phase(OUT/'set_candidates'/str(vs),phase,rr,st,end,meta,candidate,measured)
        print('SET candidate',vs,metric,meta,'seconds',time.monotonic()-start,flush=True)
    _,p,rr,st,end,meta,metric=min(options,key=lambda x:x[0]+5*x[5]['sequence_exhausted'])
    params=dict(rate_scale=math.exp(p[0]),voltage_scale_V=math.exp(p[1]),
                recovery_rate=math.exp(p[2]),recovery_voltage_scale_V=.3,circuit=asdict(CircuitParams()),
                driving_voltage='material voltage after compliance and effective contact drop')
    info=save_phase(OUT/'SET',phase,rr,st,end,meta,params,measured)
    print('SELECTED SET',json.dumps(info),flush=True)

def calibrated_initial(params):
    """Re-express the saved two electronic damage components under a trial law.

    Calibration changes coordinate units, not the physical saved q. Within a
    simulated cycle parameters and coordinates are inherited without remapping.
    """
    old=ActivityState.from_retention(OLD/'retention/HRS')
    low=old.params.activity_exponent*math.log1p(old.D/old.params.dose_scale)
    high=old.params.high_field_rate*old.H
    return ActivityState(params.dose_scale*math.expm1(low/params.activity_exponent),
                         high/params.high_field_rate,params)

def fit_cycle(extended=False):
    growth=Phase(OUT/'growth_initial',701);erosion=Phase(OUT/'erosion_cycle')
    sv,si=experimental('SET');rv,ri=experimental('RESET')
    prior=json.loads((OUT/'SET/summary.json').read_text(encoding='utf-8'))
    pos=prior['clock'];old=ActivityState.from_retention(OLD/'retention/HRS').params
    required=int(prior['final_event']);choices=[];start=time.monotonic()
    assert np.array_equal(growth.ids[required],erosion.ids[0])
    def run(p,steps=4):
        par=replace(old,recovery_rate=math.exp(p[2]),recovery_voltage_scale_V=.3)
        if len(p)>3:par=replace(par,dose_scale=math.exp(p[3]),activity_exponent=math.exp(p[4]),high_field_rate=math.exp(p[5]))
        initial=calibrated_initial(par)
        a,mid,am=simulate(growth,sv,initial,pos['rate_scale'],pos['voltage_scale_V'],substeps=steps)
        b,end,bm=simulate(erosion,rv,mid,math.exp(p[0]),math.exp(p[1]),substeps=steps)
        return a,b,initial,mid,end,am,bm
    def residual(p):
        a,b,initial,mid,end,am,bm=run(p)
        aa=np.array([x['current_A'] for x in a]);bb=np.array([x['current_A'] for x in b])
        er=np.log10(np.maximum(abs(bb[rv!=0]),1e-30)/abs(ri[rv!=0]))
        es=np.log10(np.maximum(abs(aa[sv!=0]),1e-30)/abs(si[sv!=0]))
        return np.r_[er,es,.2*np.diff(er),
                     5. if bm['sequence_exhausted'] or am['sequence_exhausted'] else 0.,
                     .2*abs(am['final_event']-required)]
    seeds=[(.1948,.28265,.5),(1.,.3,1.),(.05,.2,.3)]
    simple=OUT/'joint_simple.json'
    if extended and simple.exists():
        saved=json.loads(simple.read_text());p=saved['opt_log_parameters']
        seeds=[tuple(np.exp(p[:3])),(np.exp(p[0])*2,np.exp(p[1]),.5)]
    for scale,vs,recovery in seeds:
        p=np.log([scale,vs,recovery]).tolist()
        lo=np.log([1e-6,.04,.03]).tolist();hi=np.log([100.,1.,4.]).tolist()
        if extended:
            p+=np.log([old.dose_scale,old.activity_exponent,old.high_field_rate]).tolist()
            lo+=np.log([.05,.1,.005]).tolist();hi+=np.log([300.,4.,5.]).tolist()
        opt=least_squares(residual,p,bounds=(lo,hi),max_nfev=160,diff_step=.002)
        result=run(opt.x,16)
        a,b,initial,mid,end,am,bm=result
        sm=scores(sv,si,np.array([x['current_A'] for x in a]));rm=scores(rv,ri,np.array([x['current_A'] for x in b]))
        valid=not (am['sequence_exhausted'] or bm['sequence_exhausted']) and am['final_event']==required
        choices.append((sm['all']+rm['all']+(0 if valid else 10),opt.x,result,sm,rm))
        print('Joint candidate',extended,sm['all'],rm['all'],am,bm,'qLRS',mid.q,'seconds',time.monotonic()-start,flush=True)
    _,p,result,sm,rm=min(choices,key=lambda a:a[0])
    a,b,initial,mid,end,am,bm=result
    if am['final_event']!=required:raise ValueError('Joint fit changed SET geometry; regenerate negative sequence')
    params=dict(rate_scale=pos['rate_scale'],voltage_scale_V=pos['voltage_scale_V'],
                recovery_rate=initial.params.recovery_rate,recovery_voltage_scale_V=.3,
                circuit=asdict(CircuitParams()),activity_params=asdict(initial.params),
                driving_voltage='material voltage after compliance and effective contact drop',
                calibration='All nonzero C74 SET and RESET samples; parity subsets are descriptive, not held out')
    save_phase(OUT/'SET',growth,a,initial,mid,am,params,si)
    negative=dict(params,rate_scale=math.exp(p[0]),voltage_scale_V=math.exp(p[1]))
    save_phase(OUT/'RESET',erosion,b,mid,end,bm,negative,ri)
    info=dict(extended=extended,opt_log_parameters=p.tolist(),SET=sm,RESET=rm,
              initial_q=initial.q,LRS_q=mid.q,HRS_q=end.q,
              SET_final_event=am['final_event'],RESET_final_event=bm['final_event'],
              negative_clock=negative,calibration=negative['calibration'])
    (OUT/('joint_extended.json' if extended else 'joint_simple.json')).write_text(json.dumps(info,indent=2),encoding='utf-8')
    print('SELECTED JOINT',json.dumps(info),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--fit-set',action='store_true')
    parser.add_argument('--fit-cycle',action='store_true');parser.add_argument('--extended',action='store_true')
    args=parser.parse_args()
    if args.fit_set:fit_set()
    if args.fit_cycle:fit_cycle(args.extended)
