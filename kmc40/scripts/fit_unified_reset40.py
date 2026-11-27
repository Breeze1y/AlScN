"""Fit a unified passive network and a separable stochastic RESET clock.

Only even nonzero C74 points train parameters. Odd points check interpolation.
The activity interpolation approximates network solutions, not measured I-V.
"""
from pathlib import Path
from dataclasses import asdict, replace
import argparse, csv, hashlib, json, time
import numpy as np
from scipy.optimize import least_squares
from mixed_network40 import NetworkParams, prepare, solve

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/reset_unified40'
BASE=ROOT.parent/'baseline/diff/fit_baselines/AlScN_40nm_file074.csv'
SOURCE=ROOT/'results/vacancy_assisted40/source_cycle/assisted/states.npz'
LOGQ=np.linspace(-10,0,21)

def baseline():
    with BASE.open(encoding='utf-8-sig') as f:
        rows=[r for r in csv.DictReader(f) if r['phase']=='RESET']
    return np.array([float(r['voltage_V']) for r in rows]),np.array([float(r['current_A']) for r in rows])

def scores(v,observed,pred):
    valid=v!=0; idx=np.arange(len(v))
    e=np.log10(np.maximum(abs(pred),1e-30)/np.maximum(abs(observed),1e-30))
    out={k:float(np.sqrt(np.mean(e[valid&m]**2))) for k,m in
         [('all',np.ones(len(v),bool)),('train',idx%2==0),('interleaved',idx%2==1),('out',idx<=100),('back',idx>=101)]}
    mask=(v[1:]!=0)&(v[:-1]!=0)
    jumps=np.abs(np.diff(np.log10(np.maximum(abs(pred),1e-30))))
    out['max_adjacent_jump_dec']=float(jumps[mask].max())
    return out

def table(sequence,vn,params,cache):
    signature=hashlib.sha256(sequence.tobytes()+vn.tobytes()+json.dumps(asdict(params),sort_keys=True).encode()+LOGQ.tobytes()).hexdigest()
    if cache.exists():
        with np.load(cache) as z:
            if str(z['signature'])==signature:return z['G_S']
    ans=np.zeros((len(sequence),len(LOGQ)))
    seen={};start=time.monotonic()
    for i,metal in enumerate(sequence):
        key=metal.tobytes()
        if key in seen:ans[i]=ans[seen[key]];continue
        topo=prepare(metal,vn)
        for j,logq in enumerate(LOGQ):
            ans[i,j]=solve(topo,params,activity=10.**logq)['conductance_S']
        seen[key]=i
        if i%25==0:print(f'table {cache.stem} {i}/{len(sequence)} elapsed={time.monotonic()-start:.1f}s',flush=True)
    np.savez_compressed(cache,G_S=ans,log_activity=LOGQ,signature=signature,params=json.dumps(asdict(params)))
    return ans

def activity(v,p,field_scale):
    dose=np.cumsum(abs(v))
    high=np.cumsum(abs(v)*np.exp((abs(v)-1.5)/field_scale))
    # q is an effective participating electronic-state fraction, not VN density.
    return np.exp(-p[2]*np.log1p(dose/np.exp(p[1]))-np.exp(p[3])*high)

def event_values(gtable,q):
    x=np.clip(np.log10(np.maximum(q,1e-300)),LOGQ[0],LOGQ[-1])
    ix=np.minimum(np.searchsorted(LOGQ,x,side='right')-1,len(LOGQ)-2)
    w=(x-LOGQ[ix])/(LOGQ[ix+1]-LOGQ[ix])
    lg=np.log(np.maximum(gtable,1e-30))
    return np.exp((1-w[:,None])*lg[:,ix].T+w[:,None]*lg[:,ix+1].T)

def clock_average(values,clock,start,end):
    """Exact dwell average for piecewise constant event states and fixed q per dwell."""
    widths=np.diff(clock)
    integ=np.c_[np.zeros(len(values)),np.cumsum(values[:,:-1]*widths[None,:],axis=1)]
    def primitive(d):
        ix=np.clip(np.searchsorted(clock,d,side='right')-1,0,len(clock)-1)
        return integ[np.arange(len(d)),ix]+values[np.arange(len(d)),ix]*(d-clock[ix]),ix
    hi,ix=primitive(end);lo,_=primitive(start)
    endpoint=values[np.arange(len(ix)),ix]
    avg=np.divide(hi-lo,end-start,out=endpoint.copy(),where=end>start)
    return avg,endpoint,ix

def evaluate(v,gtable,p,field_scale,clock=None):
    q=activity(v,p,field_scale); vals=event_values(gtable,q)
    vals*=np.exp(p[0])
    if len(p)>7:
        # Localized electrode/filament contact resistance is in series with
        # each instantaneous material network, before acquisition averaging.
        vals=vals/(1+np.exp(p[7])*vals)
    if clock is None:
        g=vals[np.arange(len(v)),np.arange(len(v))];endg=g;ix=np.arange(len(v))
    else:
        raw=np.expm1(abs(v)/np.exp(p[6]))
        # Parameterize the final exposure as a fraction of the simulated clock.
        # No frozen continuation beyond the last simulated event is permitted.
        rate=np.exp(p[5])*clock[-1]*raw/max(raw.sum(),1e-30)
        end=np.cumsum(rate);start=end-rate
        g,endg,ix=clock_average(vals,clock,start,end)
    g=g+np.exp(p[4]);endg=endg+np.exp(p[4])
    return np.clip(v*g,-.01,.01),q,g,endg,ix

def fit(v,measured,gtable,clock=None,contact=False):
    train=(v!=0)&(np.arange(len(v))%2==0); idx=np.flatnonzero(train)
    choices=[]
    for field in [.04,.12]:
        for d0 in [3.,25.]:
            p=[0.,np.log(d0),2.,-3.,np.log(1e-6)]
            lo=[np.log(.1),-7,.05,-14,np.log(1e-9)]
            hi=[np.log(10),9,10,5,np.log(1e-4)]
            if clock is not None:
                scale=float(clock[min(40,len(clock)-1)])/max(np.expm1(abs(v[:30])/.3).sum(),1e-30)
                fraction=np.clip(scale*np.expm1(abs(v)/.3).sum()/clock[-1],1e-4,.8)
                p += [np.log(fraction),np.log(.3)]
                lo += [np.log(1e-5),np.log(.08)]
                hi += [np.log(.98),np.log(1.5)]
                if contact:
                    p += [np.log(100.)];lo += [np.log(.01)];hi += [np.log(1500.)]
            def residual(x):
                pred=evaluate(v,gtable,x,field,clock)[0]
                e=np.log10(np.maximum(abs(pred[idx]),1e-30)/abs(measured[idx]))
                # Penalize mismatched change between TRAIN points only.
                return np.r_[e,.3*np.diff(e)]
            opt=least_squares(residual,p,bounds=(lo,hi),max_nfev=250,diff_step=1e-3)
            pred,q,g,endg,ix=evaluate(v,gtable,opt.x,field,clock)
            metric=scores(v,measured,pred)
            choices.append((metric['train'],opt.x,field,metric,pred,q,g,endg,ix))
    return min(choices,key=lambda t:t[0])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sequence');ap.add_argument('--name',default='fixed_morphology')
    ap.add_argument('--trap',type=float,default=7.748e-5);ap.add_argument('--beta',type=float,default=2.)
    ap.add_argument('--contact',action='store_true')
    args=ap.parse_args();dest=OUT/args.name;dest.mkdir(parents=True,exist_ok=True)
    v,measured=baseline()
    if args.sequence:
        source=Path(args.sequence)
        with np.load(source) as z:
            key='metal_id' if 'metal_id' in z else 'sequence'
            ids=z[key].copy();seq=ids>0;vn=z['VN'].copy();clock=z['cumulative_hazard_clock'].copy()
    else:
        source=SOURCE
        with np.load(source) as z:
            mask=(z['phase'].astype(str)=='RESET')&(z['cycle']==1)
            ids=z['metal_id'][mask].copy();seq=ids>0;vn=z['VN'].copy();clock=None
            np.testing.assert_allclose(v,z['voltage_V'][mask],atol=1e-14)
    par=NetworkParams(trap_prefactor_S=args.trap,trap_beta_nm=args.beta,leakage_S=0.)
    gt=table(seq,vn,par,dest/'network_activity_table.npz')
    _,p,field,metric,pred,q,g,endg,ix=fit(v,measured,gt,clock,args.contact)
    print(args.name,metric,flush=True)
    np.savez_compressed(dest/'fit.npz',voltage_V=v,measured_A=measured,current_A=pred,activity=q,
                        conductance_S=g,endpoint_G_S=endg,event_index=ix,metal_id=ids[ix],VN=vn,
                        parameters=p,field_scale_V=field)
    with (dest/'iv.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(['voltage_V','measured_A','current_A','endpoint_current_A','activity','event_index','Ag_total'])
        w.writerows(zip(v,measured,pred,v*endg,q,ix,seq[ix].sum((1,2))))
    info=dict(name=args.name,source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              network=asdict(par),parameters=dict(network_scale=np.exp(p[0]),dose_scale=np.exp(p[1]),
              activity_exponent=p[2],high_field_rate=np.exp(p[3]),leakage_S=np.exp(p[4]),field_scale_V=field,
              series_contact_ohm=np.exp(p[7]) if len(p)>7 else 0.),
              metrics=metric,measurement='end of sample' if clock is None else 'full dwell average; effective acquisition assumption',
              activity='state participation fraction, updated throughout RESET; fixed VN coordinates',
              clock='frozen original frames' if clock is None else dict(
                  rate_scale=np.exp(p[5])*clock[-1]/np.expm1(abs(v)/np.exp(p[6])).sum(),
                  voltage_scale_V=np.exp(p[6]),final_exposure_fraction=np.exp(p[5])),
              final_event=int(ix[-1]),available_events=len(seq)-1,final_read_ohm=float(1/endg[-1]),
              validation='even points train; odd points same-curve interpolation, not independent device validation')
    (dest/'summary.json').write_text(json.dumps(info,indent=2),encoding='utf-8')

if __name__=='__main__':main()
