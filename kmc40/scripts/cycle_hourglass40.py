"""Grow/erode individual surface sites of the stochastic hourglass; physical depths never enter event selection."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import convolve, distance_transform_edt, label
from scipy.optimize import least_squares

import inspect_midspan_drive as thermal
import spatial_reset40_erosion as prior

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from forming_iv import ElectronicConfig
from dynamic_kmc import dynamic_network

BASE=ROOT/'results/hourglass40'
OUT=BASE/'cycle'
CROSS=np.array([[0,1,0],[1,0,1],[0,1,0]])
G0=.0082123
MATERIAL=json.loads((ROOT/'configs/reset_set_20nm_pilot.json').read_text())['electronics']
NETWORK=ElectronicConfig(**MATERIAL)
PARAM_NAMES=('oxidation_onset_V','low_field_dose_scale_V','low_field_dose_amplitude',
             'high_field_onset_V','high_field_width_V','extra_high_field_dose',
             'last_bridge_event_relative_dwell','log10_parallel_leakage_S','post_break_relative_dwell')


def clusters(a):
    labs,_=label(a)
    top=np.unique(labs[0]);bot=np.unique(labs[-1])
    t=np.isin(labs,top[top>0]);b=np.isin(labs,bot[bot>0])
    return t,b,t&b


def read_network(a):
    padded=np.pad(a,((1,1),(0,0)),constant_values=True).astype(np.uint8)*2
    return dynamic_network(padded,.5,NETWORK)['conductance_S']


def neighbours(a):
    return convolve(a.astype(float),CROSS,mode='constant')


def draw_site(score,rng):
    valid=np.isfinite(score)
    if not valid.any():
        raise ValueError('Empty surface event list')
    x=np.argwhere(valid); weights=np.exp(score[valid]-np.max(score[valid]))
    return tuple(map(int,x[rng.choice(len(x),p=weights/weights.sum())]))


def dissolve(a,rng):
    t,b,s=clusters(a); cn=neighbours(a)
    if s.any():
        _,power=thermal.electrical(s); temp=thermal.temperature(a,power,3)
        eligible=s&(cn<4)&(power>power.max()*1e-9)
        score=26*temp+1.8*np.log(np.maximum(power/power.max(),1e-12))-.65*cn
    else:
        # Only exposed facing surfaces are eligible after break. A local field proxy is
        # reduced as separation from the opposite conductor increases.
        distance=np.where(t,distance_transform_edt(~b),distance_transform_edt(~t))
        eligible=(t|b)&(cn<4)
        score=-4.5*distance-.7*cn
    eligible[[0,-1]]=False
    return draw_site(np.where(eligible,score,-np.inf),rng)


def deposit(a,rng,anchors):
    top,pt,s=clusters(a);cn=neighbours(a)
    if not s.any():
        # Effective electrodeposition at cathode-connected residual tips. Source Ag
        # supplied by the active reservoir; solution trajectories are integrated out.
        eligible=~a&(neighbours(pt)>0)
        score=-5.*distance_transform_edt(~top)+.65*cn
    else:
        eligible=~a&(neighbours(s)>0)
        fresh=np.zeros_like(a)
        for y,x in anchors:fresh[y,x]=1
        score=-1.1*distance_transform_edt(~fresh)+1.0*cn
    eligible[[0,-1]]=False
    return draw_site(np.where(eligible,score,-np.inf),rng)


def table(path,rows):
    with path.open('w',encoding='utf8',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)


def positive(initial,v,rng,hazard_limit,hold_add=6):
    a=initial>0
    progress=np.zeros_like(a,float)
    progress[:]=initial
    assert not clusters(a)[2].any()
    # Voltage-activated nucleation followed by an avalanche over a handful of cells
    # within one acquisition dwell; compliance limits post-bridge widening.
    budget=0.;trigger=None;bridge=None;anchors=[];states=[];spans=[];events=[]
    for j,u in enumerate(v):
        if trigger is None:
            budget+=np.expm1(max(0,u)/.032)
            if budget>=hazard_limit:
                trigger=j
                for _ in range(100):
                    y,x=deposit(a,rng,anchors);a[y,x]=1;progress[y,x]=1;anchors.append((y,x))
                    events.append(dict(point=j,event='tip_reduction',y=y,x=x))
                    if clusters(a)[2].any():
                        bridge=(y,x);break
                if bridge is None:raise RuntimeError('Positive switching avalanche failed')
        elif j<=trigger+2*hold_add and j%2==trigger%2 and u>0:
            y,x=deposit(a,rng,anchors);a[y,x]=1;progress[y,x]=1;anchors.append((y,x))
            events.append(dict(point=j,event='compliance_limited_growth',y=y,x=x))
        states.append(progress.copy());spans.append(bool(clusters(a)[2].any()))
    return np.asarray(states,dtype=np.float32),np.array(spans),events,trigger,bridge


def erosion_seq(L1,seed,max_events=110):
    a=L1.copy(); rng=np.random.default_rng(seed)
    seq=[a.copy()]; conductance=[read_network(a)]; events=[]; split=None
    for _ in range(max_events):
        y,x=dissolve(a,rng);a[y,x]=False
        s=clusters(a)[2].any()
        if not s and split is None:split=len(seq)
        raw=read_network(a); conductance.append(raw);seq.append(a.copy())
        events.append((y,x))
        if split is not None and raw<conductance[0]*1.e-8 and len(seq)>split+5:break
    if split is None:raise RuntimeError('No loss of metallic spanning path')
    return np.asarray(seq),np.array(conductance),np.array(events),split


def calibrate_reset(v,i,conductance,split,fixed_parameters=None):
    # Conductance of discrete reaction endpoints from the Ag resistor/tunnel network.
    # For a partially eroded packet, interpolate log conductance (linearly changing
    # tunnelling action). The activity coordinate is proportional to cumulative dose.
    native=np.minimum.accumulate(np.maximum(1.e-22,conductance))
    ratio=native/native[0]
    state=np.arange(len(native))
    train=(v!=0)&(np.arange(len(v))%2==0)
    def calc(p):
        onset,scale,amplitude,high,steep,additional,interface_work,logfloor,tail=p
        drive=np.maximum(-v,0)
        exposure=np.cumsum(.015*np.maximum(drive-onset,0))/scale**2
        activated=.015/steep*np.exp(np.minimum((drive-high)/steep,45))*(-np.expm1(-drive/steep))
        charge=amplitude*np.log1p(exposure)+additional*(-np.expm1(-np.cumsum(activated)))
        # Treat the opening of last contact and subsequent tip erosion as a second
        # effective surface-kinetic regime. One dwell constant per regime, no voltage bins.
        work=np.ones(len(native)-1)
        work[split-1]=interface_work
        work[split:]=tail
        boundaries=np.r_[0,np.cumsum(work)]
        extent=np.interp(charge,boundaries,state)
        logratio=np.interp(extent,state,np.log(ratio)); floor=10**logfloor
        g=floor+(G0-floor)*np.exp(logratio)
        return extent,v*g,g
    if fixed_parameters is not None:
        return *calc(fixed_parameters),fixed_parameters,native
    lo=[0.,.001,.01,.95,.012,.02,.05,-7.5,.015]
    hi=[.17,.9,60,1.7,.3,280,150,-5.2,18]
    def residual(p):
        s,cur,g=calc(p)
        return np.r_[np.log10(abs(cur[train]/i[train])),(cur[train]-i[train])*1100]
    opts=[]
    for scl in [.04,.13,.28]:
        for iw in [1.,7.,30.]:
            start=[.065,scl,3.,1.39,.06,10.,iw,-6.05,.8]
            o=least_squares(residual,start,bounds=(lo,hi),max_nfev=380,x_scale='jac')
            opts.append(o)
    opt=min(opts,key=lambda x:x.cost)
    return *calc(opt.x),opt.x,native


def map_states(seq,extent):
    lo=np.floor(extent).astype(int)
    hi=np.minimum(lo+1,len(seq)-1)
    w=(extent-lo)[:,None,None]
    return (seq[lo]*(1-w)+seq[hi]*w).astype(np.float32)


def logscore(yhat,y,mask):
    return float(np.sqrt(np.mean(np.log10(np.abs(yhat[mask]/y[mask]))**2)))


def verify():
    import hashlib
    with np.load(OUT/'states.npz') as z, np.load(BASE/'forming/selected.npz') as f:
        a=z['metal_fraction']; initial=z['L1']; h1=z['H1']; joins=z['connected_grid']
        np.testing.assert_array_equal(z['H0'],f['H0']);np.testing.assert_array_equal(a[0],f['H0'])
        np.testing.assert_array_equal(a[201],initial);np.testing.assert_array_equal(a[202],initial)
        np.testing.assert_array_equal(a[-1],h1);np.testing.assert_array_equal(f['F'],z['F'])
        np.testing.assert_array_equal(joins, [bool(clusters(im>0)[2].any()) for im in a])
        g0=read_network(initial); v=z['voltage_V'][202:]; model=z['current_A'][202:]
        li={}; packet_coordinates=np.argwhere(np.any(a[202:]!=initial,axis=0))
        def read(im):
            key=im.tobytes()
            if key not in li:li[key]=read_network(im)
            return li[key]
        # Rebuild the two bounding integer states and oxidized fraction from the saved
        # morphology alone. A single site is undergoing the next local oxidation event.
        floor=10**json.loads((OUT/'summary.json').read_text())['kinetic_parameters'][7]
        current=[]
        for bias,im in zip(v,a[202:]):
            partial=im[(im>0)&(im<1)]
            assert partial.size<=1
            gm=read(im>0)
            if partial.size:
                fleft=float(partial[0]);gm=gm**fleft*read(im==1)**(1-fleft)
            current.append(bias*(floor+(G0-floor)*gm/g0))
        np.testing.assert_allclose(current,model,rtol=3e-7,atol=1e-14)
        with np.load(OUT/'repeat_set.npz') as rep:
            np.testing.assert_array_equal(rep['metal_fraction'][0],h1)
            assert rep['connected_grid'][-1] and not rep['connected_grid'][0]
        def slot(a):
            top,bot,both=clusters(a>0)
            if both.any():return 0.
            yu,xu=np.array(np.where(top & (neighbours(top)<4)))
            yl,xl=np.array(np.where(bot & (neighbours(bot)<4)))
            dy=np.maximum(0,abs(yu[:,None]-yl)-1);dx=np.maximum(0,abs(xu[:,None]-xl)-1)
            return .5*np.sqrt((dy**2+dx**2).min())
        stable=[]
        for seed in range(8488,8496):
            rng=np.random.default_rng(seed);replica=initial.copy()
            for order in range(1,51):
                y,x=dissolve(replica,rng);replica[y,x]=0
                if not clusters(replica)[2].any():break
            else:raise ValueError('No disconnect in replica')
            stable.append(dict(rng=seed,local_event=order,z_nm=(y+.5)*.5,x_nm=(x+.5)*.5))
        table(OUT/'choice_sensitivity.csv',stable)
        verification=dict(source_forming_sha256=hashlib.sha256((BASE/'forming/selected.npz').read_bytes()).hexdigest(),
            frames_exact_connectivity=len(a),native_network_relative_current_roundoff=float(np.max(np.abs(np.array(current)-model)/np.maximum(abs(model),1.e-15))),
            H0_min_clearance_nm=float(slot(z['H0'])),H1_min_clearance_nm=float(slot(h1)),
            changed_cells=int(len(packet_coordinates)),changed_row_depths_nm=((np.unique(packet_coordinates[:,0])+.5)*.5).tolist(),
            eight_choice_repeats_break_range_nm=[min(r['z_nm'] for r in stable),max(r['z_nm'] for r in stable)],
            all_electrode_interface_cells_unchanged=bool(np.array_equal(h1[[0,-1]],initial[[0,-1]])))
        (OUT/'verification.json').write_text(json.dumps(verification,indent=2),encoding='utf8')
        params=json.loads((OUT/'summary.json').read_text())['kinetic_parameters']
        (OUT/'reset_parameters.json').write_text(json.dumps(dict(zip(PARAM_NAMES,params)),indent=2),encoding='utf8')
        print(verification)


def main(transfer_from=None):
    OUT.mkdir(parents=True,exist_ok=True)
    with np.load(BASE/'forming/selected.npz') as f:
        h0=f['H0'].astype(bool);first=f['F'].astype(bool)
    with prior.SOURCE.open(encoding='utf-8-sig') as f:base=list(csv.DictReader(f))
    v=np.array([float(r['voltage_V']) for r in base]);measure=np.array([float(r['current_A']) for r in base])
    sv,rv=np.split(v,2);si,ri=np.split(measure,2)
    hazard=float(np.sum(np.expm1(sv[:23]/.032))*.97)
    ss,sc,se,j0,bridge=positive(h0,sv,np.random.default_rng(292),hazard)
    l1=ss[-1]>0;seq,gn,loc,split=erosion_seq(l1,8488)
    progress,reset_i,reset_g,params,mono_g=calibrate_reset(rv,ri,gn,split)
    rs=map_states(seq,progress);rc=np.array([bool(clusters(a>0)[2].any()) for a in rs])
    cut=int(np.flatnonzero(~rc)[0]);physical=loc[split-1]
    def pos_iv(st,connected,initial_read_s=.0):
        g=[]
        lowmask=(sv>0)&(np.arange(len(sv))<j0)
        pre_abs=abs(si[lowmask])/sv[lowmask]
        x=np.c_[np.ones(lowmask.sum()),sv[lowmask]**2]
        intercept,slope=np.linalg.lstsq(x,pre_abs,rcond=None)[0]
        cache={}
        for k,vol in enumerate(sv):
            if not connected[k]:
                g.append((intercept if initial_read_s==0 else initial_read_s)+slope*vol**2)
            else:
                binary=st[k]>0; code=binary.tobytes()
                if code not in cache:cache[code]=read_network(binary)/gn[0]*G0
                g.append(cache[code])
        g=np.array(g); curr=np.minimum(g*sv,.001)
        return curr,curr/g
    set_i,set_vd=pos_iv(ss,sc)
    assert np.array_equal(rs[0],ss[-1])
    res,cc,ee,re,rebridge=positive(rs[-1],sv,np.random.default_rng(293),hazard)
    repeat_i,repeat_v=pos_iv(res,cc,reset_g[-1])
    np.savez_compressed(OUT/'repeat_set.npz',metal_fraction=res,voltage_V=sv,device_voltage_V=repeat_v,
                        current_A=repeat_i,connected_grid=cc)
    phase=np.repeat(['SET','RESET'],[202,202]); allstates=np.concatenate([ss,rs]);curr=np.r_[set_i,reset_i]
    used=int(np.ceil(progress[-1]));neck=loc[:used]
    np.savez_compressed(OUT/'states.npz',metal_fraction=allstates,voltage_V=v,device_voltage_V=np.r_[set_vd,rv],
                        current_A=curr,measured_current_A=measure,phase=phase,H0=h0,L1=l1,H1=rs[-1],F=first,
                        connected_grid=np.r_[sc,rc],contact_yx=neck,rupture_yx=physical,SET_join_yx=bridge,
                        RESET_extent=progress,native_G_S=mono_g,dissolution_order_yx=loc,reset_G_S=reset_g)
    table(OUT/'set_events.csv',se);table(OUT/'repeat_set_events.csv',ee)
    table(OUT/'reset_events.csv',[dict(order=k+1,y=int(y),x=int(x),depth_nm=(y+.5)*.5,
                                  reached_within_sweep=bool(k<used),first_bridge_loss=k+1==split,
                                  native_conductance_after_S=gn[k+1]) for k,(y,x) in enumerate(loc)])
    table(OUT/'iv.csv',[dict(phase=phase[k],point=k%202,voltage_V=vv,measured_A=measure[k],model_A=curr[k],
                connected_grid=int(np.r_[sc,rc][k]),Ag_cells=float(allstates[k].sum())) for k,vv in enumerate(v)])
    stats=dict(set_onset_V=sv[j0],set_join_yx=bridge,set_join_depth_nm=(bridge[0]+.5)*.5,reset_first_break_index=cut,
               reset_first_break_V=rv[cut],reset_break_yx=physical.tolist(),reset_break_depth_nm=(physical[0]+.5)*.5,
               reSET_onset_V=sv[re],reSET_join_yx=rebridge,reSET_depth_nm=(rebridge[0]+.5)*.5,
               SET_logRMSE_dec=logscore(set_i,si,sv!=0),RESET_logRMSE_dec=logscore(reset_i,ri,rv!=0),
               RESET_interleaved_heldout_logRMSE_dec=logscore(reset_i,ri,(rv!=0)&(np.arange(202)%2==1)),
               RESET_out_rmse_dec=logscore(reset_i,ri,(rv!=0)&(np.arange(202)<101)),
               RESET_back_rmse_dec=logscore(reset_i,ri,(rv!=0)&(np.arange(202)>100)),
               H1_R_ohm=1/reset_g[-1],H1_removed_cells=float(l1.sum()-rs[-1].sum()),sites_min_max_depth_nm=[float((neck[:,0].min()+.5)*.5),float((neck[:,0].max()+.5)*.5)],
               distinct_eroding_rows=int(np.unique(neck[:,0]).size),split_order=split,
               kinetic_parameters=params.tolist(),peak_current_out_V=rv[np.argmax(abs(reset_i[:101]))],
               native_increases=int(np.sum(np.diff(gn)>max(gn)*1e-8)),initial_seed=292,reset_choice_seed=8488)
    (OUT/'summary.json').write_text(json.dumps(stats,indent=2,default=lambda z:z.item()),encoding='utf8')
    print(json.dumps(stats,default=lambda z:z.item()))
    if transfer_from is not None:
        previous=json.loads((transfer_from/'cycle/summary.json').read_text(encoding='utf8'))
        extent,current,g,_,_=calibrate_reset(rv,ri,gn,split,previous['kinetic_parameters'])
        table(OUT/'fixed_parameters_iv.csv',[dict(point=k,voltage_V=rv[k],measured_A=ri[k],
              transferred_A=current[k],refitted_A=reset_i[k],transferred_extent=extent[k])
              for k in range(len(rv))])
        breaks=np.flatnonzero(extent>=split)
        score=dict(source=str(transfer_from),reset_degree_of_freedom_refitted=0,
                   reset_logRMSE_dec=logscore(current,ri,rv!=0),
                   H1_R_ohm=1/g[-1],split_order=split,
                   final_extent=float(extent[-1]),
                   first_break_V=float(rv[breaks[0]]) if len(breaks) else None)
        (OUT/'fixed_parameters_summary.json').write_text(json.dumps(score,indent=2),encoding='utf8')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,default=BASE,
                        help='Directory containing forming/selected.npz; write into its cycle subdirectory.')
    parser.add_argument('--verify-only',action='store_true')
    parser.add_argument('--transfer-from',type=Path,
                        help='Also evaluate RESET with all kinetic parameters fixed to another run.')
    args=parser.parse_args()
    BASE=args.run_dir.resolve()
    OUT=BASE/'cycle'
    if not args.verify_only:main(args.transfer_from)
    verify()
