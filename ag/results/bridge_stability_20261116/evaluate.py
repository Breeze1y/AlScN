"""Read-only comparison of saved full loops; never resumes legacy physics."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))
from run_kmc import hashes, read_forming, write_json
from forming_iv import ElectronicConfig, circuit_point

def read_checkpoint(folder):
    keys=('state','site','vacancy_mask','scan_states','first_bridge_state',
          'hop_layer_counts','net_deposition','ion_residence_s',
          'vacancy_ion_residence_s','model_snapshots')
    with np.load(folder/'checkpoint.npz',allow_pickle=False) as z:
        text=str(z['metadata']);meta=json.loads(text)
        arrays={k:z[k].copy() for k in keys}
        digest=hashlib.sha256(text.encode()+b''.join(v.tobytes() for v in arrays.values())).hexdigest()
        assert digest==str(z['checksum']), 'Checkpoint checksum mismatch'
    current=hashes()
    old={name:hashlib.sha256((HERE/'baseline_sources'/name).read_bytes()).hexdigest() for name in current}
    assert meta['physics_hashes'] in (current,old), 'Unrecognized physics revision'
    return meta,arrays,'corrected' if meta['physics_hashes']==current else 'baseline'

def assess(folder):
    folder=Path(folder).resolve()
    result=json.loads((folder/'result.json').read_text(encoding='utf-8'))
    assert result['status']=='loop_and_zero_hold_complete'
    meta,arrays,revision=read_checkpoint(folder)
    wave=[r for r in meta['records'] if r['phase']=='waveform']
    base=json.loads((HERE/'baseline_settings.json').read_text(encoding='utf-8'))
    data=read_forming(ROOT/base['protocol']['source'],'full')
    assert len(wave)==202 and [r['point_index'] for r in wave]==list(range(202))
    v=np.array([r['command_voltage_V'] for r in wave])
    np.testing.assert_allclose(v,[p['voltage_V'] for p in data['points']],rtol=0,atol=0)
    np.testing.assert_allclose([r['time_s'] for r in wave],np.arange(1,203)*.02,rtol=0,atol=1e-12)
    assert np.count_nonzero(arrays['scan_states'][0,1:-1])==0
    assert result['ag_balance_error']==result['reaction_balance_error']==0
    assert meta['electronics']==base['electronics']
    allowed={'attempt_Hz','bond_eV','source_barrier_eV','seed','voltage','max_time_s','max_events'}
    assert all(meta['config'][k]==value for k,value in base['model'].items() if k not in allowed)
    assert sum(result['counts'].values())==result['events']
    assert meta['events']==result['events']
    final=arrays['state'].reshape(arrays['scan_states'].shape[1:])
    assert np.all(final[[0,-1]]==2)
    assert int(np.count_nonzero(final[1:-1]))==result['net_source']
    assert int(np.sum(final[1:-1]==2))==result['counts']['reduce']-result['counts']['oxidize']

    # Between metal events G is constant; at each voltage boundary recalculate I.
    changes=[r for r in meta['iv'] if r['kind']=='metal_change']
    cursor=0;start=0.;gprev=meta['iv'][0]['network_conductance_S'];means=[]
    electronics=ElectronicConfig(**meta['electronics'])
    for row in wave:
        end=row['time_s'];t=start;charge=0.
        current=circuit_point(row['command_voltage_V'],gprev,meta['config']['nx'],electronics)['current_A']
        while cursor<len(changes) and changes[cursor]['time_s']<=end:
            event=changes[cursor];cursor+=1
            if event['time_s']<start:continue
            charge+=current*(event['time_s']-t);t=event['time_s'];current=event['current_A']
        charge+=current*(end-t);means.append(charge/(end-start))
        start=end;gprev=row['network_conductance_S']
    sampled=np.array([r['current_A'] for r in wave]);means=np.asarray(means)
    def first_on(values,threshold=.0005):
        indices=np.flatnonzero(np.abs(values[:101])>=threshold)
        return int(indices[0]) if len(indices) else None
    idx=first_on(sampled);avgidx=first_on(means)
    in_window=lambda value:value is not None and 1.320-1e-12<=value<=1.335+1e-12
    drop=next((i for i in range(idx+1,101) if abs(sampled[i])<.0005),None) if idx is not None else None
    last_time=meta['time_s'];start=last_time-.1
    continuous=[r for r in meta['iv'] if start-1e-12<=r['time_s']<=last_time+1e-12]
    before=[r for r in meta['iv'] if r['time_s']<start-1e-12]
    if before:continuous.insert(0,dict(before[-1],time_s=start))
    g=np.array([r['network_conductance_S'] for r in continuous]);metal=np.array([r['metal'] for r in continuous])
    grange=float(np.ptp(g)/max(np.median(g),1e-30));mrange=float(np.ptp(metal)/max(np.median(metal),1))
    zero_start=wave[-1]['time_s']
    ruptures=[r for r in meta['topology_log'] if not r['connected']]
    # Identifies event class from metal-count change, without claiming a site location.
    iv_by_event={r['events']:i for i,r in enumerate(meta['iv']) if r['kind']=='metal_change'}
    causes={}
    for r in ruptures:
        i=iv_by_event[r['events']]
        delta=meta['iv'][i]['metal']-meta['iv'][i-1]['metal']
        kind={-1:'oxidize',0:'migrate',1:'reduce'}.get(delta,'unknown')
        causes[kind]=causes.get(kind,0)+1
    first=result['first_bridge']
    stats=dict(folder=str(folder),solver=revision,seed=meta['config']['seed'],
        attempt_Hz=meta['config']['attempt_Hz'],bond_eV=meta['config']['bond_eV'],
        source_barrier_eV=meta['config']['source_barrier_eV'],
        full_protocol_complete=True,source_sha256=data['source_sha256'],first_bridge=first,
        first_endpoint_opening_V=float(v[idx]) if idx is not None else None,
        first_whole_dwell_mean_opening_V=float(v[avgidx]) if avgidx is not None else None,
        target_window_V=[1.32,1.335],
        inclusive_window_pass=bool(first and in_window(first['command_voltage_V']) and idx is not None and in_window(v[idx]) and avgidx is not None and in_window(v[avgidx])),
        exact_experimental_bin_pass=bool(idx==89 and abs(sampled[88])<.0005 and abs(sampled[90])>=.0005),
        first_up_dropout_V=float(v[drop]) if drop is not None else None,
        endpoint_up_plateau_retained=idx is not None and drop is None,
        geometric_ruptures=len(ruptures),rupture_event_classes=causes,
        up_ruptures=sum(r['time_s']<=wave[100]['time_s']+1e-12 for r in ruptures),
        zero_ruptures=sum(r['time_s']>zero_start+1e-12 for r in ruptures),
        zero_hold_s=result['zero_hold_s'],final_connected=result['connected'],
        final_metal=result['endpoint']['metal'],final_ions=result['endpoint']['ions'],
        vertex_bottleneck_4=result['morphology']['vertex_bottleneck_4'],
        continuous_last100ms=dict(records=len(continuous),conductance_relative_range=grange,
            metal_relative_range=mrange,all_connected=all(r['connected'] for r in continuous)),
        continuous_filament_stable=bool(all(r['connected'] for r in continuous) and grange<=.05 and mrange<=.01),
        endpoint_stability=result['stability'],ag_balance_error=0,reaction_balance_error=0,
        point_1p320=dict(endpoint_A=sampled[88],whole_dwell_mean_A=means[88]),
        point_1p335=dict(endpoint_A=sampled[89],whole_dwell_mean_A=means[89]),
        electrical_threshold_A=.0005,
        stability_definition='Last 100 ms, all metal events and endpoint observations plus interval-start state: geometrically connected; G range / record median <=5%; metal range / record median <=1%. Short observation, not long-term retention.',
        averaging_scope='Average over full nominal 20 ms, without instrument settling/integration timing.',
        physics_hashes=meta['physics_hashes'])
    write_json(folder/'stability_assessment.json',stats)
    write_json(folder/'comparison_series.json',dict(voltage_V=v.tolist(),endpoint_A=sampled.tolist(),
        whole_dwell_mean_A=means.tolist(),measured_A=[p['measured_A'] for p in data['points']]))
    return stats,meta,arrays

def plot_comparison(folders):
    cases=[assess(f) for f in folders]
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    colors=['#777777','#b76526','#147e95','#6b5595']
    data=json.loads((Path(folders[-1])/'comparison_series.json').read_text())
    v=np.array(data['voltage_V'])
    for ax in axes[0]:
        ax.semilogy(v[:101],np.maximum(np.abs(data['measured_A'][:101]),1e-14),'k--',lw=1,label='Measured up')
        ax.axvspan(1.32,1.335,color='#76ad63',alpha=.18)
        ax.axhline(.0005,color='#999999',lw=.7,ls=':')
    for (a,m,arrays),color in zip(cases,colors):
        series=json.loads((Path(a['folder'])/'comparison_series.json').read_text())
        label=f"{a['solver']}; J={a['bond_eV']:.2f}; Bsrc={a['source_barrier_eV']:.3f}; nu={a['attempt_Hz']/1e12:g}e12; s{a['seed']}"
        for ax in axes[0]:ax.semilogy(v[:101],np.maximum(np.abs(series['endpoint_A'][:101]),1e-14),color=color,lw=1.3,label=label)
        t=np.array([r['time_s'] for r in m['iv']]);g=np.array([r['network_conductance_S'] for r in m['iv']])
        axes[1,0].step(t,g*1e3,where='post',color=color,lw=1,label=label)
        mask=t>=4.04-1e-12
        axes[1,1].step((t[mask]-4.04)*1e3,g[mask]*1e3,where='post',color=color,lw=1)
    axes[0,0].set(title='Up-sweep endpoint current',xlim=(0,1.5),xlabel='Command voltage (V)',ylabel='|I| (A)')
    axes[0,1].set(title='Opening window',xlim=(1.26,1.5),ylim=(1e-8,.002),xlabel='Command voltage (V)',ylabel='|I| (A)')
    axes[1,0].set(title='Event-resolved network conductance',xlabel='Time (s)',ylabel='Network G (mS)')
    axes[1,1].set(title='Zero-bias observation after loop',xlabel='Time after return to 0 V (ms)',ylabel='Network G (mS)')
    for ax in axes.flat:ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7,loc='lower right')
    fig.savefig(HERE/'stability_comparison.png',dpi=180);plt.close(fig)

    a,m,arrays=cases[-1];labels=m['state_labels'];s=arrays['scan_states']
    states=[arrays['first_bridge_state'],s[labels.index('waveform_201')],arrays['state'].reshape(s.shape[1:])]
    titles=[f"First bridge: {a['first_bridge']['command_voltage_V']:.3f} V",'Return to 0 V',f"Zero bias: {a['zero_hold_s']*1000:.0f} ms"]
    fig,axes=plt.subplots(3,1,figsize=(10,10),layout='constrained')
    cmap=ListedColormap(['#edf2fa','#2eb1a5','#d89319'])
    for ax,state,title in zip(axes,states,titles):
        ax.imshow(state,cmap=cmap,vmin=0,vmax=2,extent=(0,150,61,-1),interpolation='none')
        for z in (20,40):ax.axhline(z,color='#596586',ls='--',lw=.7)
        for low,high in [(-1,0),(60,61)]:ax.axhspan(low,high,color='#8797a5')
        ax.set(xlabel='x (nm)',ylabel='Depth from Ag (nm)',title=f"{title} | Ag {np.sum(state[1:-1]==2)}; Ag+ {np.sum(state[1:-1]==1)}")
    fig.legend(handles=[Patch(color=c,label=t) for c,t in [('#d89319','Metal Ag'),('#2eb1a5','Ag+'),('#8797a5','Fixed electrodes')]],loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle(f"Fresh forming, seed {a['seed']}, J={a['bond_eV']:.2f} eV, Bsrc={a['source_barrier_eV']:.3f} eV, nu={a['attempt_Hz']/1e12:g}e12 Hz")
    fig.savefig(HERE/'selected_morphology.png',dpi=180);plt.close(fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--folder',type=Path,nargs='+',required=True);parser.add_argument('--plot',action='store_true')
    args=parser.parse_args()
    for folder in args.folder:
        a,_,_=assess(folder)
        print(json.dumps({k:a[k] for k in ('folder','first_endpoint_opening_V','inclusive_window_pass','exact_experimental_bin_pass','endpoint_up_plateau_retained','geometric_ruptures','continuous_filament_stable','continuous_last100ms')},ensure_ascii=False))
    if args.plot:plot_comparison(args.folder)
