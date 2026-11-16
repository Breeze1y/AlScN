"""Current entry point for the restored infinite-Ag forming baseline."""
import argparse
from dataclasses import asdict, replace
import csv
import hashlib
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np

from kmc_model import Config
from forming_iv import ElectronicConfig, run_forming
from fast_kmc import CachedFormingModel
from morphology_metrics import metrics

ROOT=Path(__file__).resolve().parent
PHYSICS_FILES=('kmc_model.py','forming_iv.py','fast_kmc.py','morphology_metrics.py','dynamic_kmc.py')


def hashes():
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in PHYSICS_FILES}


def write_json(path,data):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temporary.replace(path)


def read_forming(path,sweep='up'):
    points=[]; columns=None; compliance=None
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:
        for line,row in enumerate(csv.reader(handle),1):
            row=[value.strip() for value in row]
            if not row:continue
            if row[0]=='TestParameter' and row[1]=='Measurement.Primary.Compliance':compliance=float(row[2])
            elif row[0]=='DataName':columns=row[1:]
            elif row[0]=='DataValue':
                if columns is None:raise ValueError('Missing measured data columns')
                values=dict(zip(columns,row[1:]))
                points.append(dict(voltage_V=float(values['V1']),measured_A=float(values['I1']),source_row=line))
    data=np.array([[p['voltage_V'],p['measured_A']] for p in points])
    if len(points)<3 or not np.isfinite(data).all() or compliance is None or compliance<=0:
        raise ValueError('Select a finite B1500 forming export with positive compliance')
    peak=int(np.argmax(data[:,0]))
    if peak==0 or np.any(np.diff(data[:peak+1,0])<0):raise ValueError('Expected positive forming up sweep')
    return dict(source=str(Path(path).resolve()),source_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                original_points=len(points),points=points if sweep=='full' else points[:peak+1],compliance_A=compliance)


def save_dynamic(folder,m,position,records,states,labels,status,audit):
    metadata=dict(schema='infinite_ag_dynamic_v1',physics_hashes=hashes(),config=asdict(m.c),electronics=asdict(m.electronics),
                  command_voltage=m.command_voltage,time_s=m.time,events=m.steps,net_source=m.net_source,counts=m.counts,
                  rng=m.rng.bit_generator.state,pending_event=m.pending_event,first_bridge=m.first_bridge,
                  topology_log=m.topology_log,position=position,records=records,state_labels=labels,
                  snapshot_target=m.snapshot_target,trace=m.trace,iv=m.iv,event_log=m.event_log,electronic_updates=m.electronic_updates)
    arrays=dict(state=m.state,site=m.site,vacancy_mask=m.vacancy_mask,scan_states=np.asarray(states),
                first_bridge_state=m.first_bridge_state if m.first_bridge_state is not None else np.array([],np.uint8),
                hop_layer_counts=m.hop_layer_counts,net_deposition=m.net_deposition,ion_residence_s=m.ion_residence_s,
                vacancy_ion_residence_s=m.vacancy_ion_residence_s,model_snapshots=np.asarray(m.snapshots))
    text=json.dumps(metadata,sort_keys=True,allow_nan=False)
    checksum=hashlib.sha256(text.encode()+b''.join(v.tobytes() for v in arrays.values())).hexdigest()
    temporary=folder/'checkpoint.npz.tmp'
    with temporary.open('wb') as h:np.savez_compressed(h,metadata=text,checksum=checksum,**arrays)
    temporary.replace(folder/'checkpoint.npz')
    write_json(folder/'progress.json',dict(status=status,events=m.steps,time_s=m.time,connected=m.connected,
               first_bridge=m.first_bridge,position=position,checkpoint='checkpoint.npz',endpoint=m.iv[-1]))
    write_json(folder/'manifest.json',audit)


def load_dynamic(folder):
    from dynamic_kmc import DynamicCachedModel
    with np.load(folder/'checkpoint.npz',allow_pickle=False) as a:
        text=str(a['metadata']);meta=json.loads(text)
        keys=('state','site','vacancy_mask','scan_states','first_bridge_state','hop_layer_counts','net_deposition','ion_residence_s','vacancy_ion_residence_s','model_snapshots')
        arrays={k:a[k].copy() for k in keys}
        if hashlib.sha256(text.encode()+b''.join(v.tobytes() for v in arrays.values())).hexdigest()!=str(a['checksum']):raise ValueError('Checkpoint checksum mismatch')
    if meta['schema']!='infinite_ag_dynamic_v1' or meta['physics_hashes']!=hashes():raise ValueError('Dynamic checkpoint version differs')
    m=DynamicCachedModel(Config(**meta['config']),ElectronicConfig(**meta['electronics']))
    np.testing.assert_array_equal(m.site,arrays['site']);np.testing.assert_array_equal(m.vacancy_mask,arrays['vacancy_mask'])
    m.state=arrays['state'];m.time=meta['time_s'];m.steps=meta['events'];m.net_source=meta['net_source'];m.counts=meta['counts']
    m.command_voltage=meta['command_voltage'];m.first_bridge=meta['first_bridge']
    m.first_bridge_state=arrays['first_bridge_state'] if arrays['first_bridge_state'].size else None
    m.update_structure()
    m.topology_log=meta['topology_log'];m.iv=meta['iv'];m.trace=meta['trace'];m.snapshot_target=meta['snapshot_target']
    m.event_log=meta['event_log'];m.electronic_updates=meta['electronic_updates'];m.snapshots=list(arrays['model_snapshots'])
    for key in ('hop_layer_counts','net_deposition','ion_residence_s','vacancy_ion_residence_s'):setattr(m,key,arrays[key])
    m.pending_event=meta['pending_event'];m.rng.bit_generator.state=meta['rng']
    return m,meta['position'],meta['records'],list(arrays['scan_states']),meta['state_labels']


def run_dynamic(args,settings,folder,audit,resume=False):
    from dynamic_kmc import DynamicCachedModel,stability
    measurement=read_forming(ROOT/settings['protocol']['source'],'full');dwell=settings['protocol']['dwell_s']
    settling=settings['protocol'].get('settling',dict(min_s=.1,max_s=.2,window_s=.02))
    if resume:
        m,pos,records,states,labels=load_dynamic(folder)
    else:
        cfg=replace(Config(**settings['model']),max_events=args.max_events or 20000000,
                    max_time_s=len(measurement['points'])*dwell+settling['max_s']+1e-12)
        e=replace(ElectronicConfig(**settings['electronics']),compliance_A=measurement['compliance_A'])
        m=DynamicCachedModel(cfg,e)
        pos=dict(index=0,deadline_s=None,phase='waveform',zero_start_s=None)
        records=[];states=[m.state.reshape(m.shape).copy()];labels=['initial']
    audit.update(measurement=measurement,dwell_s=dwell,post_bridge_mode='dynamic',settling=settling)
    limit=m.steps+(args.max_events or 20000000);wall=time.perf_counter();saved=m.steps
    status='running';assessment=dict(assessed=False,filament_stable=False,ions_stationary=False)
    while True:
        phase=pos['phase'];i=pos['index']
        if phase=='settled':status='loop_and_zero_hold_complete';break
        if pos['deadline_s'] is None:
            voltage=measurement['points'][i]['voltage_V'] if phase=='waveform' else 0.
            duration=dwell if phase=='waveform' else settling['window_s']
            m.apply_command(voltage);pos['deadline_s']=m.time+duration
        result=m.advance(pos['deadline_s'],min(limit,m.steps+1000))
        if result['status']=='point_complete':
            m.record('waveform_end' if phase=='waveform' else 'zero_hold_end')
            row=dict(phase=phase,point_index=i,**m.iv[-1]);records.append(row)
            states.append(m.state.reshape(m.shape).copy());labels.append(f'{phase}_{i:03d}')
            print('stage',phase,'point',i,'V',m.command_voltage,'events',m.steps,'metal',row['metal'],'ions',row['ions'],'connected',m.connected,flush=True)
            pos['deadline_s']=None;pos['index']+=1
            if phase=='waveform' and pos['index']==len(measurement['points']):
                pos.update(phase='zero_hold',index=0,zero_start_s=m.time)
            elif phase=='zero_hold':
                zero=[r for r in records if r['phase']=='zero_hold' or (r['phase']=='waveform' and r['point_index']==len(measurement['points'])-1)]
                assessment=stability(zero,m.topology_log)
                held=m.time-pos['zero_start_s']
                if held>=settling['max_s']-1e-12 or (held>=settling['min_s']-1e-12 and assessment['filament_stable'] and assessment['ions_stationary']):
                    pos['phase']='settled';status='loop_and_zero_hold_complete';break
        if m.steps>=limit:status='event_budget_unknown';break
        if time.perf_counter()-wall>=args.max_wall_s:status='wall_budget_unknown';break
        if m.steps-saved>=50000:
            save_dynamic(folder,m,pos,records,states,labels,'running_checkpoint',audit);saved=m.steps
    zero=[r for r in records if r['phase']=='zero_hold' or (r['phase']=='waveform' and r['point_index']==len(measurement['points'])-1)]
    assessment=stability(zero,m.topology_log)
    save_dynamic(folder,m,pos,records,states,labels,status,audit)
    summary=dict(status=status,fresh_forming_run=True,connected=m.connected,first_bridge=m.first_bridge,
                 events=m.steps,time_s=m.time,counts=m.counts,net_source=m.net_source,
                 ag_balance_error=int(np.count_nonzero(m.state[m.interior]))-m.net_source,
                 reaction_balance_error=int(np.sum(m.state[m.interior]==2))-(m.counts['reduce']-m.counts['oxidize']),
                 completed_points=sum(r['phase']=='waveform' for r in records),total_points=len(measurement['points']),
                 actual_end_voltage_V=m.command_voltage,post_bridge_mode='dynamic',frozen_shape=False,
                 zero_hold_s=m.time-pos['zero_start_s'] if pos['zero_start_s'] is not None else 0.,
                 stability=assessment,endpoint=m.iv[-1],morphology=metrics(m.state.reshape(m.shape),m.c.spacing_nm,m.c.layer_cells),
                 topology_log=m.topology_log,config=asdict(m.c),electronics=asdict(m.electronics))
    write_json(folder/'result.json',summary)
    np.savez_compressed(folder/'states.npz',state=m.state.reshape(m.shape),site_energy_eV=m.site.reshape(m.shape),
                        vacancy_mask=m.vacancy_mask.reshape(m.shape),scan_states=np.asarray(states),labels=labels,
                        first_bridge_state=m.first_bridge_state if m.first_bridge_state is not None else np.array([],np.uint8))
    if records:
        with (folder/'forming_trace.csv').open('w',encoding='utf-8-sig',newline='') as h:
            writer=csv.DictWriter(h,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    draw_morphology(folder,m)
    draw_dynamic_results(folder,m,records,states,labels)
    return summary


def draw_dynamic_results(folder,m,records,states,labels):
    n=m.c.layer_cells;shape=m.shape;final=m.state.reshape(shape)
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
    cmap=ListedColormap(['#f4dfbd','#d3d5af','#efdab5','#287b93','#263d5e','#f08731'])
    for ax,ions,title in zip(axes,(True,False),('Final state: Ag metal and Ag+','Final Ag metal only (ions hidden)')):
        pixels=np.zeros(shape,int);pixels[1+n:1+2*n]=1;pixels[1+2*n:-1]=2
        pixels[final==2]=3;pixels[[0,-1]]=4
        if ions:pixels[final==1]=5
        ax.imshow(pixels,cmap=cmap,vmin=0,vmax=5,extent=(0,m.c.nx*m.c.spacing_nm,(3*n+1)*m.c.spacing_nm,-m.c.spacing_nm),interpolation='none')
        for d in (0,n,2*n,3*n):ax.axhline(d*m.c.spacing_nm,color='#777',ls='--',lw=.5)
        ax.set(xlabel='x (nm)',ylabel='Depth (nm)',title=title)
    fig.suptitle(f'After executed waveform/hold: command {m.command_voltage:g} V; no ions removed from simulation')
    fig.savefig(folder/'final_filament.png',dpi=180);plt.close(fig)
    if records:
        fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
        t=[r['time_s'] for r in records]
        axes[0,0].plot(t,[r['command_voltage_V'] for r in records],label='Command');axes[0,0].plot(t,[r['device_voltage_V'] for r in records],label='Device');axes[0,0].set(ylabel='Voltage (V)');axes[0,0].legend()
        axes[0,1].plot(t,[r['metal'] for r in records],label='Metal Ag');axes[0,1].plot(t,[r['ions'] for r in records],label='Ag+');axes[0,1].set(ylabel='Effective packets');axes[0,1].legend()
        axes[1,0].semilogy(t,np.maximum([r['network_conductance_S'] for r in records],1e-16));axes[1,0].set(ylabel='Network conductance (S)')
        wave=[r for r in records if r['phase']=='waveform']
        axes[1,1].semilogy([r['command_voltage_V'] for r in wave],np.maximum([abs(r['current_A']) for r in wave],1e-15));axes[1,1].set(xlabel='Command voltage (V)',ylabel='Current (A)')
        for ax in axes.ravel():ax.grid(alpha=.2)
        for ax in (axes[0,0],axes[0,1],axes[1,0]):ax.set(xlabel='Model time (s)')
        fig.savefig(folder/'loop_evolution.png',dpi=180);plt.close(fig)


def restore_filament(path):
    """Recover a historical shape and frozen readout, not an RNG continuation."""
    path=Path(path).resolve(); digest=hashlib.sha256(path.read_bytes()).hexdigest()
    sidecar=path.parent/'manifest.json'
    if sidecar.exists():
        archived=json.loads(sidecar.read_text(encoding='utf-8'))
        for row in archived.get('states',[]):
            if row['file']==path.name and row['sha256']!=digest:raise ValueError('Historical filament checksum mismatch')
    with np.load(path,allow_pickle=False) as archive:
        required={'state','site_energy_eV','vacancy_mask','metadata'}
        if not required.issubset(archive.files):raise ValueError('Expected an archived infinite-source forming shape')
        state=archive['state'].copy(); site=archive['site_energy_eV'].copy(); vacancy=archive['vacancy_mask'].copy()
        metadata=json.loads(str(archive['metadata']))
    cfg=Config(**metadata['forming_config']); e=ElectronicConfig(**metadata['electronics'])
    m=CachedFormingModel(cfg,e)
    if state.shape!=m.shape or not np.isin(state,(0,1,2)).all() or not np.all(state[[0,-1]]==2):
        raise ValueError('Historical state has incompatible geometry or fixed electrodes')
    np.testing.assert_array_equal(site.ravel(),m.site)
    np.testing.assert_array_equal(vacancy.ravel(),m.vacancy_mask)
    m.state=state.ravel().copy();m.time=metadata['forming_time_s'];m.steps=metadata['forming_events']
    m.net_source=metadata['net_source'];m.counts=metadata['counts'].copy()
    inside=m.state[m.interior]
    if (np.count_nonzero(inside)!=m.net_source or m.net_source!=m.counts['inject']-m.counts['absorb']
            or np.count_nonzero(inside==2)!=m.counts['reduce']-m.counts['oxidize']):
        raise ValueError('Historical source/reaction inventory mismatch')
    m.command_voltage=metadata['forming_readout']['command_voltage_V']
    m.first_bridge=None;m.update_structure()
    if not m.connected:raise ValueError('Archived first-bridge shape is not connected')
    np.testing.assert_allclose(m.operating['current_A'],metadata['forming_readout']['current_A'],rtol=1e-11,atol=1e-14)
    np.testing.assert_allclose(m.operating['device_voltage_V'],metadata['forming_readout']['device_voltage_V'],rtol=1e-11,atol=1e-13)
    m.pending_event=None
    # No historical event RNG or post-bridge field is available in this archive.
    # This object is used only for shape diagnostics and its frozen network.
    audit=dict(source=str(path),source_sha256=digest,exact_state_match=bool(np.array_equal(m.state.reshape(m.shape),state)),
               mode='historical_shape_restoration_not_fresh_simulation',historical_dwell_s=metadata.get('protocol_dwell_s',.01),
               rng_continuation_available=False,historical_metadata=metadata)
    return m,audit


def draw_morphology(folder,m,restored=False):
    state=m.state.reshape(m.shape)
    blank=np.zeros(m.shape,np.uint8);blank[[0,-1]]=2
    shown=[(blank,'Blank dielectric + fixed Ag/Pt'),(state,'Restored historical filament' if restored else 'Actual final morphology')]
    cmap=ListedColormap(['#f4dfbd','#d3d5af','#efdab5','#287b93','#263d5e','#f08731'])
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
    for ax,(snapshot,title) in zip(axes,shown):
        pixels=np.zeros(m.shape,int)
        n=m.c.layer_cells;pixels[1+n:1+2*n]=1;pixels[1+2*n:-1]=2
        pixels[snapshot==2]=3;pixels[[0,-1]]=4;pixels[snapshot==1]=5
        ax.imshow(pixels,cmap=cmap,vmin=0,vmax=5,extent=(0,m.c.nx*m.c.spacing_nm,
                  (3*n+1)*m.c.spacing_nm,-m.c.spacing_nm),interpolation='none')
        for d in (0,n,2*n,3*n):ax.axhline(d*m.c.spacing_nm,color='#777',ls='--',lw=.5)
        ax.set(xlabel='x (nm)',ylabel='Depth from Ag interface (nm)',title=title)
    fig.legend([Patch(facecolor=c) for c in cmap.colors],['AlN top','AlScN','AlN bottom','Ag metal','Fixed electrode','Ag+'],
               loc='outside lower center',ncol=3,frameon=False,fontsize=9)
    fig.suptitle('Infinite Ag reservoir: fixed one-row boundary; no dissolving electrode volume')
    fig.savefig(folder/'morphology.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=ROOT/'configs/current_kmc.json')
    parser.add_argument('--resume',type=Path,help='Resume a dynamic loop in the same result directory')
    parser.add_argument('--restore-filament',nargs='?',const='reference',help='Recover saved historical shape, without rerunning forming')
    parser.add_argument('--max-events',type=int,default=None)
    parser.add_argument('--max-wall-s',type=float,default=300,help='Checked at completed voltage segments')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.max_wall_s<=0 or (args.max_events is not None and args.max_events<=0):parser.error('Positive budgets required')
    if args.resume:
        if args.output or args.restore_filament:parser.error('Resume uses its saved output and dynamic protocol')
        folder=args.resume.resolve();audit=json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
        summary=run_dynamic(args,audit['configuration'],folder,audit,resume=True)
        print(json.dumps({k:summary[k] for k in ('status','connected','events','actual_end_voltage_V','stability')},indent=2));return
    if args.output is None:parser.error('Select a new output directory')
    settings=json.loads(args.config.read_text(encoding='utf-8'))
    if settings.get('backend')!='infinite_ag_v1_restored':parser.error('Select the restored infinite-source configuration')
    folder=args.output.resolve();folder.mkdir(parents=True,exist_ok=False)
    audit=dict(physics_source_sha256=hashes(),runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),configuration=settings)
    if args.restore_filament is not None:
        path=ROOT/settings['restore_reference'] if args.restore_filament=='reference' else Path(args.restore_filament)
        m,restoration=restore_filament(path)
        summary=dict(status='historical_filament_restored',fresh_forming_run=False,connected=m.connected,
                     metal=int(np.sum(m.state[m.interior]==2)),ions=int(np.sum(m.state==1)),net_source=m.net_source,
                     counts=m.counts,ag_balance_error=int(np.count_nonzero(m.state[m.interior]))-m.net_source,
                     original_first_bridge=m.first_bridge,frozen_readout=m.operating,
                     morphology=metrics(m.state.reshape(m.shape),m.c.spacing_nm,m.c.layer_cells),**restoration)
        np.savez_compressed(folder/'restored_filament.npz',state=m.state.reshape(m.shape),site_energy_eV=m.site.reshape(m.shape),
                            vacancy_mask=m.vacancy_mask.reshape(m.shape),metadata=json.dumps(restoration['historical_metadata']))
        draw_morphology(folder,m,restored=True)
    elif settings['protocol'].get('evolution')=='dynamic':
        summary=run_dynamic(args,settings,folder,audit)
    else:
        measurement=read_forming(ROOT/settings['protocol']['source'])
        dwell=settings['protocol']['dwell_s']
        if not np.isfinite(dwell) or dwell<=0:raise ValueError('Invalid point duration')
        protocol=[dict(voltage_V=p['voltage_V'],duration_s=dwell) for p in measurement['points']]
        cfg=Config(**settings['model']);cfg=replace(cfg,max_time_s=sum(x['duration_s'] for x in protocol)+1e-12,
                      max_events=args.max_events or cfg.max_events)
        e=replace(ElectronicConfig(**settings['electronics']),compliance_A=measurement['compliance_A'])
        audit.update(measurement=measurement,dwell_s=dwell,timing_status='20 ms nominal from 50 Hz / assumed PLC=1; actual timestamps unavailable')
        states=[];records=[];holder={};started=time.perf_counter()
        class WallStop(Exception):pass
        def observe(i,model,status):
            states.append(model.state.reshape(model.shape).copy())
            records.append(dict(segment=i,scheduled_end_s=(i+1)*dwell,**model.iv[-1]))
            holder['model']=model
            if i%10==0 or model.connected:print('point',i,'V',model.command_voltage,'events',model.steps,'connected',model.connected,flush=True)
            if time.perf_counter()-started>=args.max_wall_s and not model.connected and i<len(protocol)-1:raise WallStop
        try:
            m,summary,_,_=run_forming(cfg,e,protocol,on_segment=observe,model_class=CachedFormingModel)
        except WallStop:
            m=holder['model'];summary=dict(m.result('wall_budget_unknown'),first_bridge=m.first_bridge,
                         electronics=asdict(e),post_bridge_mode='frozen_morphology_readout')
        summary.update(fresh_forming_run=True,connected=m.connected,completed_points=sum(r['time_s']>=r['scheduled_end_s']-1e-12 for r in records),
                       total_points=len(protocol),source_mode='fixed_infinite_ag_reservoir',
                       morphology=metrics(m.state.reshape(m.shape),m.c.spacing_nm,m.c.layer_cells))
        np.savez_compressed(folder/'states.npz',state=m.state.reshape(m.shape),site_energy_eV=m.site.reshape(m.shape),
                            vacancy_mask=m.vacancy_mask.reshape(m.shape),scan_states=np.asarray(states))
        if records:
            with (folder/'forming_trace.csv').open('w',encoding='utf-8-sig',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
            fig,ax=plt.subplots(figsize=(7,4.5),constrained_layout=True)
            ax.semilogy([r['command_voltage_V'] for r in records],np.maximum([abs(r['current_A']) for r in records],1e-15),label='Restored model')
            ax.semilogy([p['voltage_V'] for p in measurement['points']],np.maximum([abs(p['measured_A']) for p in measurement['points']],1e-15),label='Measured')
            ax.set(xlabel='Command voltage (V)',ylabel='Absolute current (A)',title='Fresh run: executed model points only');ax.legend();ax.grid(alpha=.2)
            fig.savefig(folder/'forming_curve.png',dpi=180);plt.close(fig)
        draw_morphology(folder,m)
    write_json(folder/'manifest.json',audit);write_json(folder/'result.json',summary)
    print(json.dumps({k:summary[k] for k in ('status','fresh_forming_run','connected','ag_balance_error')},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
