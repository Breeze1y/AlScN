"""Run independent SET/RESET cycles with field-selected, supply-limited Ag growth."""
from pathlib import Path
from dataclasses import asdict, replace
import argparse
import hashlib
import json
import shutil
import time

import numpy as np

import run_ten_cycles as old
import run_full_chain as full
import cycle_unified40 as cycle
import run_local_retreat as local
import constrained_retreat40 as retreat_kernel
from cycle_activity40 import ActivityState
from form_hourglass40 import span
from field_growth40 import FieldGrowthParams, FieldGrowthSequence, checkpoint_at_event

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/ten_cycles_local_retreat40'
DEFAULT=ROOT/'results/ten_cycles_ag_constrained40'


def run_set(root,number,checkpoint,positive,network,grid,steps,growth,*,checkpoint_writer=checkpoint_at_event):
    target=root/f'cycle_{number:02d}'/'SET'
    folder=target.parent/'catalogues/SET'
    if (target/'summary.json').exists() and (target/'checkpoint/restart.npz').exists():
        return target/'checkpoint'
    start=time.perf_counter()
    if (folder/'sequence.npz').exists():
        sequence=FieldGrowthSequence.restore(folder)
        assert asdict(sequence.params)==asdict(growth), 'Growth parameters changed on resume'
    else:
        sequence=FieldGrowthSequence.from_checkpoint(checkpoint,params=growth)
        sequence.extend(96)
        sequence.save(folder)
    old.check_initial(sequence,checkpoint)
    initial=ActivityState.from_dict(full.js(checkpoint/'restart.json')['electronic_state'])
    v,reference=full.measured('SET')
    while True:
        phase=old.ContinuousPhase(folder,network,grid)
        records,end,info=cycle.simulate(phase,v,initial,positive['rate_scale'],
            positive['voltage_scale_V'],cycle.CircuitParams(**positive['circuit']),substeps=steps)
        cache_info=phase.lookup.save()
        if not info['sequence_exhausted']:
            break
        if len(sequence.events)>12000:
            raise RuntimeError('SET reached the declared catalogue budget without completing the sweep')
        print(f'cycle {number} SET extends field catalogue from {len(sequence.events)}',flush=True)
        sequence.extend(128)
        sequence.save(folder)
    s=cycle.save_phase(target,phase,records,initial,end,info,positive,reference)
    a=full.load(target/'states.npz')
    a['metal_connected']=np.array([bool(span(m>0).any()) for m in a['metal_id']])
    a['branch']=np.where(np.arange(len(v))<=100,'outward','return')
    a['global_elapsed_units']=initial.elapsed_units+np.arange(1,len(v)+1)
    full.save_arrays(target,a)
    cp=target/'checkpoint'
    checkpoint_writer(folder,info['final_event'],cp,info['event_exposure'],end.to_dict())
    full.attach_electronics(cp,end.to_dict())
    s.update(old.source_info(checkpoint,folder,number,'SET'),network_params=asdict(network),
        growth=asdict(growth),integration_substeps=steps,
        parameter_role='Same frozen growth, circuit and RESET parameters across all newly simulated cycles',
        final_Ag=int(a['Ag_total'][-1]),final_read_ohm=float(1/a['read_G_S'][-1]),
        final_metal_connected=bool(a['metal_connected'][-1]),lazy_network=cache_info,
        elapsed_seconds=time.perf_counter()-start)
    full.dump(target/'summary.json',s)
    print(f'cycle {number} SET complete Ag={s["final_Ag"]} R={s["final_read_ohm"]:.6g} connected={s["final_metal_connected"]}',flush=True)
    return cp


def baseline_hashes():
    return {str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(BASE.rglob('*')) if p.is_file()}


def run_reset(root,number,checkpoint,negative,network,grid,retreat,steps,*,wait_for_return_route=False,
              share_residual_retreat=False,checkpoint_writer=old.erosion_checkpoint):
    """Use the existing retreat law with actual live top-source route checks."""
    target=root/f'cycle_{number:02d}'/'RESET'
    folder=target.parent/'catalogues/RESET'
    if (target/'summary.json').exists() and (target/'checkpoint/restart.npz').exists():
        saved_mode=full.js(target/'summary.json').get('retreat_mode')
        if (saved_mode=='shared_residual')!=share_residual_retreat:
            raise ValueError('Existing RESET output belongs to another retreat mode')
        return target/'checkpoint'
    start=time.perf_counter()
    if 'background_oxidation_weight' in negative and not (folder/'sequence.npz').exists():
        sequence = full.new_sequence('erosion', checkpoint, network)
        sequence.params = replace(sequence.params, background_weight=negative['background_oxidation_weight'])
        sequence.extend(320).save(folder)
    sequence=old.prepare_sequence(folder,'erosion',checkpoint,network)
    while True:
        cat=old.ResetCatalogue(folder,checkpoint,negative,network,grid)
        cat.share_residual_retreat=share_residual_retreat
        try:
            a,s,actions,snaps=retreat_kernel.simulate(cat,retreat,steps,wait_for_return_route=wait_for_return_route)
        except RuntimeError as exc:
            cat.lookup.save()
            if str(exc)!='Spatial catalogue exhausted':
                raise
            sequence.extend(160).save(folder)
            continue
        cache_info=cat.lookup.save()
        break
    assert not s['pending_microsteps']
    n=retreat.subunits
    assert np.all((a['Ag_units'][-1]==0)|(a['Ag_units'][-1]==n))
    a['global_elapsed_units']=s['initial_electronic_state']['elapsed_units']+a['elapsed_units']
    full.save_arrays(target,a)
    full.dump(target/'snapshots.json',snaps[1])
    np.savez_compressed(target/'snapshots.npz',Ag_units=snaps[0],VN=cat.vn,subunits=n)
    np.savez_compressed(target/'actual_transport.npz',
        offsets=np.r_[0,np.cumsum([len(r) for r in cat.actual_routes])],
        nodes=np.concatenate(cat.actual_routes) if cat.actual_routes else np.array([],np.int16))
    if actions:
        import csv
        with (target/'events.csv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=actions[0]);w.writeheader();w.writerows(actions)
    cp=target/'checkpoint'
    checkpoint_writer(folder,s['final_event'],cp,s['final_event_exposure'])
    full.attach_electronics(cp,s['final_electronic_state'])
    saved=full.load(cp/'restart.npz')
    np.testing.assert_array_equal(saved['metal_id'],a['metal_id'][-1])
    saved.update(Ag_units=a['Ag_units'][-1],subunits=n)
    np.savez_compressed(cp/'restart.npz',**saved)
    full.dump(cp/'local_retreat.json',dict(parameters=asdict(retreat),trigger=s['trigger'],
        pending_microsteps=[],elapsed_units=s['elapsed_units'],
        observation_exposure=s['final_event_exposure'],final_event=s['final_event'],
        **({'retreat_mode':'shared_residual'} if share_residual_retreat else {})))
    s.update(old.source_info(checkpoint,folder,number,'RESET'),source_chain=str(root.relative_to(ROOT)),
        network_params=asdict(network),final_Ag=float(a['Ag_total'][-1]),
        final_read_ohm=float(1/a['read_G_S'][-1]),final_metal_connected=bool(a['metal_connected'][-1]),
        lazy_network=cache_info,elapsed_seconds=time.perf_counter()-start,
        return_kernel='constrained_retreat40: live top-source endpoint and empty-route checks')
    if share_residual_retreat:
        s['retreat_mode']='shared_residual'
        s['residual_cohort']=getattr(cat,'residual_cohort',None)
    full.dump(target/'summary.json',s)
    print(f'cycle {number} RESET complete Ag={s["final_Ag"]} R={s["final_read_ohm"]:.6g}',flush=True)
    return cp


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cycles',type=int,default=10)
    ap.add_argument('--out',type=Path,default=DEFAULT)
    ap.add_argument('--start-checkpoint',type=Path)
    ap.add_argument('--field-exponent',type=float,default=2.)
    ap.add_argument('--crowding-strength',type=float,default=2.)
    ap.add_argument('--soft-volume',type=float,default=1000.)
    ap.add_argument('--soft-scale',type=float,default=500.)
    ap.add_argument('--no-supply-limit',action='store_true')
    ap.add_argument('--reset-rate-factor',type=float,default=1.)
    ap.add_argument('--substeps',type=int,default=16)
    args=ap.parse_args()
    assert 1<=args.cycles<=10
    root=args.out.resolve()
    if root==BASE.resolve():
        raise ValueError('The baseline result directory is read-only for this experiment')
    root.mkdir(parents=True,exist_ok=True)
    growth=FieldGrowthParams(field_exponent=args.field_exponent,
        crowding_strength=args.crowding_strength,soft_volume=args.soft_volume,soft_scale=args.soft_scale,
        supply_enabled=not args.no_supply_limit)
    positive,negative,network=full.common_parameters()
    assert args.reset_rate_factor>0
    negative=dict(negative,rate_scale=negative['rate_scale']*args.reset_rate_factor)
    grid=np.array(full.js(old.BASE/'parameters.json')['logq_grid'])
    retreat=local.Params(**full.js(old.FIRST_RESET/'summary.json')['parameters'])
    configuration=dict(schema='field_constrained_cycles40/v1',positive=positive,negative=negative,
        network=asdict(network),local_retreat=asdict(retreat),growth=asdict(growth),logq_grid=grid.tolist(),
        substeps=args.substeps,requested_cycles=args.cycles,
        start_checkpoint=str(args.start_checkpoint.resolve()) if args.start_checkpoint else None,
        first_cycle='Every SET/RESET cycle is newly simulated; FORMING and INIT_RESET inherited unchanged',
        units='0.5 nm coarse Ag volumes; effective dwell units, not calibrated seconds')
    params=root/'parameters.json'
    if params.exists():
        assert full.js(params)==configuration, 'Resume parameters do not match'
    else:
        full.dump(params,configuration)
    protected=root/'baseline_sha256.json'
    if not protected.exists():
        full.dump(protected,baseline_hashes())
    if args.start_checkpoint:
        checkpoint=args.start_checkpoint.resolve()
    else:
        for name in ['FORMING','INIT_RESET']:
            if not (root/name).exists():
                shutil.copytree(BASE/name,root/name)
        checkpoint=root/'INIT_RESET/checkpoint'
    full.dump(root/'inherited_sources.json',dict(initial_checkpoint=str(checkpoint),
              initial_checkpoint_sha256=full.sha(checkpoint/'restart.npz'),baseline=str(BASE)))
    try:
        for number in range(1,args.cycles+1):
            full.dump(root/'progress.json',dict(status='running',cycle=number,phase='SET',requested_cycles=args.cycles))
            checkpoint=run_set(root,number,checkpoint,positive,network,grid,args.substeps,growth)
            full.dump(root/'progress.json',dict(status='running',cycle=number,phase='RESET',requested_cycles=args.cycles))
            checkpoint=run_reset(root,number,checkpoint,negative,network,grid,retreat,args.substeps)
            old.update_summary(root,args.cycles)
        rows=old.update_summary(root,args.cycles,complete=True)
        summary=full.js(root/'chain_summary.json')
        summary.update(schema='field_constrained_cycles40/chain/v1',
            sequence=f'Inherited initial state -> (new SET -> new RESET) x {args.cycles}',
            growth=asdict(growth),parameter_fit=False,
            initial_checkpoint=full.js(root/'inherited_sources.json')['initial_checkpoint'])
        full.dump(root/'chain_summary.json',summary)
        assert baseline_hashes()==full.js(protected), 'Baseline artifact changed'
        full.dump(root/'progress.json',dict(status='complete',completed_cycles=len(rows),baseline_unchanged=True))
        print(json.dumps(dict(complete=True,cycles=len(rows),out=str(root))),flush=True)
    except Exception as exc:
        full.dump(root/'failure.json',dict(error_type=type(exc).__name__,message=str(exc),
                                         progress=full.js(root/'progress.json')))
        raise


if __name__=='__main__':
    main()
