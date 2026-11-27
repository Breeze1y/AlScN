"""Local Ag retreat with live source validation after each volume removal.

Integration and contact physics match multi_cycle_retreat.simulate. Before
using a return path, its host nodes and its connection to a live top reservoir
are checked on a geometry with the oxidizing site opened. Invalid inherited
or cached paths are recomputed on the current geometry. Event timing is
unchanged by default. The long-cycle run can delay an event or a microstep
until its return route reopens. Effective transport routes are saved for
every volume action.
"""
from pathlib import Path
from dataclasses import asdict
import heapq
import json
import math

import numpy as np
from scipy.ndimage import binary_propagation
from cycle_activity40 import ActivityState
from cycle_unified40 import CircuitParams, read_circuit
from form_hourglass40 import span
from rate_reset40 import RateParams, return_path
from run_local_retreat import Params, metrics
import local_fraction_network as net


class BlockedReturnRoute(RuntimeError):
    """No currently empty route exists; no Ag is removed or bypassed."""
    def __init__(self, event_id, y, x, time, pending):
        self.event_id=event_id;self.y=y;self.x=x;self.time=time;self.pending=pending
        super().__init__(f'No empty return route: event={event_id}, y={y}, x={x}, time={time}, pending={pending}')


def simulate(cat,p,substeps=16,*,wait_for_return_route=False,retreat_offsets=None,
             scale_retreat_by_faces=False,network_read=None):
    """Optional offsets schedule microsteps after the initial unit removal.

    When requested, each initially exposed in-plane face contributes equal
    flux, shortening that site's retreat by its exposed face count.
    """
    if retreat_offsets is not None:
        retreat_offsets=np.asarray(retreat_offsets,float)
        if (retreat_offsets.shape!=(p.subunits-1,) or not np.isfinite(retreat_offsets).all()
                or np.any(np.diff(np.r_[0.,retreat_offsets])<=0)
                or not np.isclose(retreat_offsets[-1],p.retreat_duration)):
            raise ValueError('Retreat offsets must increase from >0 to retreat_duration.')
    cfg=cat.summary['clock'];circ=CircuitParams(**cfg['circuit']);state=ActivityState.from_dict(cat.summary['initial_electronic_state'])
    units=((cat.ids[0]>0)*p.subunits).astype(np.uint8);initial_units=int(units.sum())
    index=0;exposure=0.;now=0.;queue=[];serial=0;trigger=None;records=[];states=[];actions=[];snapshots=[];snapshot_meta=[];exact_solves=0
    rate_params=RateParams(**json.loads((Path(cat.folder)/'restart.json').read_text(encoding='utf-8-sig'))['params'])
    cat.actual_routes=[]
    rerouted={}
    blocked_since=None
    paused_steps={}
    def route_valid(route,y,x):
        if len(route)==0 or int(route[0])!=y*units.shape[1]+x:
            return False
        opened=units>0
        opened[y,x]=False
        if np.any(opened.ravel()[route[1:]]):
            return False
        ey,ex=divmod(int(route[-1]),units.shape[1])
        if ey==0:
            return True
        seed=np.zeros_like(opened)
        seed[0]=opened[0]
        top=binary_propagation(seed,mask=opened)
        return any(top[yy,xx] for yy,xx in [(ey-1,ex),(ey+1,ex),(ey,ex-1),(ey,ex+1)]
                   if 0<=yy<units.shape[0] and 0<=xx<units.shape[1])
    def resolve_route(y,x,event_id):
        route=np.asarray(cat.routes[event_id-1],dtype=np.int16)
        if route_valid(route,y,x):
            return route,'inherited'
        cached=rerouted.get(event_id)
        if cached is not None and route_valid(cached,y,x):
            return cached,'rerouted'
        alternative,travel=return_path(units>0,y*units.shape[1]+x,cat.vn,rate_params)
        if alternative is None:
            raise BlockedReturnRoute(event_id,y,x,now,len(queue))
        route=np.asarray(alternative,dtype=np.int16)
        assert route_valid(route,y,x),'Recomputed route lacks a live top reservoir'
        rerouted[event_id]=route
        return route,'rerouted'
    connected=bool(span(units>0).any());actual_break_recorded=False
    def read(v):
        nonlocal exact_solves
        if network_read is not None:
            g=network_read(units,state.q,len(actions))
        elif not np.array_equal(units,(cat.ids[index]>0)*p.subunits):
            g=net.solve(cat.topology,units.astype(float)/p.subunits,cat.network,activity=state.q,contact_exponent=p.contact_exponent)['conductance_S'];exact_solves+=1
        else:g=cat.binary_g(index,state.q)
        return read_circuit(g,state.q,float(v),circ)
    def snap(label,point,v):
        snapshots.append(units.copy());snapshot_meta.append(dict(label=label,point=point,time=now,voltage_V=float(v),q=state.q,event_index=index,current_A=read(v)['current_A']))
    def remove(y,x,amount,point,v,kind,event_id,step,*,resolved_route=None,route_delay=0.):
        nonlocal connected,actual_break_recorded
        before=connected;previous=int(units[y,x]);assert previous>=amount
        route,route_kind=resolve_route(y,x,event_id) if resolved_route is None else resolved_route
        assert int(route[0])==y*units.shape[1]+x
        assert not np.any(units.ravel()[route[1:]])
        route_index=len(cat.actual_routes)
        cat.actual_routes.append(route.copy())
        trial=units.copy();trial[y,x]-=amount;after=bool(span(trial>0).any())
        if before and not after and not actual_break_recorded:snap('metal_break_before',point,v)
        units[y,x]-=amount;connected=after
        actions.append(dict(time=now,point=point,voltage_V=float(v),event_index=event_id,y=y,x=x,kind=kind,
            unit_step=step,removed_units=amount,remaining_units=int(units[y,x]),Ag_total=float(units.sum()/p.subunits),
            returned_effective=float(cat.returned+(initial_units-int(units.sum()))/p.subunits),q=state.q,
            connected_before=before,connected_after=after,route_index=route_index,
            route_kind=route_kind,route_nodes=len(route),route_delay_units=route_delay))
        if before and not after and not actual_break_recorded:
            actual_break_recorded=True;snap('metal_break_after',point,v)
            if getattr(cat,'share_residual_retreat',False) and queue:
                # The neck surfaces already dissolving at the first disconnection
                # share a single site's remaining flux. Intrinsic progress u
                # takes dt/du=N_active; each site's last due step frees its share.
                remaining_by_site={}
                for due,_,_,_,site,_ in queue:
                    remaining_by_site[site]=max(0.,due-now,remaining_by_site.get(site,0.))
                remaining_times=tuple(remaining_by_site.values())
                queue[:]=[(now+sum(min(max(0.,due-now),tail) for tail in remaining_times),
                           order,yy,xx,site,step) for due,order,yy,xx,site,step in queue]
                heapq.heapify(queue)
                cat.residual_cohort=dict(sites=len(remaining_times),remaining_times=remaining_times)
    def resume_microsteps(point,v):
        # A site can reopen another site's route only when it becomes empty.
        cleared=True
        while cleared:
            cleared=False
            for ev,tasks in list(paused_steps.items()):
                due,_,y,x,_,step=tasks[0]
                delay=float(now-due)
                try:
                    remove(y,x,1,point,v,'retreat_microstep',ev,step,route_delay=delay)
                except BlockedReturnRoute:
                    continue
                del paused_steps[ev]
                for scheduled,order,yy,xx,event_id,next_step in tasks[1:]:
                    heapq.heappush(queue,(scheduled+delay,order,yy,xx,event_id,next_step))
                cleared |= units[y,x]==0
    for point,v in enumerate(cat.base['voltage_V']):
        remaining=1.;mean_i=mean_terminal=mean_bulk=0.;before=index;micro_before=len(actions);count=0
        while remaining>1e-10:
            rd=read(v);speed=cfg['rate_scale']*math.expm1(min(abs(rd['material_voltage_V'])/cfg['voltage_scale_V'],60.))
            if blocked_since is not None:
                speed=0.  # Hold this catalogue event at threshold while earlier retreats finish.
            if index+1>=len(cat.clock) and speed>0:raise RuntimeError('Spatial catalogue exhausted')
            to_event=(cat.clock[index+1]-exposure)/speed if speed>0 and index+1<len(cat.clock) else math.inf
            to_micro=max(0.,queue[0][0]-now) if queue else math.inf
            cap=min(remaining,1/substeps);dt=max(0.,min(cap,to_event,to_micro))
            event_hit=to_event<=dt+1e-12;micro_hit=to_micro<=dt+1e-12
            mean_i+=rd['current_A']*dt;mean_terminal+=rd['terminal_voltage_V']*dt;mean_bulk+=rd['material_voltage_V']*dt
            state.advance(rd['material_voltage_V'],dt,current_A=max(rd['material_current_A'],0.))
            exposure+=speed*dt;remaining-=dt;now+=dt;count+=1
            cleared_site=False
            if micro_hit:
                while queue and queue[0][0]<=now+1e-10:
                    task=heapq.heappop(queue)
                    _,_,y,x,ev,step=task
                    # Keep binary/partial read mode valid while recording a last microstep.
                    try:
                        remove(y,x,1,point,v,'retreat_microstep',ev,step)
                    except BlockedReturnRoute:
                        if not wait_for_return_route:
                            raise
                        # Shift the remaining schedule together when this step can resume.
                        paused_steps[ev]=[task]+sorted(item for item in queue if item[4]==ev)
                        queue[:]=[item for item in queue if item[4]!=ev]
                        heapq.heapify(queue)
                        continue
                    cleared_site |= units[y,x]==0
                if cleared_site and paused_steps:
                    resume_microsteps(point,v)
            if event_hit or (blocked_since is not None and cleared_site):
                event_id=index+1
                exposure=float(cat.clock[event_id]);e=cat.events[index];y,x=int(e['y']),int(e['x'])
                try:
                    resolved_route=resolve_route(y,x,event_id)
                except BlockedReturnRoute:
                    if not wait_for_return_route or not queue:
                        raise
                    if blocked_since is None:
                        blocked_since=now
                    continue
                route_delay=0. if blocked_since is None else float(now-blocked_since)
                blocked_since=None
                index=event_id
                assert units[y,x]==p.subunits,'A catalogue site may be initiated only once'
                if p.enabled and trigger is None and connected:
                    trial=units.copy();trial[y,x]=0
                    if not span(trial>0).any():
                        trigger=dict(event_index=index,point=point,time=now,voltage_V=float(v),y=y,x=x,until=now+p.local_window)
                        snap('retreat_onset_before',point,v)
                local=p.enabled and trigger is not None and now<=trigger['until']+1e-12 and math.hypot(y-trigger['y'],x-trigger['x'])*.5<=p.radius_nm+1e-12
                postbreak_duration=0.
                if not local and actual_break_recorded and hasattr(cat,'postbreak_retreat_duration'):
                    postbreak_duration=cat.postbreak_retreat_duration(trigger,y,x,units)
                    local=postbreak_duration>0
                if local:
                    faces=1
                    if scale_retreat_by_faces:
                        neighbors=((y-1,x),(y+1,x),(y,x-1),(y,x+1))
                        faces=max(1,sum(1 for yy,xx in neighbors
                                        if 0<=yy<units.shape[0] and 0<=xx<units.shape[1]
                                        and units[yy,xx]==0))
                    for step in range(2,p.subunits+1):
                        offset=(p.retreat_duration*(step-1)/(p.subunits-1) if retreat_offsets is None
                                else retreat_offsets[step-2])
                        if postbreak_duration:
                            offset=postbreak_duration*(step-1)/(p.subunits-1)
                        serial+=1;heapq.heappush(queue,(now+offset/faces,serial,y,x,index,step))
                    remove(y,x,1,point,v,'retreat_start',index,1,resolved_route=resolved_route,route_delay=route_delay)
                    if trigger['event_index']==index:snap('retreat_onset_after',point,v)
                else:
                    remove(y,x,p.subunits,point,v,'whole_cell',index,1,resolved_route=resolved_route,route_delay=route_delay)
                    if paused_steps:
                        resume_microsteps(point,v)
            if count>20000:raise RuntimeError('Integration did not advance')
        last=read(v);total=float(units.sum()/p.subunits);returned=cat.returned+(initial_units-int(units.sum()))/p.subunits
        assert abs(total-cat.injected+returned)<1e-12
        records.append(dict(point=point,voltage_V=float(v),current_A=mean_i,endpoint_current_A=last['current_A'],
            terminal_voltage_V=mean_terminal,material_voltage_V=mean_bulk,read_G_S=last['read_G_S'],activity=state.q,D=state.D,H=state.H,
            event_index=index,events=index-before,microevents=len(actions)-micro_before,event_exposure=exposure,
            Ag_total=total,Ag_occupied_cells=int(np.count_nonzero(units)),injected=cat.injected,returned=returned,
            metal_connected=connected,pending_retreats=len({item[4] for item in queue})+len(paused_steps),elapsed_units=now))
        states.append(units.copy())
    arrays={k:np.array([r[k] for r in records]) for k in records[0]}
    arrays.update(Ag_units=np.asarray(states),VN=cat.vn,subunits=np.array(p.subunits),initial_metal_id=cat.ids[0],
        metal_id=np.where(np.asarray(states)>0,cat.ids[0],0),measured_A=cat.base['measured_A'],branch=cat.base['branch'])
    summary=dict(schema='constrained_retreat40/run/v1',parameters=asdict(p),clock=cfg,network=asdict(cat.network),integration_substeps=substeps,
        source_chain='results/full_chain_refined',source_checkpoint='results/full_chain_refined/SET/checkpoint',
        initial_electronic_state=cat.summary['initial_electronic_state'],final_electronic_state=state.to_dict(),
        trigger=trigger,final_event=index,final_event_exposure=exposure,elapsed_units=now,
        pending_microsteps=[list(item) for item in queue]+[list(item) for tasks in paused_steps.values() for item in tasks],
        microsteps=len(actions),fractional_sites=len({r['event_index'] for r in actions if r['kind']=='retreat_start'}),exact_fractional_solves=exact_solves,
        metrics=metrics(arrays,cat.base['measured_A']),
        semantics='Catalogue entries initiate dissolution; Ag effective volume returns only on actual microsteps. Fixed spatial q=1 event order; dynamic q and event timing recomputed. Every volume removal checks its actual dielectric return route; a blocked inherited path is recomputed without traversing remaining Ag.',
        time_units='effective sampling units; no conversion to seconds',Ag_units='0.5 nm coarse-volume packets, subdivided into effective volumes, not atoms',
        route_policy='Before every volume action, the inherited or cached dielectric route must have empty host nodes and end at a live top reservoir with the oxidizing site opened. An invalid path is recomputed on the current geometry at the original reference bias. Every actual route is saved.',
        route_kernel='constrained_retreat40/v1',
        rerouted_catalogue_events=len(rerouted),rerouted_volume_actions=sum(a['route_kind']=='rerouted' for a in actions),
        delayed_return_actions=sum(1 for a in actions if a['route_delay_units']>0))
    if wait_for_return_route:
        summary['blocked_event_policy']='A blocked catalogue event waits at its threshold. A blocked microstep delays the remaining steps at its site. Every removal still requires a live return route.'
    return arrays,summary,actions,(np.asarray(snapshots),snapshot_meta)
