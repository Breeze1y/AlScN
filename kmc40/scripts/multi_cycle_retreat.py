"""Local Ag retreat with actual dielectric return routes for continued cycles.

The integration and contact physics match run_local_retreat.simulate. Only a
return route blocked by remaining fractional Ag is recomputed, using the same
reference-voltage route kernel and kinetic parameters. The actual route for
each removal is available as cat.actual_routes[action['route_index']].
"""
from pathlib import Path
from dataclasses import asdict
import heapq
import json
import math

import numpy as np
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


def simulate(cat,p,substeps=16):
    cfg=cat.summary['clock'];circ=CircuitParams(**cfg['circuit']);state=ActivityState.from_dict(cat.summary['initial_electronic_state'])
    units=((cat.ids[0]>0)*p.subunits).astype(np.uint8);initial_units=int(units.sum())
    index=0;exposure=0.;now=0.;queue=[];serial=0;trigger=None;records=[];states=[];actions=[];snapshots=[];snapshot_meta=[];exact_solves=0
    rate_params=RateParams(**json.loads((Path(cat.folder)/'restart.json').read_text(encoding='utf-8-sig'))['params'])
    cat.actual_routes=[]
    rerouted={}
    def resolve_route(y,x,event_id):
        route=np.asarray(cat.routes[event_id-1],dtype=np.int16)
        if np.any(units.ravel()[route[1:]]):
            if event_id not in rerouted:
                alternative,travel=return_path(units>0,y*units.shape[1]+x,cat.vn,rate_params)
                if alternative is None:
                    raise BlockedReturnRoute(event_id,y,x,now,len(queue))
                rerouted[event_id]=np.asarray(alternative,dtype=np.int16)
            route=rerouted[event_id]
            assert not np.any(units.ravel()[route[1:]])
            return route,'rerouted'
        return route,'inherited'
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
        route,route_kind=resolve_route(y,x,event_id)
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
            route_kind=route_kind,route_nodes=len(route),route_delay_units=0.))
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
        semantics='Catalogue entries initiate dissolution; Ag effective volume returns only on actual microsteps. Fixed spatial q=1 event order; dynamic q and event timing recomputed. Every volume removal checks its actual dielectric return route; a blocked inherited path is recomputed without traversing remaining Ag.',
        time_units='effective sampling units; no conversion to seconds',Ag_units='0.5 nm coarse-volume packets, subdivided into effective volumes, not atoms',
        route_policy='Inherited dielectric route is used when empty; otherwise a shortest-time route is recomputed on the actual remaining-Ag mask at the original reference bias. Every volume removal stores its actual route and checks that its non-source nodes are empty.',
        rerouted_catalogue_events=len(rerouted),rerouted_volume_actions=sum(a['route_kind']=='rerouted' for a in actions),
        delayed_return_actions=0)
    return arrays,summary,actions,(np.asarray(snapshots),snapshot_meta)
