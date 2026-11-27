"""Geometry-field SET growth with local crowding and soft Ag supply feedback.

The spatial catalogue remains at a fixed reference bias. A full host-grid
quasistatic conduction-field surrogate (metal/host conductivity contrast)
replaces the old empirical depth and downward-growth multipliers. This field
is used for site selection; terminal I-V still comes from the existing mixed
electronic network. Ion delivery retains the established empty-host routes.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import argparse
import copy
import csv
import json
import math
import tempfile

import numpy as np
from scipy.ndimage import convolve
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.sparse.linalg import spsolve

import form_hourglass40 as geom
import source_cycle_vn40 as source
from set_growth40 import GrowthParams, GrowthSequence, validate_routes, _hash


@dataclass(frozen=True)
class FieldGrowthParams(GrowthParams):
    use_field: bool = True
    field_exponent: float = 2.
    field_floor: float = .05
    metal_host_contrast: float = 1.e4
    crowding_strength: float = 2.
    crowding_radius_nm: float = 2.
    supply_enabled: bool = True
    soft_volume: float = 1000.
    soft_scale: float = 500.


def geometry_field(metal, contrast=1.e4):
    """Unit terminal-bias Laplace field; electrode faces at y=0 and y=40 nm.

    Cell-centred finite volumes, harmonic edge conductivities, half-cell
    electrode links, and insulating lateral boundaries. Host sigma is one;
    finite high metallic sigma keeps even bridged geometries well posed.
    Returned enhancement is the maximum host-side normal field at an exposed
    metal/electrode contact divided by the uniform empty-host field.
    """
    ny, nx = metal.shape
    nodes = np.arange(ny*nx).reshape(ny, nx)
    a = np.r_[nodes[:-1].ravel(), nodes[:, :-1].ravel()]
    b = np.r_[nodes[1:].ravel(), nodes[:, 1:].ravel()]
    sigma = np.where(metal.ravel(), contrast, 1.)
    edge = 2*sigma[a]*sigma[b]/(sigma[a]+sigma[b])
    top, bottom = nodes[0], nodes[-1]
    gt, gb = 2*sigma[top], 2*sigma[bottom]
    n = ny*nx
    diagonal = np.bincount(np.r_[a,b,top,bottom], np.r_[edge,edge,gt,gb], minlength=n)
    rhs = np.bincount(top, gt, minlength=n)
    matrix = coo_matrix((np.r_[diagonal,-edge,-edge],
                        (np.r_[np.arange(n),a,b], np.r_[np.arange(n),b,a])),
                       shape=(n,n)).tocsc()
    phi = spsolve(matrix, rhs, permc_spec='MMD_AT_PLUS_A')
    if not np.isfinite(phi).all():
        raise FloatingPointError('Nonfinite growth geometry field')
    if phi.min() < -1.e-8 or phi.max() > 1+1.e-8:
        raise FloatingPointError('Growth potential outside boundary range')
    normal = np.zeros(n)
    crossing = metal.ravel()[a] != metal.ravel()[b]
    host = np.where(metal.ravel()[a[crossing]],b[crossing],a[crossing])
    np.maximum.at(normal,host,np.abs(edge[crossing]*(phi[a[crossing]]-phi[b[crossing]])))
    normal[top] = np.maximum(normal[top],np.abs(gt*(1-phi[top])))
    normal[bottom] = np.maximum(normal[bottom],np.abs(gb*phi[bottom]))
    return phi.reshape(ny,nx), normal.reshape(ny,nx)*ny


def supply_factor(mass, params):
    return math.exp(-max(float(mass)-params.soft_volume,0.)/params.soft_scale) if params.supply_enabled else 1.


class FieldPilot(source.Pilot):
    def __init__(self, seed, params, vn_fraction):
        super().__init__(seed,params.eta,params.assist_eV,params.pt_mode,vn_fraction)
        self.field_params = params
        radius = int(math.ceil(params.crowding_radius_nm/.5))
        yy,xx=np.mgrid[-radius:radius+1,-radius:radius+1]
        self.crowding_kernel=((yy*.5)**2+(xx*.5)**2 <= params.crowding_radius_nm**2+1.e-12).astype(float)

    def deposit(self,v,phase,cycle,point):
        metal=self.metallic();src,reservoir=self.choose_source()
        if src is None:return False
        graph,cn,pair=source.source_routes(metal,v,self.vn,self.eta,self.assist)
        dist,pred=dijkstra(graph,indices=src,return_predecessors=True)
        front=(~metal)&(cn>0);front[0]=~metal[0]
        if self.pt_mode=='arrival_pt':front[-1]|=~metal[-1]
        y,x=np.where(front & np.isfinite(dist.reshape(source.NY,source.NX)))
        if not len(y):return False
        p=self.field_params
        if p.use_field:
            phi,enhancement=geometry_field(metal,p.metal_host_contrast)
        else:
            phi=np.broadcast_to(1-(np.arange(source.NY)[:,None]+.5)/source.NY,metal.shape)
            enhancement=np.ones(metal.shape)
        rho=convolve(metal.astype(float),self.crowding_kernel,mode='constant')
        rho/=convolve(np.ones(metal.shape),self.crowding_kernel,mode='constant')
        focus=np.maximum(enhancement[y,x],p.field_floor)**p.field_exponent
        crowd=np.exp(-p.crowding_strength*rho[y,x])
        adhesion=.08*((y==0)|(y==source.NY-1))
        dg=.075-.035*cn[y,x]+self.eta*(-.5-2.48)*pair[y,x]-.08*v-adhesion
        barrier=.30+np.maximum(dg,0.)
        capture=np.exp(-.5*((x-src%source.NX)*.5/3.)**2)
        redrate=1e3*np.exp(-barrier/source.KT)
        travel=dist[y*source.NX+x]
        w=focus*crowd*capture/(1/redrate+travel)
        if not np.isfinite(w).all() or not w.sum()>0:
            raise FloatingPointError('Invalid field-growth candidate weights')
        j=int(self.rng.choice(len(y),p=w/w.sum()))
        target=int(y[j]*source.NX+x[j]);route=[target]
        while route[-1]!=src:
            parent=int(pred[route[-1]])
            if parent<0:raise RuntimeError('Unreachable reaction site')
            route.append(parent)
        route=route[::-1]
        self.injected+=1;agid=self.injected
        assert not np.any(metal.ravel()[route])
        acceptance=-np.expm1(-1000./(travel[j]+1/redrate[j]))
        accepted=bool(self.rng.random()<acceptance)
        if accepted:self.ids[y[j],x[j]]=agid
        else:
            self.returned+=1
            route=route+route[-2::-1]
        self.paths.append(np.array(route,np.int16));self.reservoir_paths.append(reservoir)
        self.events.append(dict(event=len(self.events),phase=phase,cycle=cycle,point=point,voltage_V=v,
            kind='inject_transport_reduce' if accepted else 'inject_unreacted_return',Ag_id=agid,
            source_x=src%source.NX,source_y=src//source.NX,y=int(y[j]),x=int(x[j]),
            path_nodes=len(route),mean_path_us=float(travel[j]),redox_delta_eV=float(dg[j]),
            forward_barrier_eV=float(barrier[j]),reverse_barrier_eV=float(.30+max(-dg[j],0)),
            source_kind='top_connected_Ag_reservoir',delivery_probability=float(acceptance),
            injected=self.injected,returned=self.returned,
            field_enhancement=float(enhancement[y[j],x[j]]),field_weight=float(focus[j]),
            local_Ag_fraction=float(rho[y[j],x[j]]),crowding_weight=float(crowd[j]),
            local_unit_potential=float(phi[y[j],x[j]])))
        return accepted


class FieldGrowthSequence(GrowthSequence):
    def __init__(self,ids,vn,provenance,params=FieldGrowthParams(),seed=7521,
                 rng_state=None,electronic_state=None):
        if params.soft_scale<=0 or params.metal_host_contrast<1 or params.crowding_radius_nm<0:
            raise ValueError('Positive supply scale and nonnegative crowding radius required')
        super().__init__(ids,vn,provenance,params,seed,rng_state,electronic_state)
        previous=self.model
        self.model=FieldPilot(seed,params,float(np.asarray(vn).mean()))
        self.model.ids=previous.ids.copy();self.model.vn=previous.vn.copy()
        self.model.injected=previous.injected;self.model.returned=previous.returned
        self.model.mu=params.source_mu_nm;self.model.sigma=params.source_sigma_nm
        self.model.rng.bit_generator.state=copy.deepcopy(previous.rng.bit_generator.state)

    @classmethod
    def from_checkpoint(cls,path,seed=None,params=FieldGrowthParams(),inherit_rng=True):
        return super().from_checkpoint(path,seed=seed,params=params,inherit_rng=inherit_rng)

    def extend(self,attempts=1000):
        for _ in range(int(attempts)):
            factor=supply_factor(np.count_nonzero(self.ids),self.params)
            self.time+=float(self.model.rng.exponential(1.))/factor
            number=len(self.events)+1;before=len(self.model.events)
            accepted=self.model.deposit(self.params.reference_voltage_V,'SET',0,number)
            if len(self.model.events)>before:
                record=copy.deepcopy(self.model.events[-1])
                route=self.model.paths[-1].copy();reservoir=self.model.reservoir_paths[-1].copy()
            else:
                self.null_attempts+=1
                record=dict(event=before,phase='SET',cycle=0,point=number,
                    voltage_V=self.params.reference_voltage_V,kind='no_accessible_source',Ag_id=-1,
                    source_x=-1,source_y=-1,y=-1,x=-1,path_nodes=0,mean_path_us=0.,
                    redox_delta_eV=0.,forward_barrier_eV=0.,reverse_barrier_eV=0.,
                    source_kind='top_connected_Ag_reservoir',delivery_probability=0.,
                    injected=self.model.injected,returned=self.model.returned,
                    field_enhancement=0.,field_weight=0.,local_Ag_fraction=0.,crowding_weight=0.,local_unit_potential=0.)
                route=np.array([],np.int16);reservoir=np.array([],np.int16)
            record.update(supply_factor=factor,attempt=number,cumulative_hazard_clock=self.time,accepted=bool(accepted))
            self.events.append(record);self.paths.append(route);self.reservoir_paths.append(reservoir)
            self._assert_mass()
            self.states.append(self.ids.copy());self.clock.append(self.time)
            self.injected.append(self.model.injected);self.returned.append(self.model.returned)
            self.connected.append(bool(geom.span(self.ids>0).any()))
        return self

    def save(self,path):
        meta=super().save(path)
        meta.update(kernel='field_growth40/v1',params=asdict(self.params),
            field_semantics='Full host-grid quasistatic conductivity surrogate; metal/host contrast; top=1, bottom=0, insulating sides. Field changes candidate probabilities only, distinct from mixed-network terminal I-V.',
            selection='No exp(y/13) depth factor or downward-growth multiplier. Field-enhancement power times local crowding suppression replaces both. Original top-connected supply access and Gaussian lateral source distribution retained.',
            supply='State-dependent attempt rate h(M)=exp(-max(M-soft_volume,0)/soft_scale), or one when disabled. Waiting exposure Exp(1)/h(M); actual material voltage still drives the external clock.',
            clock_units='voltage-driven source-attempt exposure including smooth volume feedback',
            approximation='Geometry field recomputed per attempt at reference bias; existing dielectric transport route kernel is retained.')
        (Path(path)/'restart.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
        return meta

    @classmethod
    def restore(cls,path):
        path=Path(path);meta=json.loads((path/'restart.json').read_text(encoding='utf-8'))
        with np.load(path/'sequence.npz') as z:
            obj=cls(z['metal_id'][0],z['VN'],meta['initial_source'],FieldGrowthParams(**meta['params']),
                    meta['seed'],meta['initial_rng_state'],meta.get('electronic_state',{}))
            obj.states=[s.copy() for s in z['metal_id']];obj.clock=z['cumulative_hazard_clock'].tolist()
            obj.injected=z['injected'].tolist();obj.returned=z['returned'].tolist();obj.connected=z['connected'].tolist()
        obj.model.ids=obj.states[-1].copy();obj.model.injected=meta['injected'];obj.model.returned=meta['returned']
        obj.model.rng.bit_generator.state=meta['rng_state'];obj.time=meta['final_clock'];obj.null_attempts=meta['null_attempts']
        if (path/'events.csv').exists():
            with (path/'events.csv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
            ints={'event','cycle','point','Ag_id','source_x','source_y','y','x','path_nodes','injected','returned','attempt'}
            strings={'phase','kind','source_kind'}
            obj.events=[{k:(v if k in strings else v=='True' if k=='accepted' else int(v) if k in ints else float(v))
                         for k,v in row.items()} for row in rows]
        with np.load(path/'transport.npz') as z:
            obj.paths=[z['nodes'][a:b].copy() for a,b in zip(z['offsets'][:-1],z['offsets'][1:])]
            obj.reservoir_paths=[z['reservoir_nodes'][a:b].copy() for a,b in zip(z['reservoir_offsets'][:-1],z['reservoir_offsets'][1:])]
        obj.model.events=[{k:v for k,v in e.items() if k not in ('attempt','cumulative_hazard_clock','accepted','supply_factor')}
                          for e in obj.events if e['kind']!='no_accessible_source']
        obj._assert_mass()
        return obj


def checkpoint_at_event(source_dir,event_index,destination,observation_exposure=None,electronic_state=None):
    full=FieldGrowthSequence.restore(source_dir);index=int(event_index)
    if not 0<=index<len(full.states):raise ValueError('Event index outside the source sequence')
    part=FieldGrowthSequence(full.states[0],full.vn,full.provenance,full.params,full.seed,
                             full.initial_rng_state,full.initial_electronic_state).extend(index)
    np.testing.assert_array_equal(part.states,full.states[:index+1])
    np.testing.assert_array_equal(part.clock,full.clock[:index+1])
    if observation_exposure is not None:
        if observation_exposure<part.time:raise ValueError('Read exposure precedes selected event')
        if index+1<len(full.clock) and observation_exposure>=full.clock[index+1]:
            raise ValueError('Read exposure includes another event')
    if electronic_state is not None:part.initial_electronic_state=copy.deepcopy(electronic_state)
    meta=part.save(destination)
    read=dict(event_index=index,observation_exposure=observation_exposure,last_event_exposure=part.time,
              next_event_exposure=full.clock[index+1] if index+1<len(full.clock) else None,
              source_sequence=str(Path(source_dir).resolve()),exact_prefix_reproduction=True,
              mass_balance_error=int(np.count_nonzero(part.ids))-part.model.injected+part.model.returned)
    (Path(destination)/'read_state.json').write_text(json.dumps(read,indent=2),encoding='utf-8')
    return meta,read


def self_check(initial):
    empty=np.zeros((80,100),bool)
    phi,enhance=geometry_field(empty)
    np.testing.assert_allclose(phi,np.broadcast_to(1-(np.arange(80)[:,None]+.5)/80,phi.shape),atol=1.e-11)
    np.testing.assert_allclose(enhance[[0,-1]],1.,atol=1.e-10)
    # Deliberately activate supply feedback in this small persistence test.
    p=FieldGrowthParams(soft_volume=0.,soft_scale=500.)
    a=FieldGrowthSequence.from_checkpoint(initial,params=p).extend(8)
    check=validate_routes(a)
    with tempfile.TemporaryDirectory(prefix='field_growth_check_') as temp:
        temp=Path(temp);a.save(temp/'sequence')
        b=FieldGrowthSequence.restore(temp/'sequence')
        a.extend(3);b.extend(3)
        np.testing.assert_array_equal(a.states,b.states);np.testing.assert_array_equal(a.clock,b.clock)
        assert a.events==b.events
        a.save(temp/'sequence')
        checkpoint_at_event(temp/'sequence',5,temp/'checkpoint',a.clock[5])
        c=FieldGrowthSequence.restore(temp/'checkpoint').extend(2)
        np.testing.assert_array_equal(c.states,a.states[:8]);np.testing.assert_array_equal(c.clock,a.clock[:8])
        check.update(restored_continuation_exact=True,checkpoint_prefix_exact=True,
                     empty_host_field_exact=True,minimum_supply_factor=min(e['supply_factor'] for e in a.events))
    return check


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check',type=Path,required=True,help='Initial saved checkpoint for the persistence/field check')
    args=ap.parse_args()
    print(json.dumps(self_check(args.check),indent=2))
