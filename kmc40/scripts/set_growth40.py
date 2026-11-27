"""Top-source-conserving SET growth sequences with a separable attempt clock.

Every attempted packet follows the existing top-connected reservoir and host
routes. The spatial reaction distribution is evaluated at a fixed reference
bias. Applied/device voltage and compliance determine the attempt-clock speed
in the continuous-cycle driver, not in this reusable geometry sequence.
"""
from dataclasses import asdict, dataclass
from pathlib import Path
import argparse, copy, csv, hashlib, json, time
import numpy as np

import source_cycle_vn40 as source
import form_hourglass40 as geom

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_INITIAL=ROOT/'results/reset_unified40/retention/HRS'
OUT=ROOT/'results/cycle_unified40/growth_initial'


@dataclass(frozen=True)
class GrowthParams:
    reference_voltage_V: float=.33
    eta: float=.03
    assist_eV: float=.05
    pt_mode: str='top_front'
    source_mu_nm: float=25.
    source_sigma_nm: float=6.


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class GrowthSequence:
    def __init__(self, ids, vn, provenance, params=GrowthParams(), seed=7521,
                 rng_state=None, electronic_state=None):
        ids=np.asarray(ids,dtype=np.int32);vn=np.asarray(vn,dtype=bool)
        if ids.shape!=(80,100) or vn.shape!=ids.shape:
            raise ValueError('The current growth kernel expects an 80x100 lattice')
        self.params=params;self.seed=int(seed);self.provenance=copy.deepcopy(provenance)
        self.initial_electronic_state=copy.deepcopy(electronic_state or {})
        self.model=source.Pilot(self.seed,params.eta,params.assist_eV,params.pt_mode,
                                float(vn.mean()))
        self.model.ids=ids.copy();self.model.vn=vn.copy()
        self.model.injected=int(provenance['injected']);self.model.returned=int(provenance['returned'])
        self.model.mu=params.source_mu_nm;self.model.sigma=params.source_sigma_nm
        if rng_state is not None:self.model.rng.bit_generator.state=copy.deepcopy(rng_state)
        self.initial_rng_state=copy.deepcopy(self.model.rng.bit_generator.state)
        self.states=[ids.copy()];self.clock=[0.];self.events=[]
        self.injected=[self.model.injected];self.returned=[self.model.returned]
        self.connected=[bool(geom.span(ids>0).any())]
        self.paths=[];self.reservoir_paths=[]
        self.time=0.;self.null_attempts=0
        self._assert_mass()

    @classmethod
    def from_arrays(cls, ids, vn, provenance, **kwargs):
        return cls(ids,vn,provenance,**kwargs)

    @classmethod
    def from_checkpoint(cls, path=DEFAULT_INITIAL, seed=None, params=GrowthParams(), inherit_rng=True):
        path=Path(path);meta=json.loads((path/'restart.json').read_text(encoding='utf-8'))
        with np.load(path/'restart.npz') as z:
            ids=z['metal_id'].copy();vn=z['VN'].copy()
        provenance=dict(path=str(path.resolve()),restart_sha256=_hash(path/'restart.npz'),
                        injected=int(meta['injected']),returned=int(meta['returned']))
        return cls(ids,vn,provenance,params,meta.get('seed',7521) if seed is None else seed,
                   meta['rng_state'] if inherit_rng and seed is None else None,
                   meta.get('electronic_state',{}))

    @property
    def ids(self):return self.model.ids

    @property
    def vn(self):return self.model.vn

    def _assert_mass(self):
        positive=self.ids[self.ids>0]
        assert len(positive)==len(np.unique(positive))
        assert len(positive)==self.model.injected-self.model.returned
        assert positive.max(initial=0)<=self.model.injected

    def extend(self, attempts=1000):
        for _ in range(int(attempts)):
            self.time+=float(self.model.rng.exponential(1.))
            number=len(self.events)+1
            before=len(self.model.events)
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
                    injected=self.model.injected,returned=self.model.returned)
                route=np.array([],np.int16);reservoir=np.array([],np.int16)
            record.update(attempt=number,cumulative_hazard_clock=self.time,accepted=bool(accepted))
            self.events.append(record);self.paths.append(route);self.reservoir_paths.append(reservoir)
            self._assert_mass()
            self.states.append(self.ids.copy());self.clock.append(self.time)
            self.injected.append(self.model.injected);self.returned.append(self.model.returned)
            self.connected.append(bool(geom.span(self.ids>0).any()))
        return self

    def save(self,path):
        path=Path(path);path.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(path/'sequence.npz',metal_id=np.asarray(self.states),VN=self.vn,
            cumulative_hazard_clock=self.clock,injected=self.injected,returned=self.returned,
            connected=self.connected,initial_injected=self.injected[0],initial_returned=self.returned[0])
        np.savez_compressed(path/'transport.npz',
            offsets=np.r_[0,np.cumsum([len(p) for p in self.paths])],
            nodes=np.concatenate(self.paths) if self.paths else np.array([],np.int16),
            reservoir_offsets=np.r_[0,np.cumsum([len(p) for p in self.reservoir_paths])],
            reservoir_nodes=np.concatenate(self.reservoir_paths) if self.reservoir_paths else np.array([],np.int16))
        if self.events:
            with (path/'events.csv').open('w',encoding='utf-8',newline='') as f:
                w=csv.DictWriter(f,fieldnames=self.events[0].keys());w.writeheader();w.writerows(self.events)
        split=np.flatnonzero(self.connected)
        meta=dict(params=asdict(self.params),seed=self.seed,initial_source=self.provenance,
            initial_rng_state=self.initial_rng_state,rng_state=self.model.rng.bit_generator.state,
            electronic_state=self.initial_electronic_state,
            attempts=len(self.events),successful_depositions=int(np.count_nonzero(self.ids)-np.count_nonzero(self.states[0])),
            null_attempts=self.null_attempts,injected=self.model.injected,returned=self.model.returned,
            initial_Ag=int(np.count_nonzero(self.states[0])),final_Ag=int(np.count_nonzero(self.ids)),
            final_clock=self.time,first_bridge_event=int(split[0]) if len(split) else None,
            clock_units='unit-rate source-attempt exposure; exponential thresholds include unreacted returns',
            approximation='Fixed reference-voltage spatial growth probabilities; actual voltage/compliance scales the attempt clock externally.',
            electronic_state_note='Inherited initial electronic state only; geometry precomputation does not evolve q. Continuous-cycle driver owns q.',
            energies_eV=dict(unbound=.075,coordination=.035,adhesion=.08,redox_reference=.30,
                             vacancy_pair_effect=params_pair(self.params)))
        (path/'restart.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
        np.savez_compressed(path/'restart.npz',metal_id=self.ids,VN=self.vn)
        return meta

    @classmethod
    def restore(cls,path):
        path=Path(path);meta=json.loads((path/'restart.json').read_text(encoding='utf-8'))
        with np.load(path/'sequence.npz') as z:
            obj=cls(z['metal_id'][0],z['VN'],meta['initial_source'],GrowthParams(**meta['params']),
                    meta['seed'],meta['initial_rng_state'],meta.get('electronic_state',{}))
            obj.states=[s.copy() for s in z['metal_id']];obj.clock=z['cumulative_hazard_clock'].tolist()
            obj.injected=z['injected'].tolist();obj.returned=z['returned'].tolist()
            obj.connected=z['connected'].tolist()
        obj.model.ids=obj.states[-1].copy();obj.model.injected=meta['injected'];obj.model.returned=meta['returned']
        obj.model.rng.bit_generator.state=meta['rng_state'];obj.time=meta['final_clock']
        obj.null_attempts=meta['null_attempts']
        if (path/'events.csv').exists():
            with (path/'events.csv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
            ints={'event','cycle','point','Ag_id','source_x','source_y','y','x','path_nodes','injected','returned','attempt'}
            strings={'phase','kind','source_kind'}
            obj.events=[{k:(v if k in strings else v=='True' if k=='accepted' else int(v) if k in ints else float(v))
                         for k,v in row.items()} for row in rows]
        with np.load(path/'transport.npz') as z:
            obj.paths=[z['nodes'][a:b].copy() for a,b in zip(z['offsets'][:-1],z['offsets'][1:])]
            obj.reservoir_paths=[z['reservoir_nodes'][a:b].copy() for a,b in zip(z['reservoir_offsets'][:-1],z['reservoir_offsets'][1:])]
        # The source kernel's event counter is diagnostic. Keep it continuous,
        # while its live path arrays only need to hold newly generated routes.
        obj.model.events=[{k:v for k,v in e.items() if k not in ('attempt','cumulative_hazard_clock','accepted')}
                          for e in obj.events if e['kind']!='no_accessible_source']
        obj._assert_mass()
        return obj


def params_pair(params):return params.eta*(-.5-2.48)


def checkpoint_at_event(source_dir,event_index,destination,observation_exposure=None,electronic_state=None):
    full=GrowthSequence.restore(source_dir);index=int(event_index)
    if not 0<=index<len(full.states):raise ValueError('Event index outside the source sequence')
    part=GrowthSequence(full.states[0],full.vn,full.provenance,full.params,full.seed,
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


def validate_routes(sequence):
    nx=100;checked=0
    for i,(e,route,reservoir) in enumerate(zip(sequence.events,sequence.paths,sequence.reservoir_paths)):
        before=sequence.states[i];after=sequence.states[i+1]
        if e['kind']=='no_accessible_source':continue
        assert len(route)>0 and len(reservoir)>0
        assert not np.any(before.ravel()[route]>0)
        assert reservoir[0]//nx==0 and reservoir[-1]==route[0]
        assert np.all(before.ravel()[reservoir[:-1]]>0)
        for path in [route,reservoir]:
            y,x=np.divmod(path.astype(int),nx)
            assert np.all(abs(np.diff(y))+abs(np.diff(x))==1)
        if e['accepted']:
            assert np.count_nonzero(after!=before)==1
            assert after[e['y'],e['x']]==e['Ag_id']
        else:np.testing.assert_array_equal(before,after)
        assert np.count_nonzero(after)==e['injected']-e['returned']
        checked+=1
    return dict(attempts_checked=checked,source_and_host_routes_contiguous=True,
                unique_source_Ag_ids_and_mass_balance=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--initial',type=Path,default=DEFAULT_INITIAL)
    p.add_argument('--events',type=int,default=1000);p.add_argument('--seed',type=int)
    p.add_argument('--out',type=Path,default=OUT);args=p.parse_args()
    start=time.perf_counter();seq=GrowthSequence.from_checkpoint(args.initial,seed=args.seed).extend(args.events)
    meta=seq.save(args.out);check=validate_routes(seq)
    bridge=meta['first_bridge_event']
    if bridge and bridge>0:
        at=seq.events[bridge-1];check['first_bridge_site']=[at['y'],at['x']]
        check['first_bridge_depth_nm']=(at['y']+.5)*.5
    check.update(elapsed_s=time.perf_counter()-start,first_bridge_event=bridge,
                 initial_Ag=meta['initial_Ag'],final_Ag=meta['final_Ag'])
    (args.out/'validation.json').write_text(json.dumps(check,indent=2),encoding='utf-8')
    print(json.dumps(check,indent=2),flush=True)

if __name__=='__main__':main()
