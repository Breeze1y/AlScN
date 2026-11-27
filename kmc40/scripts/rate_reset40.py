"""Source-conserving RESET erosion with an absolute, separable hazard clock.

The event catalogue has rates r_i(V)=a(V)*w_i(state).  w_i includes a paired
oxidation barrier and present-state mixed-network dissipation.  It is never
normalized before computing the total hazard.  A voltage exposure replay can
therefore fit a(V) without resampling which Ag packets leave the device.

Fields are evaluated at unit bias and fixed electronic activity.  This is a
separable mesoscopic approximation, not a self-consistent thermal KMC model.
"""
from dataclasses import asdict, dataclass
from pathlib import Path
import argparse
import copy
import csv
import hashlib
import json
import time

import numpy as np
from scipy.ndimage import convolve
from scipy.sparse.csgraph import dijkstra

import source_cycle_vn40 as source
import form_hourglass40 as geom
import cycle_hourglass40 as oldcycle

ROOT = Path(__file__).resolve().parents[1]
INITIAL = ROOT / 'results/vacancy_assisted40/source_cycle/assisted/states.npz'
OUT = ROOT / 'results/reset_unified40/kinetic_pilot'


@dataclass(frozen=True)
class RateParams:
    eta: float = .03
    assist_eV: float = .05
    kT_eV: float = .025852
    reference_barrier_eV: float = .30
    unbound_energy_eV: float = .075
    coordination_eV: float = .035
    adhesion_eV: float = .08
    background_weight: float = .002
    power_exponent: float = 1.0
    path_reference_voltage_V: float = -.6


@dataclass(frozen=True)
class ClockParams:
    rate_scale: float = .03
    voltage_scale_V: float = .30
    dwell_units: float = 1.0


def load_initial(path=INITIAL, cycle=1):
    """Recover an actual SET return-to-zero state, preserving every Ag ID."""
    path = Path(path)
    with np.load(path) as z:
        pick = np.flatnonzero((z['phase'].astype(str) == 'SET') & (z['cycle'] == cycle))[-1]
        ids, vn = z['metal_id'][pick].copy(), z['VN'].copy()
        injected, returned = int(z['injected'][pick]), int(z['returned'][pick])
    assert (ids > 0).sum() == injected - returned
    provenance = dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      frame=int(pick), cycle=cycle, injected=injected, returned=returned)
    return ids, vn, provenance


def common_exposure(voltage, params=ClockParams()):
    if params.rate_scale < 0 or params.voltage_scale_V <= 0 or params.dwell_units <= 0:
        raise ValueError('Clock parameters must be nonnegative/positive.')
    v = np.asarray(voltage, float)
    if np.any(~np.isfinite(v)):
        raise ValueError('Non-finite voltage')
    return np.cumsum(params.rate_scale * params.dwell_units *
                     np.expm1(np.abs(v) / params.voltage_scale_V))


def replay_indices(cumulative_hazard_clock, voltage, params=ClockParams()):
    times = np.asarray(cumulative_hazard_clock, float)
    if times.ndim != 1 or len(times) == 0 or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError('Sequence clock must start at zero and increase strictly.')
    exposure = common_exposure(voltage, params)
    indices = np.searchsorted(times, exposure, side='right') - 1
    return np.clip(indices, 0, len(times) - 1), exposure


def barrier_and_weights(metal, vn, power, reference_power, params=RateParams()):
    """The same vacancy energy changes both site preference and total rate."""
    cn = convolve(metal.astype(int), geom.KERNEL, mode='constant')
    surface = metal & (cn < 4)
    dg = params.unbound_energy_eV - params.coordination_eV * cn
    dg = dg + params.eta * (-.5 - 2.48) * (vn & (cn > 0))
    dg[[0, -1]] -= params.adhesion_eV
    barrier = params.reference_barrier_eV + np.maximum(-dg, 0.)
    # Reference is frozen once at the SET state, unlike per-frame maximum
    # normalization.  Actual power lost after a break lowers the event hazard.
    relative_power = np.maximum(np.nan_to_num(power), 0.) / max(reference_power, 1e-30)
    drive = params.background_weight + relative_power ** params.power_exponent
    weights = surface * np.exp(-(barrier - params.reference_barrier_eV) / params.kT_eV) * drive
    return barrier, weights


def return_path(metal, target, vn, params):
    """Compute a continuous dielectric return route to the top Ag reservoir."""
    y, x = divmod(int(target), metal.shape[1])
    opened = metal.copy()
    opened[y, x] = False
    graph, _, _ = source.source_routes(opened, params.path_reference_voltage_V,
                                       vn, params.eta, params.assist_eV)
    dist, pred = dijkstra(graph, indices=int(target), return_predecessors=True)
    top = oldcycle.clusters(opened)[0]
    destinations = (~opened) & (convolve(top.astype(int), geom.KERNEL, mode='constant') > 0)
    destinations[0] |= ~opened[0]
    destinations[y, x] = False
    valid = np.flatnonzero(np.isfinite(dist) & destinations.ravel())
    if not len(valid):
        return None, None
    end = int(valid[np.argmin(dist[valid])])
    route = [end]
    while route[-1] != target:
        node = int(pred[route[-1]])
        if node < 0:
            raise RuntimeError('Broken return route')
        route.append(node)
    route = np.array(route[::-1], np.int16)
    assert not opened.ravel()[route].any()
    yy, xx = np.divmod(route.astype(int), metal.shape[1])
    assert np.all(np.abs(np.diff(yy)) + np.abs(np.diff(xx)) == 1)
    return route, float(dist[end])


def mixed_field_solver(network_params=None, activity=1.):
    import mixed_network40 as net
    network_params = network_params or net.NetworkParams()
    def field(metal, vn):
        return net.solve(net.prepare(metal, vn), network_params, activity=activity, fields=True)
    field.description = dict(kind='mixed_network40', params=asdict(network_params), activity=activity)
    return field


def legacy_field_solver(metal, vn):
    power, conductance = geom.current_field(metal)
    return dict(unit_power=power, conductance_S=.04 * conductance,
                unit_phi=np.zeros(metal.shape))


class ErosionSequence:
    def __init__(self, ids, vn, provenance, params=RateParams(), seed=7421, field_solver=None,
                 rng_state=None, reference_power=None):
        self.ids = ids.copy()
        self.vn = vn.copy()
        self.provenance = provenance
        self.params = params
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        if rng_state is not None:
            self.rng.bit_generator.state = copy.deepcopy(rng_state)
        self.initial_rng_state = copy.deepcopy(self.rng.bit_generator.state)
        self.field_solver = field_solver or mixed_field_solver()
        self.solver_description = getattr(self.field_solver, 'description', {'kind': self.field_solver.__name__})
        self.injected = int(provenance['injected'])
        self.returned = int(provenance['returned'])
        self.events = []
        self.paths = []
        self.states = [self.ids.copy()]
        self.clock = [0.]
        self.null_events = 0
        self.time = 0.
        self.field = self.field_solver(self.ids > 0, self.vn)
        self.reference_power = (float(np.max(np.nan_to_num(self.field['unit_power'])))
                                if reference_power is None else float(reference_power))
        if not np.isfinite(self.reference_power) or self.reference_power < 0:
            raise ValueError('reference_power must be finite and nonnegative')
        self.initial_reference_power = self.reference_power
        self.conductance = [float(self.field['conductance_S'])]
        self._assert_mass()

    def _assert_mass(self):
        assert (self.ids > 0).sum() == self.injected - self.returned
        positive = self.ids[self.ids > 0]
        assert len(positive) == len(np.unique(positive))

    def site_weights(self, metal):
        return barrier_and_weights(metal, self.vn, self.field['unit_power'],
                                   self.reference_power, self.params)

    def extend(self, events=300, max_null_attempts=5000):
        target_count = len(self.events) + int(events)
        attempts = 0
        while len(self.events) < target_count:
            metal = self.ids > 0
            barrier, weights = self.site_weights(metal)
            total = float(weights.sum())
            if total <= 0:
                break
            # Exact thinning: impossible-return attempts consume waiting time
            # but leave state unchanged, rather than silently renormalizing.
            wait = float(self.rng.exponential(1. / total))
            self.time += wait
            node = int(self.rng.choice(weights.size, p=weights.ravel() / total))
            route, travel = return_path(metal, node, self.vn, self.params)
            if route is None:
                self.null_events += 1
                attempts += 1
                if attempts > max_null_attempts:
                    raise RuntimeError('Too many inaccessible surface events')
                continue
            y, x = divmod(node, metal.shape[1])
            agid = int(self.ids[y, x])
            self.ids[y, x] = 0
            self.returned += 1
            self._assert_mass()
            self.events.append(dict(event=len(self.events) + 1, Ag_id=agid, y=y, x=x,
                cumulative_hazard_clock=self.time, attempted_wait=wait, total_weight=total,
                site_weight=float(weights[y, x]), oxidation_barrier_eV=float(barrier[y, x]),
                unit_power=float(self.field['unit_power'][y, x]),
                adjacent_vacancy=bool(self.vn[y, x]), route_nodes=len(route),
                reference_path_mean_us=travel, injected=self.injected, returned=self.returned))
            self.paths.append(route)
            self.states.append(self.ids.copy())
            self.clock.append(self.time)
            self.field = self.field_solver(self.ids > 0, self.vn)
            self.conductance.append(float(self.field['conductance_S']))
        return self

    def save(self, path):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path / 'sequence.npz', metal_id=np.asarray(self.states), VN=self.vn,
            cumulative_hazard_clock=self.clock, conductance_S=self.conductance,
            total_weight=np.r_[np.nan, [e['total_weight'] for e in self.events]],
            initial_injected=self.injected, initial_returned=self.provenance['returned'])
        np.savez_compressed(path / 'transport.npz', offsets=np.r_[0, np.cumsum([len(p) for p in self.paths])],
            nodes=np.concatenate(self.paths) if self.paths else np.array([], np.int16),
            Ag_id=np.array([e['Ag_id'] for e in self.events], dtype=np.int32))
        if self.events:
            with (path / 'events.csv').open('w', encoding='utf8', newline='') as f:
                w = csv.DictWriter(f, fieldnames=self.events[0]); w.writeheader(); w.writerows(self.events)
        metadata = dict(params=asdict(self.params), seed=self.seed, initial_source=self.provenance,
            field_solver=self.solver_description, reference_power=self.reference_power,
            initial_reference_power=self.initial_reference_power, initial_rng_state=self.initial_rng_state,
            events=len(self.events), null_events=self.null_events, final_clock=self.time,
            initial_Ag=int(np.count_nonzero(self.states[0])), final_Ag=int(np.count_nonzero(self.ids)),
            injected=self.injected, returned=self.returned, rng_state=self.rng.bit_generator.state,
            clock_units='dimensionless integrated voltage exposure; dwell in sampling-point units',
            approximation='Separable rates; unit-voltage mixed fields at fixed electronic activity; no thermal memory',
            return_paths='Source-conserving continuous dielectric paths at a fixed reference bias')
        (path / 'restart.json').write_text(json.dumps(metadata, indent=2), encoding='utf8')
        np.savez_compressed(path / 'restart.npz', metal_id=self.ids, VN=self.vn)
        return metadata

    @classmethod
    def restore(cls, path, field_solver=None):
        path = Path(path)
        meta = json.loads((path / 'restart.json').read_text(encoding='utf8'))
        if field_solver is None:
            description = meta['field_solver']
            if description['kind'] == 'mixed_network40':
                import mixed_network40 as net
                field_solver = mixed_field_solver(net.NetworkParams(**description['params']), description['activity'])
            elif description['kind'] == 'legacy_field_solver':
                field_solver = legacy_field_solver
            else:
                raise ValueError('Supply the custom field solver to restore this run')
        with np.load(path / 'restart.npz') as z:
            provenance = dict(meta['initial_source'], injected=meta['injected'], returned=meta['returned'])
            obj = cls(z['metal_id'], z['VN'], provenance, RateParams(**meta['params']), meta['seed'], field_solver,
                      rng_state=meta.get('initial_rng_state'),
                      reference_power=meta.get('initial_reference_power', meta['reference_power']))
        with np.load(path / 'sequence.npz') as z:
            obj.states = [state.copy() for state in z['metal_id']]
            obj.clock = z['cumulative_hazard_clock'].tolist()
            obj.conductance = z['conductance_S'].tolist()
        if (path / 'events.csv').exists():
            with (path / 'events.csv').open(encoding='utf8', newline='') as f:
                rows = list(csv.DictReader(f))
            integer_keys = ('event', 'Ag_id', 'y', 'x', 'route_nodes', 'injected', 'returned')
            obj.events = [{k: int(v) if k in integer_keys else v == 'True' if k == 'adjacent_vacancy'
                           else float(v) for k, v in row.items()} for row in rows]
        with np.load(path / 'transport.npz') as z:
            obj.paths = [z['nodes'][a:b].copy() for a, b in zip(z['offsets'][:-1], z['offsets'][1:])]
        obj.rng.bit_generator.state = meta['rng_state']
        obj.time = meta['final_clock']
        obj.reference_power = meta['reference_power']
        obj.null_events = meta['null_events']
        obj.provenance = meta['initial_source']
        return obj


def checkpoint_at_event(source_dir, event_index, destination, observation_exposure=None):
    """Save the actual fitted read state and a deterministic continuation RNG.

    A fitted sampling point may lie between two events.  Its material state is
    that of event_index, while observation_exposure records the accumulated
    drive.  The generator retains the last event's clock/RNG so its next draw
    reproduces the original future event, including the waiting interval.
    """
    full = ErosionSequence.restore(source_dir)
    event_index = int(event_index)
    if event_index < 0 or event_index >= len(full.states):
        raise ValueError('Event index outside the stored sequence')
    provenance = full.provenance
    replayed = ErosionSequence(full.states[0], full.vn, provenance, full.params,
                              full.seed, full.field_solver, rng_state=full.initial_rng_state,
                              reference_power=full.initial_reference_power).extend(event_index)
    np.testing.assert_array_equal(replayed.states, full.states[:event_index + 1])
    np.testing.assert_allclose(replayed.clock, full.clock[:event_index + 1], rtol=0, atol=0)
    if observation_exposure is not None:
        if observation_exposure < replayed.time:
            raise ValueError('Observation precedes the selected event')
        if event_index + 1 < len(full.clock) and observation_exposure >= full.clock[event_index + 1]:
            raise ValueError('Observation includes an event after the selected state')
    meta = replayed.save(destination)
    read_state = dict(source_sequence=str(Path(source_dir).resolve()), event_index=event_index,
        last_event_exposure=replayed.time, observation_exposure=observation_exposure,
        next_event_exposure=full.clock[event_index + 1] if event_index + 1 < len(full.clock) else None,
        Ag_packets=int((replayed.ids > 0).sum()), mass_balance_error=int((replayed.ids > 0).sum())-
            (replayed.injected-replayed.returned),
        continuation='Restore restart.json/npz and keep the original accumulated exposure when replaying voltage.',
        exact_prefix_reproduction=True)
    (Path(destination) / 'read_state.json').write_text(json.dumps(read_state, indent=2), encoding='utf8')
    return meta, read_state


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--events', type=int, default=300)
    p.add_argument('--seed', type=int, default=7421)
    p.add_argument('--cycle', type=int, default=1, help='SET cycle in the inherited source trajectory')
    p.add_argument('--network-params', type=Path, help='JSON NetworkParams for mixed-network fields')
    p.add_argument('--out', type=Path, default=OUT)
    p.add_argument('--legacy-fields', action='store_true', help='Explicit diagnostic only')
    args = p.parse_args()
    start = time.perf_counter()
    ids, vn, provenance = load_initial(cycle=args.cycle)
    field_solver = legacy_field_solver if args.legacy_fields else None
    if args.network_params:
        if args.legacy_fields:
            raise ValueError('Network parameters require mixed-network fields')
        import mixed_network40 as net
        network_params = net.NetworkParams(**json.loads(args.network_params.read_text(encoding='utf8')))
        field_solver = mixed_field_solver(network_params)
    sequence = ErosionSequence(ids, vn, provenance, seed=args.seed,
                               field_solver=field_solver).extend(args.events)
    meta = sequence.save(args.out)
    bridges = np.array([bool(geom.span(x > 0).any()) for x in sequence.states])
    first_break = np.flatnonzero(~bridges)
    print(json.dumps(dict(events=len(sequence.events), elapsed_s=time.perf_counter() - start,
        initial_Ag=meta['initial_Ag'], final_Ag=meta['final_Ag'],
        first_break_event=int(first_break[0]) if len(first_break) else None,
        final_clock=sequence.clock[-1]), indent=2), flush=True)


if __name__ == '__main__':
    main()
