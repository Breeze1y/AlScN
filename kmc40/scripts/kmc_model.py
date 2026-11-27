"""Morphology-first, coarse-grained Ag/AlN/AlScN/AlN/Pt CT-KMC.

All energies eV, potentials V, distances nm. One occupied voxel is an
effective Ag packet, NOT one resolved atom. No prescribed growth path.
"""
from dataclasses import dataclass, asdict, replace
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy import sparse, ndimage
from scipy.sparse.linalg import spsolve

KB = 8.617333262e-5
EPS0 = 8.8541878128e-12
CROSS = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
CHECKPOINT_SCHEMA = 2


@dataclass
class Config:
    nx: int = 80
    layer_cells: int = 20
    spacing_nm: float = 1.0
    voltage: float = 3.0
    temperature_K: float = 300.0
    attempt_Hz: float = 1e13
    eps_aln: float = 8.5
    eps_alscn: float = 15.0
    polarization_C_m2: float = -0.004
    hop_aln_eV: float = 0.48
    hop_alscn_eV: float = 0.50
    interface_eV: float = 0.025
    source_barrier_eV: float = 0.68
    reaction_barrier_eV: float = 0.52
    pt_nucleation_eV: float = 0.08
    metal_standard_eV: float = -0.10
    bond_eV: float = 0.06
    pt_adhesion_eV: float = 0.06
    surface_redox_only: bool = True
    migration_eV: float = 0.80
    disorder_eV: float = 0.025
    max_events: int = 350000
    max_time_s: float = 2.0
    seed: int = 7
    disorder_seed: int | None = None
    record_events: bool = False
    vacancy_affected_fraction_by_layer: tuple[float, ...] = (0.0, 0.0, 0.0)
    vacancy_ag_ion_penalty_eV: float = 0.0
    vacancy_seed: int = 0
    # Explicit Ag-to-Pt material/cell-count pairs. None preserves the historical
    # AlN/AlScN/AlN stack with layer_cells cells in each of the three layers.
    stack: tuple[tuple[str, int], ...] | None = None

    def __post_init__(self):
        for name in ('nx', 'layer_cells'):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value <= 0:
                raise ValueError(f'{name} must be a positive integer')
            setattr(self, name, int(value))
        if not np.isfinite(self.spacing_nm) or self.spacing_nm <= 0:
            raise ValueError('spacing_nm must be finite and positive')
        if self.stack is not None:
            if not isinstance(self.stack, (list, tuple)) or not self.stack:
                raise ValueError('stack must contain Ag-to-Pt [material, cell_count] pairs')
            normalized = []
            for entry in self.stack:
                if not isinstance(entry, (list, tuple)) or len(entry) != 2:
                    raise ValueError('Each stack layer must be [material, cell_count]')
                material, cells = entry
                if material not in ('AlN', 'AlScN'):
                    raise ValueError('Stack material must be AlN or AlScN')
                if isinstance(cells, bool) or not isinstance(cells, (int, np.integer)) or cells <= 0:
                    raise ValueError('Stack cell counts must be positive integers')
                normalized.append((material, int(cells)))
            self.stack = tuple(normalized)
        n_layers = len(self.resolved_stack)
        fractions = self.vacancy_affected_fraction_by_layer
        if not isinstance(fractions, (list, tuple)):
            raise ValueError('Specify vacancy-affected fractions in Ag-to-Pt layer order')
        fractions = tuple(float(fraction) for fraction in fractions)
        if len(fractions) != n_layers:
            # Permit the legacy all-zero default with an explicit new stack,
            # but never silently transplant nonzero middle-layer defects.
            if self.stack is not None and fractions == (0.0, 0.0, 0.0):
                fractions = (0.0,) * n_layers
            else:
                raise ValueError('Specify one vacancy-affected fraction per physical layer')
        if not all(0.0 <= fraction <= 1.0 for fraction in fractions):
            raise ValueError('Vacancy-affected fractions must be finite and in [0, 1]')
        penalty = float(self.vacancy_ag_ion_penalty_eV)
        if not np.isfinite(penalty) or penalty < 0.0:
            raise ValueError('Ag ion vacancy penalty must be finite and nonnegative')
        if (isinstance(self.vacancy_seed, bool) or
                not isinstance(self.vacancy_seed, (int, np.integer)) or self.vacancy_seed < 0):
            raise ValueError('Vacancy seed must be a nonnegative integer')
        self.vacancy_affected_fraction_by_layer = fractions
        self.vacancy_ag_ion_penalty_eV = penalty
        self.vacancy_seed = int(self.vacancy_seed)

    @property
    def resolved_stack(self):
        if self.stack is not None:
            return self.stack
        return (('AlN', self.layer_cells), ('AlScN', self.layer_cells), ('AlN', self.layer_cells))


class Model:
    def __init__(self, cfg):
        self.c = cfg
        self.rng = np.random.default_rng(cfg.seed)
        self.layer_materials = tuple(material for material, _ in cfg.resolved_stack)
        self.layer_cell_counts = tuple(cells for _, cells in cfg.resolved_stack)
        self.n_layers = len(self.layer_cell_counts)
        self.nz = sum(self.layer_cell_counts) + 2
        self.shape = (self.nz, cfg.nx)
        self.n = self.nz * cfg.nx
        self.state = np.zeros(self.n, dtype=np.uint8)  # 0 host, 1 ion, 2 metal
        self.state[:cfg.nx] = 2
        self.state[-cfg.nx:] = 2
        self.interior = np.ones(self.n, bool)
        self.interior[:cfg.nx] = self.interior[-cfg.nx:] = False
        grid = np.arange(self.n).reshape(self.shape)
        self.neighbors = np.full((self.n, 4), -1, int)
        self.neighbors[grid[1:].ravel(), 0] = grid[:-1].ravel()
        self.neighbors[grid[:-1].ravel(), 1] = grid[1:].ravel()
        self.neighbors[grid[:, 1:].ravel(), 2] = grid[:, :-1].ravel()
        self.neighbors[grid[:, :-1].ravel(), 3] = grid[:, 1:].ravel()
        self.safe = np.maximum(self.neighbors, 0)
        self.valid = self.neighbors >= 0
        self.layer = np.zeros(self.n, int)
        # layer is the material code used by field/rate solvers; zone is the
        # physical layer index used only for inventory and residence totals.
        # Keep electrode material tags at zero, as in the archived baseline.
        self.zone = np.zeros(self.n, int)
        row = 1
        for z, (material, cells) in enumerate(cfg.resolved_stack):
            self.layer.reshape(self.shape)[row:row+cells] = int(material == 'AlScN')
            self.zone.reshape(self.shape)[row:row+cells] = z
            row += cells
        self.zone[-cfg.nx:] = self.n_layers - 1
        disorder_rng = self.rng if cfg.disorder_seed is None else np.random.default_rng(cfg.disorder_seed)
        self.site = disorder_rng.normal(0, cfg.disorder_eV, self.n)
        # Static host annotation independent of Ag occupancy and the event RNG.
        # The effective ion energy is shared by injection, hopping and redox.
        self.vacancy_mask = np.zeros(self.n, dtype=bool)
        fractions = np.asarray(cfg.vacancy_affected_fraction_by_layer)
        if fractions.any():
            probability = fractions[self.zone[self.interior]]
            vacancy_rng = np.random.default_rng(cfg.vacancy_seed)
            self.vacancy_mask[self.interior] = vacancy_rng.random(probability.size) < probability
        self.vacancy_mask.setflags(write=False)
        self._vacancy_indices = np.flatnonzero(self.vacancy_mask)
        if cfg.vacancy_ag_ion_penalty_eV:
            self.site[self._vacancy_indices] += cfg.vacancy_ag_ion_penalty_eV
        self.vacancy_affected_sites_by_layer = np.bincount(self.zone[self._vacancy_indices], minlength=self.n_layers)
        self.time = 0.0
        self.steps = 0
        self.net_source = 0
        self.counts = {k: 0 for k in ('inject', 'absorb', 'hop', 'reduce', 'oxidize', 'migrate')}
        self.snapshots = []
        self.trace = []
        self.pending_event = None
        self.snapshot_target = 25
        self.event_log = []
        self.protocol_position = None
        self.hop_layer_counts = np.zeros((self.n_layers,self.n_layers),dtype=np.int64)
        self.net_deposition = np.zeros(self.n_layers,dtype=np.int64)
        self.ion_residence_s = np.zeros(self.n_layers,dtype=float)
        self.vacancy_ion_residence_s = np.zeros(self.n_layers,dtype=float)
        self.update_structure()
        self.save()

    def component(self, side):
        labels, _ = ndimage.label((self.state == 2).reshape(self.shape), CROSS)
        return (labels == labels[side, 0]).ravel()

    def update_structure(self):
        self.cathode = self.component(-1)
        self.anode = self.component(0)
        self.connected = bool(self.cathode[:self.c.nx].any())
        metal_n = (self.state[self.safe] == 2) & self.valid
        self.coord = metal_n.sum(axis=1)
        if not self.connected:
            self.phi = self.solve_field()

    def solve_field(self):
        """Finite-volume div(eps grad(phi)) = div(P).

        Normal face P/eps is distance-weighted. Lateral walls have zero
        normal flux. Cathode-connected Ag is fixed at Pt potential.
        Floating islands are high-eps inclusions (approximation).
        """
        c = self.c
        eps = EPS0 * np.where(self.layer == 1, c.eps_alscn, c.eps_aln)
        p = np.where(self.layer == 1, c.polarization_C_m2, 0.0)
        eps[self.state == 2] = EPS0 * 1e5
        p[self.state == 2] = 0
        fixed = self.cathode | self.anode
        values = np.where(self.anode, c.voltage, 0.0)
        unknown = np.flatnonzero(~fixed)
        index = np.full(self.n, -1, int)
        index[unknown] = np.arange(len(unknown))
        rows, cols, data = [], [], []
        rhs = np.zeros(len(unknown))
        diag = np.zeros(len(unknown))
        h = c.spacing_nm * 1e-9
        for direction in range(4):
            active = self.valid[unknown, direction]
            i = unknown[active]
            j = self.neighbors[i, direction]
            r = index[i]
            g = 2 * eps[i] * eps[j] / (eps[i] + eps[j]) / h**2
            diag[r] += g
            f = fixed[j]
            rhs[r[f]] += g[f] * values[j[f]]
            if direction < 2:
                pf = (eps[j]*p[i] + eps[i]*p[j]) / (eps[i]+eps[j])
                # outward displacement: eps(phi_i-phi_j)/h + P_out
                rhs[r] -= (1 if direction == 1 else -1) * pf / h
            rows.append(r[~f]); cols.append(index[j[~f]]); data.append(-g[~f])
        rows.append(np.arange(len(unknown))); cols.append(np.arange(len(unknown))); data.append(diag)
        matrix = sparse.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                                   shape=(len(unknown),)*2)
        values[unknown] = spsolve(matrix, rhs)
        self.field_residual = float(np.max(np.abs(matrix @ values[unknown] - rhs)) / max(np.max(np.abs(rhs)), 1))
        return values

    def rate(self, barrier):
        # Symmetric barriers are guaranteed positive in the explored regime;
        # clipping caps downhill rates at the attempt frequency.
        return self.c.attempt_Hz * np.exp(-np.maximum(barrier, 0) / (KB*self.c.temperature_K))

    def catalogue(self):
        c, s = self.c, self.state
        groups = []
        def add(kind, a, b, barrier):
            a = np.asarray(a, int)
            if a.size:
                groups.append((kind, a, np.broadcast_to(b, a.shape), self.rate(barrier)))
        # Source electrode is an explicit open reservoir; return events debit it.
        a = np.arange(c.nx, 2*c.nx)
        dg = self.site[a] + self.phi[a] - c.voltage
        ok = s[a] == 0
        add('inject', a[ok], -1, c.source_barrier_eV + np.maximum(dg[ok], 0))
        ok = s[a] == 1
        add('absorb', a[ok], -1, c.source_barrier_eV + np.maximum(-dg[ok], 0))
        ions = np.flatnonzero(s == 1)
        a = np.repeat(ions, 4)
        b = self.neighbors[ions].ravel()
        ok = (b >= 0) & (s[np.maximum(b, 0)] == 0)
        a, b = a[ok], b[ok]
        dg = self.site[b]-self.site[a]+self.phi[b]-self.phi[a]
        base = np.where(self.layer[a] == self.layer[b],
                        np.where(self.layer[a] == 1, c.hop_alscn_eV, c.hop_aln_eV),
                        (c.hop_aln_eV+c.hop_alscn_eV)/2+c.interface_eV)
        add('hop', a, b, base + np.maximum(dg, 0))
        reachable = ((self.cathode[self.safe]) & self.valid).any(axis=1)
        if c.surface_redox_only:
            # Symmetric admissibility: neighbors do not change when a single
            # ion becomes metal, so both directions use the same host access.
            reachable &= (self.valid & (s[self.safe]!=2)).any(axis=1)
        # Electron reservoir at Pt=0; ion energy includes +q*phi.
        candidates = np.flatnonzero(self.interior & reachable & ((s == 1) | (s == 2)))
        dg = c.metal_standard_eV - self.binding(candidates) - self.site[candidates] - self.phi[candidates]
        pt_only = self.coord[candidates] == ((self.safe[candidates] >= self.n-c.nx) & self.valid[candidates]).sum(axis=1)
        base = c.reaction_barrier_eV + c.pt_nucleation_eV*pt_only
        red = s[candidates] == 1
        add('reduce', candidates[red], -1, base[red] + np.maximum(dg[red], 0))
        add('oxidize', candidates[~red], -1, base[~red] + np.maximum(-dg[~red], 0))
        # Bond-based neutral surface diffusion; shared saddle yields local DB.
        metals = np.flatnonzero(self.interior & (s == 2) & (self.coord < 4))
        a = np.repeat(metals, 4); b = self.neighbors[metals].ravel()
        ok = (b >= 0) & self.interior[np.maximum(b, 0)] & (s[np.maximum(b, 0)] == 0)
        a, b = a[ok], b[ok]
        dg = self.binding(a)-(self.binding(b)-c.bond_eV)
        add('migrate', a, b, c.migration_eV + np.maximum(dg, 0))
        return groups

    def binding(self, sites):
        """Effective Ag-Ag cohesion and Pt-contact adhesion, in eV."""
        neighbors=self.safe[sites]
        metal=(self.state[neighbors]==2)&self.valid[sites]
        pt=neighbors>=self.n-self.c.nx
        return np.sum(metal*np.where(pt,self.c.pt_adhesion_eV,self.c.bond_eV),axis=1)

    def save(self):
        self.snapshots.append(self.state.reshape(self.shape).copy())
        self.trace.append(dict(time_s=self.time, events=self.steps,
                               metal=int(((self.state == 2)&self.interior).sum()),
                               ions=int((self.state == 1).sum()), connected=self.connected,
                               net_source=self.net_source,
                               layer_metals=[int(np.sum((self.state==2)&self.interior&(self.zone==z))) for z in range(self.n_layers)],
                               layer_ions=[int(np.sum((self.state==1)&(self.zone==z))) for z in range(self.n_layers)],
                               hop_layer_counts=self.hop_layer_counts.tolist(),
                               net_deposition=self.net_deposition.tolist(),
                               ion_residence_s=self.ion_residence_s.tolist(),
                               vacancy_ion_residence_s=self.vacancy_ion_residence_s.tolist()))

    def advance(self, until_s, max_events=None):
        """Advance a constant-voltage segment, retaining a sampled future event.

        max_events is a cumulative ceiling, not an extra budget. Observational
        pauses preserve the exact CT-KMC path; voltage changes invalidate the
        pending event because its hazard was sampled under the previous field.
        First bridging still terminates V1; post-bridge physics is not enabled.
        """
        until_s = float(until_s)
        if not np.isfinite(until_s) or until_s < self.time:
            raise ValueError('until_s must be finite and not before current time')
        limit = self.c.max_events if max_events is None else int(max_events)
        if limit < self.steps: raise ValueError('event ceiling precedes current event count')
        if self.connected: return self.result('first_bridge')
        reason = 'protocol_complete_no_bridge'
        while self.time < until_s:
            if self.steps >= limit:
                reason = 'event_budget_unknown'; break
            if self.pending_event is None:
                proposal = self.sample_event()
                if proposal is None:
                    self.integrate_residence(until_s-self.time)
                    self.time = until_s; break
                self.pending_event = proposal
            proposal = self.pending_event
            if proposal['time_s'] > until_s:
                self.integrate_residence(until_s-self.time)
                self.time = until_s; break
            self.pending_event = None
            kind, a, b = proposal['kind'], proposal['a'], proposal['b']
            self.integrate_residence(proposal['time_s']-self.time)
            if kind=='hop': self.hop_layer_counts[self.zone[a],self.zone[b]]+=1
            elif kind=='reduce': self.net_deposition[self.zone[a]]+=1
            elif kind=='oxidize': self.net_deposition[self.zone[a]]-=1
            if kind == 'inject': self.state[a] = 1; self.net_source += 1
            elif kind == 'absorb': self.state[a] = 0; self.net_source -= 1
            elif kind in ('hop', 'migrate'): self.state[b] = self.state[a]; self.state[a] = 0
            elif kind == 'reduce': self.state[a] = 2
            elif kind == 'oxidize': self.state[a] = 1
            self.steps += 1; self.time = proposal['time_s']; self.counts[kind] += 1
            if self.c.record_events:
                self.event_log.append(dict(proposal, event_index=self.steps))
            if kind in ('reduce', 'oxidize', 'migrate'):
                self.update_structure()
                if self.connected:
                    reason = 'first_bridge'; break
                metal = int(((self.state == 2)&self.interior).sum())
                if metal >= self.snapshot_target:
                    self.save()
                    self.snapshot_target = next((x for x in [25,75,150,300,500,750,1100,1600]
                                                 if x>self.snapshot_target), self.n)
        if self.trace[-1]['events'] != self.steps or self.trace[-1]['time_s'] != self.time:
            self.save()
        return self.result(reason)

    def integrate_residence(self, dt):
        self.ion_residence_s += np.bincount(self.zone[self.state==1],minlength=self.n_layers)*dt
        if self._vacancy_indices.size:
            occupied = self._vacancy_indices[self.state[self._vacancy_indices] == 1]
            self.vacancy_ion_residence_s += np.bincount(self.zone[occupied], minlength=self.n_layers)*dt

    def sample_event(self):
        groups = self.catalogue()
        totals = np.array([g[3].sum() for g in groups])
        total = totals.sum()
        if total <= 0: return None
        dt = -np.log(max(self.rng.random(), np.finfo(float).tiny))/total
        draw = self.rng.random()*total
        k = min(int(np.searchsorted(np.cumsum(totals), draw, side='right')), len(groups)-1)
        kind, aa, bb, rates = groups[k]
        local = draw - totals[:k].sum()
        e = min(int(np.searchsorted(np.cumsum(rates), local, side='right')), len(aa)-1)
        return dict(kind=kind, a=int(aa[e]), b=int(bb[e]), time_s=float(self.time+dt))

    def result(self, reason):
        occupied = int(((self.state > 0)&self.interior).sum())
        return dict(config=asdict(self.c), status=reason, time_s=self.time, events=self.steps,
                    connected=self.connected, counts=self.counts,
                    ag_balance_error=occupied-self.net_source,
                    vacancy_affected_sites_by_layer=self.vacancy_affected_sites_by_layer.tolist(),
                    field_residual=self.field_residual, trace=self.trace)

    def run(self):
        return self.advance(self.c.max_time_s)

    def set_voltage(self, voltage):
        voltage = float(voltage)
        if not np.isfinite(voltage): raise ValueError('voltage must be finite')
        if voltage != self.c.voltage:
            if self.connected: raise RuntimeError('Post-bridge electronic physics is not implemented')
            self.c = replace(self.c, voltage=voltage)
            self.pending_event = None
            self.update_structure()

    def run_protocol(self, segments, max_events=None):
        normalized = [dict(voltage_V=float(s['voltage_V']), duration_s=float(s['duration_s'])) for s in segments]
        if any(not np.isfinite(s['duration_s']) or s['duration_s']<0 or not np.isfinite(s['voltage_V']) for s in normalized):
            raise ValueError('Invalid protocol segment')
        if self.protocol_position is None:
            self.protocol_position = dict(segments=normalized, index=0, deadline_s=None)
        position = self.protocol_position
        if position['segments'] != normalized: raise ValueError('Cannot replace an in-progress protocol')
        while position['index'] < len(normalized):
            segment = normalized[position['index']]
            if position['deadline_s'] is None:
                self.set_voltage(segment['voltage_V'])
                position['deadline_s'] = self.time+segment['duration_s']
            result = self.advance(position['deadline_s'], max_events)
            if result['status'] != 'protocol_complete_no_bridge': return result
            position['index'] += 1
            position['deadline_s'] = None
        return self.result('protocol_complete_no_bridge')

    def save_checkpoint(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        metadata = dict(schema=CHECKPOINT_SCHEMA, config=asdict(self.c), time_s=self.time, events=self.steps,
                        net_source=self.net_source, counts=self.counts, rng=self.rng.bit_generator.state,
                        pending_event=self.pending_event, snapshot_target=self.snapshot_target,
                        trace=self.trace, event_log=self.event_log, protocol_position=self.protocol_position,
                        connected=self.connected, field_residual=self.field_residual,
                        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        arrays = dict(state=self.state, site=self.site, phi=self.phi, snapshots=np.array(self.snapshots),
                      hop_layer_counts=self.hop_layer_counts,net_deposition=self.net_deposition,
                      ion_residence_s=self.ion_residence_s, vacancy_mask=self.vacancy_mask,
                      vacancy_ion_residence_s=self.vacancy_ion_residence_s)
        metadata_text = json.dumps(metadata, sort_keys=True)
        checksum = hashlib.sha256(metadata_text.encode()+b''.join(v.tobytes() for v in arrays.values())).hexdigest()
        with path.open('xb') as f:
            np.savez_compressed(f, **arrays, metadata_json=metadata_text, checksum=checksum)
        return path

    @classmethod
    def load_checkpoint(cls, path):
        with np.load(path, allow_pickle=False) as data:
            text = str(data['metadata_json'])
            metadata = json.loads(text)
            if metadata['schema'] != CHECKPOINT_SCHEMA:
                raise ValueError('Unsupported checkpoint schema; use the matching archived source')
            arrays = {k:data[k].copy() for k in ['state','site','phi','snapshots','hop_layer_counts',
                      'net_deposition','ion_residence_s','vacancy_mask','vacancy_ion_residence_s']}
            checksum = hashlib.sha256(text.encode()+b''.join(v.tobytes() for v in arrays.values())).hexdigest()
            if checksum != str(data['checksum']): raise ValueError('Checkpoint checksum mismatch')
        if metadata['source_sha256'] != hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
            raise ValueError('Checkpoint source version mismatch')
        model = cls(Config(**metadata['config']))
        if any(arrays[k].shape!=(model.n,) for k in ['state','site','phi']): raise ValueError('Checkpoint shape mismatch')
        if not np.array_equal(arrays['vacancy_mask'], model.vacancy_mask):
            raise ValueError('Checkpoint vacancy mask differs from the configured static field')
        vacancy_residence = arrays['vacancy_ion_residence_s']
        expected_shapes = {'hop_layer_counts': (model.n_layers, model.n_layers),
                           'net_deposition': (model.n_layers,), 'ion_residence_s': (model.n_layers,)}
        if any(arrays[key].shape != shape for key, shape in expected_shapes.items()):
            raise ValueError('Checkpoint layer statistics shape mismatch')
        if (vacancy_residence.shape != (model.n_layers,) or not np.isfinite(vacancy_residence).all() or
                (vacancy_residence < 0).any() or (vacancy_residence > arrays['ion_residence_s'] + 1e-12).any()):
            raise ValueError('Invalid checkpoint vacancy residence totals')
        model.state, model.site, model.phi = arrays['state'], arrays['site'], arrays['phi']
        if not np.isin(model.state,[0,1,2]).all() or not (model.state[~model.interior]==2).all():
            raise ValueError('Invalid state/electrodes')
        if int(np.count_nonzero(model.state[model.interior])) != metadata['net_source']:
            raise ValueError('Checkpoint Ag ledger mismatch')
        model.time=metadata['time_s']; model.steps=metadata['events']; model.net_source=metadata['net_source']
        for key in ['counts','pending_event','snapshot_target','trace','event_log','protocol_position','field_residual']:
            setattr(model,key,metadata[key])
        model.snapshots=list(arrays['snapshots'])
        for key in ['hop_layer_counts','net_deposition','ion_residence_s']: setattr(model,key,arrays[key])
        model.vacancy_ion_residence_s = vacancy_residence
        model.rng.bit_generator.state=metadata['rng']
        model.cathode=model.component(-1); model.anode=model.component(0)
        model.connected=bool(model.cathode[:model.c.nx].any())
        if model.connected!=metadata['connected']: raise ValueError('Checkpoint topology mismatch')
        model.coord=((model.state[model.safe]==2)&model.valid).sum(axis=1)
        return model
