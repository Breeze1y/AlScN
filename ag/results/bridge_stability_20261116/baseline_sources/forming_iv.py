"""Exploratory single-forming I-V: KMC + electronic network + ideal compliance.

The original prebridge KMC remains unchanged. Electrostatics still neglects
ionic space charge and pins electrode-connected metal before bridging.
After first bridge morphology is frozen; only circuit readout continues.
"""
from dataclasses import dataclass, asdict, replace
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve
from kmc_model import Model, Config


@dataclass(frozen=True)
class ElectronicConfig:
    metal_resistivity_ohm_m: float = 1.6e-8
    effective_depth_nm: float = 1.0
    contact_resistance_ohm: float = 500.0
    tunnel_prefactor_S: float = 7.748091729e-5
    tunnel_barrier_eV: float = 1.0
    tunnel_effective_mass_ratio: float = 1.0
    max_tunnel_gap_nm: float = 2.0
    min_tunnel_conductance_S: float = 1e-15
    leakage_per_column_S: float = 1e-12
    leakage_voltage_scale_V: float = 2.0
    series_resistance_ohm: float = 100.0
    compliance_A: float = 1e-3
    leakage_model: str = 'sinh'
    leakage_saturation_A: float = 1e-6
    leakage_half_voltage_V: float = 1.0
    leakage_exponent: float = 2.0
    leakage_reference_columns: float = 150.0
    network_parallel_factor: float = 1.0

    def __post_init__(self):
        values = {k:v for k,v in asdict(self).items() if k!='leakage_model'}
        if self.leakage_model not in ['sinh','saturating_power','ohmic_saturating_power']:
            raise ValueError('Unknown empirical leakage model')
        if any(not np.isfinite(v) for v in values.values()):
            raise ValueError('Electronic parameters must be finite')
        positive = ['metal_resistivity_ohm_m', 'effective_depth_nm',
                    'tunnel_prefactor_S', 'tunnel_barrier_eV',
                    'tunnel_effective_mass_ratio', 'leakage_voltage_scale_V', 'compliance_A',
                    'leakage_saturation_A','leakage_half_voltage_V','leakage_exponent',
                    'leakage_reference_columns','network_parallel_factor']
        if any(values[k] <= 0 for k in positive):
            raise ValueError('Resistivity, depth, tunnel parameters, voltage scale and compliance must be positive')
        if any(v < 0 for v in values.values()):
            raise ValueError('Electronic parameters must be nonnegative')


def electronic_network(state, spacing_nm, c):
    """Unit-bias resistor graph. Ions are not metal conductors.

    Adjacent Ag sites have rho*l/A resistance, with extra electrode contact R.
    Axial metal-to-metal gaps have G0 exp(-2 sqrt(2 m Phi) d / hbar).
    d is the full count of intervening dielectric cells times grid spacing.
    Edges skip neither metal nor intervening conducting nodes. The effective
    cross-section depth and tunnel channel count are uncalibrated parameters.
    """
    s = np.asarray(state)
    if s.ndim != 2 or s.shape[0] < 2 or s.shape[1] < 1 or not np.isfinite(spacing_nm) or spacing_nm <= 0 or not np.isin(s, [0,1,2]).all():
        raise ValueError('Expected a 2D occupancy lattice and positive spacing')
    if not np.all(s[[0,-1]] == 2):
        raise ValueError('Electrode rows must be metal')
    metal = s == 2
    nz, nx = s.shape
    ids = np.full(s.shape, -1, int)
    ids[0] = 0; ids[-1] = 1
    inner = np.flatnonzero(metal.ravel() & (np.indices(s.shape)[0].ravel()>0)
                           & (np.indices(s.shape)[0].ravel()<nz-1))
    ids.ravel()[inner] = np.arange(2, len(inner)+2)
    n = len(inner)+2
    rr, cc, gg, tunnel = [], [], [], []
    beta_per_nm = 10.2463344*np.sqrt(c.tunnel_effective_mass_ratio*c.tunnel_barrier_eV)
    metal_R = c.metal_resistivity_ohm_m/(c.effective_depth_nm*1e-9)
    max_gap_cells = int(np.floor(c.max_tunnel_gap_nm/spacing_nm + 1e-12))
    for axis in [0,1]:
        for distance in range(1, max_gap_cells+2):
            if distance >= s.shape[axis]: break
            start = [slice(None),slice(None)]; end = start.copy()
            start[axis] = slice(None,-distance); end[axis] = slice(distance,None)
            start, end = tuple(start), tuple(end)
            valid = metal[start] & metal[end]
            for offset in range(1,distance):
                mid = [slice(None),slice(None)]
                mid[axis] = slice(offset, s.shape[axis]-distance+offset)
                valid &= ~metal[tuple(mid)]
            a, b = ids[start][valid], ids[end][valid]
            keep = a != b; a, b = a[keep], b[keep]
            if not len(a): continue
            if distance == 1:
                contacts = (a<2).astype(int)+(b<2).astype(int)
                g = 1/(metal_R + contacts*c.contact_resistance_ohm)
            else:
                g0 = c.tunnel_prefactor_S*np.exp(-beta_per_nm*(distance-1)*spacing_nm)
                if g0 < c.min_tunnel_conductance_S: continue
                g = np.full(len(a),g0)
            rr.extend(a.tolist()); cc.extend(b.tolist()); gg.extend(g.tolist())
            tunnel.extend([distance>1]*len(a))
    a, b, g = np.asarray(rr,int), np.asarray(cc,int), np.asarray(gg,float)*c.network_parallel_factor
    tunnel = np.asarray(tunnel,bool)
    if not len(g):
        return dict(conductance_S=0., kcl_residual_A_at_1V=0., tunnel_edges=0,
                    metal_edges=0, node_count=n)
    adjacency = sparse.coo_matrix((np.ones(2*len(a)),
                                  (np.r_[a,b],np.r_[b,a])),shape=(n,n)).tocsr()
    _, labels = connected_components(adjacency, directed=False)
    if labels[0] != labels[1]:
        return dict(conductance_S=0., kcl_residual_A_at_1V=0.,
                    tunnel_edges=int(tunnel.sum()),
                    metal_edges=int((~tunnel).sum()), node_count=n)
    keep = labels[a] == labels[0]
    a,b,g = a[keep],b[keep],g[keep]
    tunnel = tunnel[keep]
    active = np.flatnonzero(labels == labels[0])
    unknown = active[active >= 2]
    phi = np.zeros(n); phi[0] = 1
    L = sparse.coo_matrix((np.r_[g,g,-g,-g],
                         (np.r_[a,b,a,b],np.r_[a,b,b,a])),shape=(n,n)).tocsr()
    if len(unknown):
        phi[unknown] = spsolve(L[unknown][:,unknown], -L[unknown,0].toarray().ravel())
    # Dissipated power at one volt avoids terminal-current cancellation in HRS.
    conductance = float(np.sum(g*(phi[a]-phi[b])**2))
    residual = float(np.max(np.abs((L@phi)[unknown]))) if len(unknown) else 0.
    if not np.isfinite(conductance) or not np.isfinite(residual):
        raise RuntimeError('Nonfinite electronic network solution')
    if residual > max(1e-12, conductance*1e-6):
        raise RuntimeError('Electronic network KCL residual exceeded tolerance')
    return dict(conductance_S=conductance, kcl_residual_A_at_1V=residual,
                tunnel_edges=int(tunnel.sum()),
                metal_edges=int((~tunnel).sum()), node_count=n)


def leakage_current(voltage_V, nx, c):
    """Empirical background only; neither option identifies a transport mechanism."""
    if c.leakage_model=='sinh':
        return float(nx*c.leakage_per_column_S*c.leakage_voltage_scale_V*np.sinh(voltage_V/c.leakage_voltage_scale_V))
    magnitude=abs(voltage_V)
    if magnitude==0:return 0.
    # Stable Hill saturation, parameterized at the entire reference window.
    x=c.leakage_exponent*np.log(magnitude/c.leakage_half_voltage_V)
    fraction=1/(1+np.exp(-np.clip(x,-700,700)))
    saturation = np.sign(voltage_V)*c.leakage_saturation_A*nx/c.leakage_reference_columns*fraction
    ohmic = nx*c.leakage_per_column_S*voltage_V if c.leakage_model=='ohmic_saturating_power' else 0.
    return float(saturation + ohmic)


def circuit_point(command_V, conductance_S, nx, c):
    """Solve monotonic load line, then reduce voltage if compliance is active.

    Leakage is an odd empirical background: N*g_leak*V0*sinh(V/V0).
    Source command and actual device voltage are distinct saved quantities.
    """
    if not np.isfinite(command_V) or not np.isfinite(conductance_S) or conductance_S < 0 or nx<=0:
        raise ValueError('Invalid circuit input')
    def parts(v):
        leak = leakage_current(v,nx,c)
        return conductance_S*v, float(leak)
    def current(v): return sum(parts(v))
    sign = -1 if command_V < 0 else 1
    command = abs(command_V)
    if command/c.leakage_voltage_scale_V > 500:
        raise ValueError('Voltage outside leakage model range')
    lo, hi = 0., command
    for _ in range(65):
        mid = (lo+hi)/2
        if mid + current(mid)*c.series_resistance_ohm > command: hi=mid
        else: lo=mid
    v = (lo+hi)/2
    limited = current(v) > c.compliance_A
    if limited:
        lo, hi = 0., v
        for _ in range(65):
            mid=(lo+hi)/2
            if current(mid) > c.compliance_A: hi=mid
            else: lo=mid
        v=(lo+hi)/2
    metal, leak = parts(v)
    return dict(command_voltage_V=float(command_V), device_voltage_V=sign*v,
                current_A=sign*(metal+leak), network_current_A=sign*metal,
                leakage_current_A=sign*leak, compliance_active=bool(limited),
                source_delivered_voltage_V=sign*(v+(metal+leak)*c.series_resistance_ohm))


class FormingModel(Model):
    """Prebridge Model with event-wise electronic readout and voltage feedback."""
    def __init__(self, cfg, electronics):
        self.electronics = electronics
        self.command_voltage = cfg.voltage
        self.iv = []
        self.electronic_updates = 0
        self.first_bridge = None
        super().__init__(cfg)

    def update_structure(self):
        self.cathode = self.component(-1)
        self.anode = self.component(0)
        self.connected = bool(self.cathode[:self.c.nx].any())
        self.coord = ((self.state[self.safe] == 2)&self.valid).sum(axis=1)
        self.network = electronic_network(self.state.reshape(self.shape), self.c.spacing_nm, self.electronics)
        self.electronic_updates += 1
        self._feedback()
        if not self.connected:
            self.phi = self.solve_field()
        elif self.first_bridge is None:
            self.first_bridge = dict(time_s=self.time, events=self.steps,
                                     command_voltage_V=self.command_voltage,
                                     device_voltage_V=self.operating['device_voltage_V'])
        self.record('metal_change' if self.steps else 'initial')

    def _feedback(self):
        self.operating = circuit_point(self.command_voltage, self.network['conductance_S'],
                                       self.c.nx, self.electronics)
        voltage = self.operating['device_voltage_V']
        if voltage != self.c.voltage:
            self.c = replace(self.c, voltage=voltage)
            self.pending_event = None

    def apply_command(self, voltage):
        if not np.isfinite(voltage): raise ValueError('Invalid source voltage')
        self.command_voltage = float(voltage)
        before = self.c.voltage
        self._feedback()
        if not self.connected and self.c.voltage != before:
            self.phi = self.solve_field()

    def save_checkpoint(self, path):
        raise NotImplementedError('Forming circuit checkpoints require a dedicated schema; use archived I-V outputs')

    @classmethod
    def load_checkpoint(cls, path):
        raise NotImplementedError('Morphology-only checkpoints cannot restore a forming circuit')

    def record(self, kind, read_time_s=None):
        self.iv.append(dict(sample=len(self.iv), time_s=self.time if read_time_s is None else read_time_s,
                            morphology_time_s=self.time, events=self.steps, kind=kind,
                            connected=self.connected, metal=int(np.sum(self.state[self.interior]==2)),
                            ions=int(np.sum(self.state[self.interior]==1)),
                            net_source=self.net_source,
                            ag_balance_error=int(np.count_nonzero(self.state[self.interior]))-self.net_source,
                            network_conductance_S=self.network['conductance_S'],
                            network_kcl_residual_A_at_1V=self.network['kcl_residual_A_at_1V'],
                            **self.operating))


def run_forming(cfg, electronics, protocol, on_segment=None, allow_frozen_return=False, model_class=None):
    """Piecewise constant ramp; no invented points after event-budget truncation."""
    if not protocol: raise ValueError('Empty forming protocol')
    segments = [dict(voltage_V=float(s['voltage_V']), duration_s=float(s['duration_s'])) for s in protocol]
    if any(not np.isfinite(s['duration_s']) or s['duration_s']<=0 or
           not np.isfinite(s['voltage_V']) or s['voltage_V']<0 for s in segments):
        raise ValueError('Single-forming protocol requires finite nonnegative voltages and positive durations')
    descending=any(b['voltage_V'] < a['voltage_V'] for a,b in zip(segments,segments[1:]))
    if descending and not allow_frozen_return:
        raise ValueError('Single-forming scan must be nondecreasing; RESET is not implemented')
    if descending:
        peak=int(np.argmax([s['voltage_V'] for s in segments]))
        if any(b['voltage_V']<a['voltage_V'] for a,b in zip(segments[:peak],segments[1:peak+1])) or any(
            b['voltage_V']>a['voltage_V'] for a,b in zip(segments[peak:],segments[peak+1:])):
            raise ValueError('Frozen return permits one positive triangular sweep only')
    if sum(s['duration_s'] for s in segments) > cfg.max_time_s+1e-12:
        raise ValueError('Protocol duration exceeds configured max_time_s')
    cls=FormingModel if model_class is None else model_class
    m = cls(replace(cfg,voltage=segments[0]['voltage_V']),electronics)
    scheduled = 0.
    scan_states, scan_records = [], []
    status = 'protocol_complete_no_bridge'
    for i, segment in enumerate(segments):
        if i and segment['voltage_V'] < segments[i-1]['voltage_V'] and not m.connected:
            status='return_not_read_no_bridge';break
        m.apply_command(segment['voltage_V'])
        m.record('frozen_bridge_readout' if m.connected else 'voltage_step', scheduled)
        deadline = scheduled + segment['duration_s']
        if not m.connected:
            result = m.advance(deadline)
            status = result['status']
            m.record('segment_end')
            if m.connected and m.time < deadline:
                m.record('frozen_bridge_readout', deadline)
        else:
            m.record('frozen_bridge_readout', deadline)
        scan_states.append(m.state.reshape(m.shape).copy())
        scan_records.append(dict(segment=i, scheduled_end_s=deadline, **m.iv[-1]))
        if on_segment: on_segment(i, m, status)
        if status == 'event_budget_unknown': break
        scheduled = deadline
    if m.connected: status = 'first_bridge_frozen_readout'
    summary = dict(status=status, first_bridge=m.first_bridge,
                   events=m.steps, morphology_time_s=m.time, scan_end_time_s=m.iv[-1]['time_s'],
                   ag_balance_error=m.iv[-1]['ag_balance_error'],
                   electronic_updates=m.electronic_updates, config=asdict(cfg),
                   electronics=asdict(electronics), protocol=segments,
                   post_bridge_mode='frozen_morphology_readout',
                   allow_frozen_return=bool(allow_frozen_return),
                   calibration_status='exploratory_not_fitted_to_experiment',
                   field_scope='prebridge_ideal_electrode_connected_metal_no_ionic_space_charge',
                   vacancy_affected_sites_by_layer=m.vacancy_affected_sites_by_layer.tolist(),
                   ion_residence_s=m.ion_residence_s.tolist(),
                   vacancy_ion_residence_s=m.vacancy_ion_residence_s.tolist(),
                   counts=m.counts)
    return m, summary, np.asarray(scan_states), scan_records
