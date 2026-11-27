"""Continuous electronic-state participation across SET/RESET/re-SET.

D and H are reversible *damage coordinates* in effective voltage-sample
units.  They are not cumulative unsigned dose counters: negative device bias
adds damage, whereas positive carrier flow removes it continuously.  The
participating electronic fraction is always derived from D and H, never set
according to a stage label.  VN positions and Ag geometry are separate state.

The SET recovery law is an effective, charge-gated hypothesis to calibrate:
    k_rec = recovery_rate * I_material/I_ref * V/(V+V_rec), V > 0.
Thus small positive device voltage during current compliance can recover the
electronic state.  A disconnected filament only supplies its actual small
electronic current; there is no abrupt metallic-connectivity gate.

Pass the actual voltage across the material, and its positive electronic
current, excluding the independently modelled background leakage.  No field
solver, measurement averaging or clock-to-seconds conversion occurs here.
"""
from dataclasses import asdict, dataclass, fields
from pathlib import Path
import json
import math


@dataclass(frozen=True)
class ActivityParams:
    dose_scale: float = 4.409115796801806
    activity_exponent: float = .7778720480535988
    high_field_rate: float = .4240697032605897
    field_scale_V: float = .12
    high_field_reference_V: float = 1.5
    recovery_rate: float = .3
    recovery_voltage_scale_V: float = .3
    reference_current_A: float = .001

    def __post_init__(self):
        if not all(math.isfinite(float(v)) for v in asdict(self).values()):
            raise ValueError('Activity parameters must be finite.')
        for name in ('dose_scale', 'field_scale_V', 'recovery_voltage_scale_V',
                     'reference_current_A'):
            if getattr(self, name) <= 0:
                raise ValueError(f'{name} must be positive.')
        for name in ('activity_exponent', 'high_field_rate',
                     'high_field_reference_V', 'recovery_rate'):
            if getattr(self, name) < 0:
                raise ValueError(f'{name} must be nonnegative.')

    @classmethod
    def from_reset_summary(cls, path, **overrides):
        data = json.loads(Path(path).read_text(encoding='utf-8'))['parameters']
        allowed = {f.name for f in fields(cls)}
        selected = {k: v for k, v in data.items() if k in allowed}
        selected.update(overrides)
        return cls(**selected)


@dataclass
class ActivityState:
    D: float = 0.
    H: float = 0.
    params: ActivityParams = ActivityParams()
    elapsed_units: float = 0.
    recovery_exposure: float = 0.

    def __post_init__(self):
        for name in ('D', 'H', 'elapsed_units', 'recovery_exposure'):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'{name} must be finite and nonnegative.')
            setattr(self, name, value)
        if not isinstance(self.params, ActivityParams):
            raise TypeError('params must be ActivityParams.')

    @property
    def q(self):
        p = self.params
        return math.exp(-p.activity_exponent * math.log1p(self.D / p.dose_scale)
                        -p.high_field_rate * self.H)

    @property
    def activity(self):
        return self.q

    def recovery_hazard(self, voltage_V, current_A):
        """Recovery rate at constant positive material voltage and current."""
        v, current = float(voltage_V), float(current_A)
        if not math.isfinite(v) or not math.isfinite(current):
            raise ValueError('Voltage and material current must be finite.')
        if v <= 0:
            return 0.
        if current < 0:
            raise ValueError('Positive bias requires a nonnegative material current.')
        p = self.params
        return (p.recovery_rate * current / p.reference_current_A *
                v / (v + p.recovery_voltage_scale_V))

    def advance(self, voltage_V, dt=1., current_A=None):
        """Advance exactly for constant V and I over dt; mutate state, return q.

        A coupled circuit caller should subdivide a step when q changes the
        circuit current or compliance voltage.  This exact constant-input
        update does not itself assume that current stays fixed for a full SMU
        measurement sample.  Positive updates require the circuit current.
        """
        v, dt = float(voltage_V), float(dt)
        if not math.isfinite(v) or not math.isfinite(dt) or dt < 0:
            raise ValueError('Finite voltage and nonnegative finite dt required.')
        if dt == 0:
            return self.q
        if v < 0:
            u = -v
            p = self.params
            delta_D = u * dt
            delta_H = delta_D * math.exp((u-p.high_field_reference_V)/p.field_scale_V)
            new_D, new_H = self.D + delta_D, self.H + delta_H
            if not math.isfinite(new_D) or not math.isfinite(new_H):
                raise ValueError('Exposure overflow; check voltage and time units.')
            self.D, self.H = new_D, new_H
        elif v > 0:
            if current_A is None:
                raise ValueError('Supply actual material current for positive-bias recovery.')
            hazard = self.recovery_hazard(v, current_A) * dt
            if not math.isfinite(hazard):
                raise ValueError('Recovery exposure overflow.')
            decay = math.exp(-hazard)
            self.D *= decay
            self.H *= decay
            self.recovery_exposure += hazard
        self.elapsed_units += dt
        return self.q

    def copy(self):
        return type(self).from_dict(self.to_dict())

    def to_dict(self):
        return dict(schema='cycle_activity40/v1', D=self.D, H=self.H, q=self.q,
                    elapsed_units=self.elapsed_units,
                    recovery_exposure=self.recovery_exposure,
                    params=asdict(self.params),
                    units='effective sampling-time units; not calibrated seconds')

    @classmethod
    def from_dict(cls, data):
        if data.get('schema', 'cycle_activity40/v1') != 'cycle_activity40/v1':
            raise ValueError('Unsupported electronic-state schema.')
        obj = cls(D=data['D'], H=data['H'], params=ActivityParams(**data['params']),
                  elapsed_units=data.get('elapsed_units', 0.),
                  recovery_exposure=data.get('recovery_exposure', 0.))
        if 'q' in data and not math.isclose(obj.q, float(data['q']), rel_tol=1e-11, abs_tol=1e-15):
            raise ValueError('Stored activity is inconsistent with D/H and parameters.')
        return obj

    def save(self, path):
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding='utf-8')

    @classmethod
    def restore(cls, path):
        return cls.from_dict(json.loads(Path(path).read_text(encoding='utf-8')))

    @classmethod
    def from_retention(cls, path, params=None):
        """Import the fitted LRS/HRS without resetting q or either coordinate."""
        path = Path(path)
        if path.is_dir():
            path = path / 'read_state.json'
        data = json.loads(path.read_text(encoding='utf-8'))
        e = data['electronic_state']
        if params is None:
            allowed = {f.name for f in fields(ActivityParams)}
            params = ActivityParams(**{k: v for k, v in e.items() if k in allowed})
        obj = cls(D=e['low_field_dose_V_samples'], H=e['high_field_dose_V_samples'],
                  params=params)
        if not math.isclose(obj.q, float(e['activity']), rel_tol=1e-11, abs_tol=1e-15):
            raise ValueError('Retention q does not match the supplied damage-law parameters.')
        return obj


def replay(voltage, initial=None, dt=1., material_current=None):
    """Convenience replay for specified inputs; no self-consistent circuit solve."""
    state = ActivityState() if initial is None else initial.copy()
    values = []
    for i, v in enumerate(voltage):
        current = None if material_current is None else material_current[i]
        values.append(state.advance(float(v), dt, current))
    return values, state
