"""Morphology-first, coarse-grained Ag/AlN/AlScN/AlN/Pt CT-KMC.

All energies eV, potentials V, distances nm. One occupied voxel is an
effective Ag packet, NOT one resolved atom. No prescribed growth path.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy import sparse, ndimage
from scipy.sparse.linalg import spsolve

KB = 8.617333262e-5
EPS0 = 8.8541878128e-12
CROSS = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])


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
    migration_eV: float = 0.80
    disorder_eV: float = 0.025
    max_events: int = 350000
    max_time_s: float = 2.0
    seed: int = 7


class Model:
    def __init__(self, cfg):
        self.c = cfg
        self.rng = np.random.default_rng(cfg.seed)
        self.nz = 3 * cfg.layer_cells + 2
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
        self.layer.reshape(self.shape)[1+cfg.layer_cells:1+2*cfg.layer_cells] = 1
        self.site = self.rng.normal(0, cfg.disorder_eV, self.n)
        self.time = 0.0
        self.steps = 0
        self.net_source = 0
        self.counts = {k: 0 for k in ('inject', 'absorb', 'hop', 'reduce', 'oxidize', 'migrate')}
        self.snapshots = []
        self.trace = []
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
            rows.extend(r[~f]); cols.extend(index[j[~f]]); data.extend(-g[~f])
        rows.extend(np.arange(len(unknown))); cols.extend(np.arange(len(unknown))); data.extend(diag)
        matrix = sparse.csr_matrix((data, (rows, cols)), shape=(len(unknown),)*2)
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
        # Electron reservoir at Pt=0; ion energy includes +q*phi.
        candidates = np.flatnonzero(self.interior & reachable & ((s == 1) | (s == 2)))
        dg = c.metal_standard_eV - c.bond_eV*self.coord[candidates] - self.site[candidates] - self.phi[candidates]
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
        dg = -c.bond_eV*(self.coord[b]-1-self.coord[a])
        add('migrate', a, b, c.migration_eV + np.maximum(dg, 0))
        return groups

    def save(self):
        self.snapshots.append(self.state.reshape(self.shape).copy())
        self.trace.append(dict(time_s=self.time, events=self.steps,
                               metal=int(((self.state == 2)&self.interior).sum()),
                               ions=int((self.state == 1).sum()), connected=self.connected,
                               net_source=self.net_source))

    def run(self):
        targets = iter([25, 75, 150, 300, 500, 750, 1100, 1600])
        target = next(targets)
        reason = 'event_budget'
        while self.steps < self.c.max_events:
            groups = self.catalogue()
            totals = np.array([g[3].sum() for g in groups])
            total = totals.sum()
            if total <= 0:
                reason = 'no_events'; break
            dt = -np.log(max(self.rng.random(), np.finfo(float).tiny))/total
            if self.time+dt > self.c.max_time_s:
                self.time = self.c.max_time_s; reason = 'time_limit'; break
            # Independent waiting-time and event random draws.
            draw = self.rng.random()*total
            k = min(int(np.searchsorted(np.cumsum(totals), draw, side='right')), len(groups)-1)
            kind, aa, bb, rates = groups[k]
            local = draw - totals[:k].sum()
            e = min(int(np.searchsorted(np.cumsum(rates), local, side='right')), len(aa)-1)
            a, b = aa[e], bb[e]
            if kind == 'inject': self.state[a] = 1; self.net_source += 1
            elif kind == 'absorb': self.state[a] = 0; self.net_source -= 1
            elif kind in ('hop', 'migrate'): self.state[b] = self.state[a]; self.state[a] = 0
            elif kind == 'reduce': self.state[a] = 2
            elif kind == 'oxidize': self.state[a] = 1
            self.steps += 1; self.time += dt; self.counts[kind] += 1
            if self.steps % 50000 == 0:
                print(f'seed={self.c.seed} P={self.c.polarization_C_m2:+.4f} step={self.steps} t={self.time:.4g}s Ag={np.sum((self.state==2)&self.interior)}', flush=True)
            if kind in ('reduce', 'oxidize', 'migrate'):
                self.update_structure()
                if self.connected:
                    reason = 'first_bridge'; break
                metal = int(((self.state == 2)&self.interior).sum())
                if metal >= target:
                    self.save(); target = next(targets, self.n)
        self.save()
        occupied = int(((self.state > 0)&self.interior).sum())
        return dict(config=asdict(self.c), status=reason, time_s=self.time, events=self.steps,
                    connected=self.connected, counts=self.counts,
                    ag_balance_error=occupied-self.net_source,
                    field_residual=self.field_residual, trace=self.trace)
