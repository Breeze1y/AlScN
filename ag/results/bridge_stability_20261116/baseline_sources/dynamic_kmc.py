"""Continuous forming return and zero-bias diagnostics on the infinite Ag model."""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve
from kmc_model import Model, Config, EPS0
from forming_iv import FormingModel, ElectronicConfig, electronic_network, circuit_point
from fast_kmc import CatalogueCache

def dynamic_network(state,spacing_nm,c):
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
    phi=np.full(n,np.nan); phi[0]=1.0;phi[1]=0.0
    if len(g):
        adjacency=sparse.coo_matrix((np.ones(2*len(a)),(np.r_[a,b],np.r_[b,a])),shape=(n,n)).tocsr()
        _,labels=connected_components(adjacency,directed=False)
        source=labels==labels[0];pt=labels==labels[1]
        phi[source]=1.0;phi[pt]=0.0;phi[0]=1.0;phi[1]=0.0
        conductance=0.0;residual=0.0
        if labels[0]==labels[1]:
            unknown=np.flatnonzero(source);unknown=unknown[unknown>=2]
            L=sparse.coo_matrix((np.r_[g,g,-g,-g],(np.r_[a,b,a,b],np.r_[a,b,b,a])),shape=(n,n)).tocsr()
            if len(unknown):phi[unknown]=spsolve(L[unknown][:,unknown],-L[unknown,0].toarray().ravel())
            active=source[a]&source[b]
            conductance=float(np.sum(g[active]*(phi[a[active]]-phi[b[active]])**2))
            residual=float(np.max(np.abs((L@np.nan_to_num(phi))[unknown]),initial=0))
    else:
        conductance=residual=0.0
    if not np.isfinite(conductance) or residual>max(1e-12,conductance*1e-6):
        raise RuntimeError('Dynamic network KCL residual exceeded tolerance')
    potential=np.full(s.shape,np.nan);potential[metal]=phi[ids[metal]]
    return dict(conductance_S=conductance,kcl_residual_A_at_1V=residual,
                tunnel_edges=int(tunnel.sum()),metal_edges=int((~tunnel).sum()),node_count=n,
                metal_phi_at_1V=potential)

class DynamicModel(FormingModel):
    def __init__(self,cfg,electronics):
        self.topology_log=[];self.first_bridge_state=None
        super().__init__(cfg,electronics)

    def update_structure(self):
        previous=getattr(self,'connected',False)
        self.cathode=self.component(-1);self.anode=self.component(0)
        self.connected=bool(self.cathode[:self.c.nx].any())
        self.coord=((self.state[self.safe]==2)&self.valid).sum(axis=1)
        historical=self.first_bridge is None and not self.connected
        self.network=(electronic_network if historical else dynamic_network)(self.state.reshape(self.shape),self.c.spacing_nm,self.electronics)
        self.electronic_updates+=1;self._feedback()
        self.phi=self.solve_field()
        if previous!=self.connected:
            entry=dict(time_s=self.time,events=self.steps,connected=self.connected,**self.operating)
            self.topology_log.append(entry)
            if self.connected and self.first_bridge is None:
                self.first_bridge={k:entry[k] for k in ('time_s','events','command_voltage_V','device_voltage_V')}
                self.first_bridge_state=self.state.reshape(self.shape).copy()
        self.record('metal_change' if self.steps else 'initial')

    def apply_command(self,voltage):
        if not np.isfinite(voltage):raise ValueError('Invalid command')
        self.command_voltage=float(voltage);before=self.c.voltage;self._feedback()
        if self.c.voltage!=before:self.phi=self.solve_field()

    def solve_field(self):
        if self.first_bridge is None and not self.connected:
            return Model.solve_field(self)
        c=self.c;h=c.spacing_nm*1e-9
        metal=self.state==2
        eps=EPS0*np.where(self.layer==1,c.eps_alscn,c.eps_aln);eps[metal]=EPS0*1e5
        polarization=np.where(self.layer==1,c.polarization_C_m2,0.);polarization[metal]=0
        boundary=self.network['metal_phi_at_1V'].ravel()*self.operating['device_voltage_V']
        fixed=np.isfinite(boundary);values=np.where(fixed,boundary,0.)
        unknown=np.flatnonzero(~fixed)
        if not len(unknown):self.field_residual=0.;return values
        index=np.full(self.n,-1,int);index[unknown]=np.arange(len(unknown))
        rr=[];cc=[];data=[];rhs=np.zeros(len(unknown));diagonal=np.zeros(len(unknown))
        for direction in range(4):
            i=unknown[self.valid[unknown,direction]];j=self.neighbors[i,direction];r=index[i]
            w=2*eps[i]*eps[j]/(eps[i]+eps[j])/h**2;diagonal[r]+=w
            f=fixed[j];rhs[r[f]]+=w[f]*values[j[f]]
            if direction<2:
                pf=(eps[j]*polarization[i]+eps[i]*polarization[j])/(eps[i]+eps[j])
                rhs[r]-=(1 if direction==1 else -1)*pf/h
            rr.extend(r[~f]);cc.extend(index[j[~f]]);data.extend(-w[~f])
        rr.extend(range(len(unknown)));cc.extend(range(len(unknown)));data.extend(diagonal)
        A=sparse.csr_matrix((data,(rr,cc)),shape=(len(unknown),)*2)
        values[unknown]=spsolve(A,rhs)
        self.field_residual=float(np.max(abs(A@values[unknown]-rhs),initial=0)/max(np.max(abs(rhs),initial=0),1))
        if not np.isfinite(values).all() or self.field_residual>1e-8:raise RuntimeError('Dynamic dielectric field residual')
        return values

    def advance(self, until_s, max_events=None):
        """Advance through bridges and ruptures without resetting the random clock."""
        until_s = float(until_s)
        if not np.isfinite(until_s) or until_s < self.time:
            raise ValueError('until_s must be finite and not before current time')
        limit = self.c.max_events if max_events is None else int(max_events)
        if limit < self.steps: raise ValueError('event ceiling precedes current event count')
        reason = 'point_complete'
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
                metal = int(((self.state == 2)&self.interior).sum())
                if metal >= self.snapshot_target:
                    self.save()
                    self.snapshot_target = next((x for x in [25,75,150,300,500,750,1100,1600]
                                                 if x>self.snapshot_target), self.n)
        if self.trace[-1]['events'] != self.steps or self.trace[-1]['time_s'] != self.time:
            self.save()
        return self.result(reason)

class DynamicCachedModel(CatalogueCache,DynamicModel):
    def _build_fixed(self):
        super()._build_fixed()
        if self.first_bridge is None and not self.connected:return
        source,inject,absorb,sites,red,ox,ma,mb,migration=self._catalogue_fixed
        donor=self.valid[sites]&(self.state[self.safe[sites]]==2)&self.cathode[self.safe[sites]]
        potential=self.network['metal_phi_at_1V'].ravel()*self.c.voltage
        electron=(np.where(donor,np.nan_to_num(potential[self.safe[sites]]),0).sum(axis=1)/np.maximum(donor.sum(axis=1),1))
        dg=self.c.metal_standard_eV-self.binding(sites)-self.site[sites]-self.phi[sites]+electron
        pt_only=self.coord[sites]==((self.safe[sites]>=self.n-self.c.nx)&self.valid[sites]).sum(axis=1)
        barrier=self.c.reaction_barrier_eV+self.c.pt_nucleation_eV*pt_only
        red=self.rate(barrier+np.maximum(dg,0));ox=self.rate(barrier+np.maximum(-dg,0))
        self._catalogue_fixed=(source,inject,absorb,sites,red,ox,ma,mb,migration)


def stability(records,transitions,windows=6):
    if len(records)<windows:return dict(assessed=False,filament_stable=False,ions_stationary=False)
    last=records[-windows:];g=np.array([r['network_conductance_S'] for r in last]);metal=np.array([r['metal'] for r in last]);ions=np.array([r['ions'] for r in last])
    grange=float(np.ptp(g)/max(np.median(g),1e-30));mrange=float(np.ptp(metal)/max(np.median(metal),1));irange=float(np.ptp(ions)/max(np.median(ions),1))
    toggles=sum(last[0]['time_s']<r['time_s']<=last[-1]['time_s'] for r in transitions)
    stable=bool(all(r['connected'] for r in last) and toggles==0 and grange<=.05 and mrange<=.01)
    return dict(assessed=True,observation_span_s=last[-1]['time_s']-last[0]['time_s'],filament_stable=stable,
                ions_stationary=bool(irange<=.05 or np.ptp(ions)<=2),conductance_relative_range=grange,
                metal_relative_range=mrange,ion_relative_range=irange,topology_changes=toggles,
                criterion='last six zero-bias observations (five intervals): connected throughout, no toggles, G range <=5%, metal range <=1%; ion range <=5% assessed separately')
