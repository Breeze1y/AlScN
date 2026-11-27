"""Exact catalogue caching between metal/voltage changes; no rate approximation."""
import numpy as np
from kmc_model import Model
from forming_iv import FormingModel


class CatalogueCache:
    def update_structure(self):
        super().update_structure()
        self._catalogue_fixed=None

    def _build_fixed(self):
        c,s=self.c,self.state
        source=np.arange(c.nx,2*c.nx)
        dg=self.site[source]+self.phi[source]-c.voltage
        inject=self.rate(c.source_barrier_eV+np.maximum(dg,0))
        absorb=self.rate(c.source_barrier_eV+np.maximum(-dg,0))
        reachable=(self.cathode[self.safe]&self.valid).any(axis=1)
        if c.surface_redox_only:reachable &= (self.valid&(s[self.safe]!=2)).any(axis=1)
        sites=np.flatnonzero(self.interior&reachable)
        dg=c.metal_standard_eV-self.binding(sites)-self.site[sites]-self.phi[sites]
        pt_only=self.coord[sites]==((self.safe[sites]>=self.n-c.nx)&self.valid[sites]).sum(axis=1)
        base=c.reaction_barrier_eV+c.pt_nucleation_eV*pt_only
        reduction=self.rate(base+np.maximum(dg,0))
        oxidation=self.rate(base+np.maximum(-dg,0))
        metals=np.flatnonzero(self.interior&(s==2)&(self.coord<4))
        a=np.repeat(metals,4);b=self.neighbors[metals].ravel()
        ok=(b>=0)&self.interior[np.maximum(b,0)]&(s[np.maximum(b,0)]!=2)
        a,b=a[ok],b[ok]
        dg=self.binding(a)-(self.binding(b)-c.bond_eV)
        migration=self.rate(c.migration_eV+np.maximum(dg,0))
        # Host energies and the field are fixed until the next structural or
        # voltage update. Only occupancy filters change after an ion event.
        ha=np.repeat(np.arange(self.n),4)
        hb=self.safe.ravel()
        hdg=self.site[hb]-self.site[ha]+self.phi[hb]-self.phi[ha]
        hbase=np.where(self.layer[ha]==self.layer[hb],
            np.where(self.layer[ha]==1,c.hop_alscn_eV,c.hop_aln_eV),
            (c.hop_aln_eV+c.hop_alscn_eV)/2+c.interface_eV)
        self._hop_rates=self.rate(hbase+np.maximum(hdg,0)).reshape(self.n,4)
        self._catalogue_fixed=(source,inject,absorb,sites,reduction,oxidation,a,b,migration)
        self._catalogue_config=tuple(vars(c).values())

    def catalogue(self):
        if self._catalogue_fixed is None or self._catalogue_config!=tuple(vars(self.c).values()):
            self._build_fixed()
        c,s=self.c,self.state;groups=[]
        source,inject,absorb,sites,reduction,oxidation,ma,mb,migration=self._catalogue_fixed
        def add(kind,a,b,rates):
            if a.size:groups.append((kind,a,np.broadcast_to(b,a.shape),rates))
        ok=s[source]==0;add('inject',source[ok],-1,inject[ok])
        ok=s[source]==1;add('absorb',source[ok],-1,absorb[ok])
        ions=np.flatnonzero(s==1)
        a=np.repeat(ions,4);b=self.neighbors[ions].ravel()
        ok=(b>=0)&(s[np.maximum(b,0)]==0);a,b=a[ok],b[ok]
        add('hop',a,b,self._hop_rates[ions].ravel()[ok])
        ok=s[sites]==1;add('reduce',sites[ok],-1,reduction[ok])
        ok=s[sites]==2;add('oxidize',sites[ok],-1,oxidation[ok])
        ok=s[mb]==0;add('migrate',ma[ok],mb[ok],migration[ok])
        return groups


class CachedModel(CatalogueCache,Model):
    pass


class CachedFormingModel(CatalogueCache,FormingModel):
    pass
