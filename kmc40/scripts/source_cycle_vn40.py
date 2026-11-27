"""Source-conserving surface-event KMC with resolved dielectric supply paths.

Mobile first passage is eliminated between reactions. Stored routes give Ag
packet provenance; they are effective fastest routes, not full random walks.
"""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
from scipy.ndimage import convolve, label
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

import form_hourglass40 as geom
import cycle_hourglass40 as oldcycle

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/vacancy_assisted40/source_cycle'
NX,NY=100,80
KT=.025852
V=.015*np.r_[np.arange(101),np.arange(100,-1,-1)]
KERNEL=geom.KERNEL


def source_routes(metal,v,vn,eta,assist):
    """Positive paired hop rates on a single pre-reaction host geometry."""
    host=~metal
    cn=convolve(metal.astype(int),KERNEL,mode='constant')
    pair=vn&(cn>0)&host
    energy=eta*(-.5-2.48)*pair + v*(1-(np.arange(NY)[:,None]+.5)/NY)
    ids=np.arange(NX*NY).reshape(NY,NX)
    a=[];b=[]
    for x,y,h in [(ids[:-1],ids[1:],host[:-1]&host[1:]),(ids[:,:-1],ids[:,1:],host[:,:-1]&host[:,1:])]:
        a.extend([x[h],y[h]]);b.extend([y[h],x[h]])
    a=np.concatenate(a);b=np.concatenate(b)
    dg=energy.ravel()[b]-energy.ravel()[a]
    base=.24-assist*(pair.ravel()[a]|pair.ravel()[b])
    barrier=base+np.maximum(dg,0)
    # Each edge weight is its conditional mean transition time in microseconds.
    time=1e-6*np.exp(barrier/KT)
    graph=coo_matrix((time,(a,b)),shape=(NX*NY,NX*NY)).tocsr()
    return graph,cn,pair


class Pilot:
    def __init__(self,seed=314,eta=.03,assist=.05,pt_mode='top_front',vn_fraction=.15):
        self.rng=np.random.default_rng(seed)
        self.seed=seed;self.eta=eta;self.assist=assist;self.pt_mode=pt_mode
        self.ids=np.zeros((NY,NX),np.int32)
        self.vn_fraction=vn_fraction
        self.vn=self.rng.random((NY,NX))<vn_fraction
        self.injected=0;self.returned=0;self.events=[];self.paths=[];self.reservoir_paths=[]
        self.snapshots=[];self.records=[]
        self.mu=25.;self.sigma=6.

    @classmethod
    def restore(cls,path):
        path=Path(path)
        data=json.loads((path/'restart.json').read_text(encoding='utf8'))
        model=cls(data['seed'],data['eta'],data['assist'],data['pt_mode'],data['vn_fraction'])
        with np.load(path/'restart.npz') as saved:
            model.ids=saved['metal_id'].copy();model.vn=saved['VN'].copy()
        model.rng.bit_generator.state=data['rng_state']
        for key in ('injected','returned','mu','sigma','form_scale'):
            setattr(model,key,data[key])
        model.resumed_cycle=int(data['phase_end']['cycle'])
        model.restart_source=str(path.resolve())
        return model

    def metallic(self):return self.ids>0

    def choose_source(self):
        metal=self.metallic();top=oldcycle.clusters(metal)[0]
        source=(~metal)&(convolve(top.astype(int),KERNEL,mode='constant')>0)
        source[0]|=~metal[0]
        free=np.flatnonzero(source)
        if not len(free):return None,None
        y,x=np.divmod(free,NX)
        p=np.exp(-.5*(((x+.5)*.5-self.mu)/self.sigma)**2+y/18.)
        node=int(self.rng.choice(free,p=p/p.sum()))
        yy,xx=divmod(node,NX)
        if yy==0:return node,np.array([node],np.int16)
        # Source reservoir access consists of a metal path from the fixed top boundary.
        ids=np.arange(NX*NY).reshape(NY,NX);valid=top.copy();valid[yy,xx]=True
        aa=[];bb=[]
        for a,b,mask in [(ids[:-1],ids[1:],valid[:-1]&valid[1:]),(ids[:,:-1],ids[:,1:],valid[:,:-1]&valid[:,1:])]:
            aa.extend([a[mask],b[mask]]);bb.extend([b[mask],a[mask]])
        aa=np.concatenate(aa);bb=np.concatenate(bb)
        graph=coo_matrix((np.ones(len(aa)),(aa,bb)),shape=(NX*NY,NX*NY)).tocsr()
        dist,pred=dijkstra(graph,indices=node,return_predecessors=True,unweighted=True)
        anchors=np.flatnonzero(top[0]&np.isfinite(dist[:NX]));end=int(anchors[np.argmin(dist[anchors])])
        path=[end]
        while path[-1]!=node:path.append(int(pred[path[-1]]))
        return node,np.array(path,np.int16)

    def deposit(self,v,phase,cycle,point):
        metal=self.metallic();src,reservoir=self.choose_source()
        if src is None:return False
        graph,cn,pair=source_routes(metal,v,self.vn,self.eta,self.assist)
        dist,pred=dijkstra(graph,indices=src,return_predecessors=True)
        front=(~metal)&(cn>0);front[0]=~metal[0]
        if self.pt_mode=='arrival_pt':
            front[-1]|=~metal[-1]
        y,x=np.where(front & np.isfinite(dist.reshape(NY,NX)))
        if not len(y):return False
        # Local metal-neighbour pairing, effective electrochemical free energy.
        adhesion=.08*((y==0)|(y==NY-1))
        dg=.075-.035*cn[y,x]+self.eta*(-.5-2.48)*pair[y,x]-.08*v-adhesion
        barrier=.30+np.maximum(dg,0.)
        # Field focusing and lateral injection profile are coarse geometric factors.
        # No target y or prescribed filament waist enters candidate selection.
        forward=np.zeros(len(y))
        nz=y>0;forward[nz]=metal[y[nz]-1,x[nz]]
        directional=np.where(y==0,1.,1.+2.*forward)
        focus=np.exp(y/13.)
        capture=np.exp(-.5*((x-src%NX)*.5/3.)**2)
        redrate=1e3*np.exp(-barrier/KT)
        travel=dist[y*NX+x]
        w=directional*focus*capture/(1/redrate+travel)
        if not w.sum():return False
        j=int(self.rng.choice(len(y),p=w/w.sum()))
        target=int(y[j]*NX+x[j]);route=[target]
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
        self.paths.append(np.array(route,np.int16))
        self.reservoir_paths.append(reservoir)
        self.events.append(dict(event=len(self.events),phase=phase,cycle=cycle,point=point,voltage_V=v,
            kind='inject_transport_reduce' if accepted else 'inject_unreacted_return',Ag_id=agid,source_x=src%NX,source_y=src//NX,y=int(y[j]),x=int(x[j]),
            path_nodes=len(route),mean_path_us=float(travel[j]),redox_delta_eV=float(dg[j]),
            forward_barrier_eV=float(barrier[j]),reverse_barrier_eV=float(.30+max(-dg[j],0)),
            source_kind='top_connected_Ag_reservoir',delivery_probability=float(acceptance),injected=self.injected,returned=self.returned))
        return accepted

    def oxidize(self,v,phase,cycle,point,thermal_memory=None):
        metal=self.metallic();cn=convolve(metal.astype(int),KERNEL,mode='constant')
        mask=metal&(cn<4)
        t,b,span=oldcycle.clusters(metal)
        if span.any():
            power,_=geom.current_field(metal);memory=power/max(power.max(),1e-30)
        elif thermal_memory is not None:
            memory=thermal_memory
        else:
            from scipy.ndimage import distance_transform_edt
            sep=np.where(t,distance_transform_edt(~b),distance_transform_edt(~t))
            memory=np.exp(-sep/2.)
        dg=.075-.035*cn+self.eta*(-.5-2.48)*(self.vn&(cn>0))-.08*v
        dg[[0,-1]]-=.08
        oxide=.30+np.maximum(-dg,0)
        score=mask*np.exp(-(oxide-.30)/KT)*(.00001+memory**1.8)
        if not score.sum():return False,memory
        # Search only exposed sites that can exit through a continuous dielectric path.
        for _ in range(15):
            target=int(self.rng.choice(NX*NY,p=score.ravel()/score.sum()));y,x=divmod(target,NX)
            opened=metal.copy();opened[y,x]=False
            graph,_,_=source_routes(opened,v,self.vn,self.eta,self.assist)
            dist,pred=dijkstra(graph,indices=target,return_predecessors=True)
            top=oldcycle.clusters(opened)[0]
            target_source=(~opened)&(convolve(top.astype(int),KERNEL,mode='constant')>0)
            target_source[0]|=~opened[0]
            target_source[y,x]=False
            tops=np.flatnonzero(np.isfinite(dist)&target_source.ravel())
            if len(tops):break
            score[y,x]=0
            if not score.sum():return False,memory
        else:return False,memory
        end=int(tops[np.argmin(dist[tops])]);route=[end]
        while route[-1]!=target:route.append(int(pred[route[-1]]))
        route=route[::-1];agid=int(self.ids[y,x]);self.ids[y,x]=0;self.returned+=1
        self.paths.append(np.array(route,np.int16));self.reservoir_paths.append(np.array([],np.int16))
        self.events.append(dict(event=len(self.events),phase=phase,cycle=cycle,point=point,voltage_V=v,
            kind='oxidize_transport_return',Ag_id=agid,source_x=end%NX,source_y=end//NX,y=y,x=x,path_nodes=len(route),
            mean_path_us=float(dist[end]),redox_delta_eV=float(-dg[y,x]),
            forward_barrier_eV=float(oxide[y,x]),reverse_barrier_eV=float(.30+max(dg[y,x],0)),
            source_kind='existing_Ag',delivery_probability=1.,injected=self.injected,returned=self.returned))
        return True,memory

    def checkpoint(self,v,phase,cycle,point,counts):
        metal=self.metallic()
        assert metal.sum()==self.injected-self.returned
        _,unit=geom.current_field(metal)
        if unit:
            g=.04*unit
        else:
            gap=geom.min_gap(metal)
            g=1.e-6+1.e-3*np.exp(-10.25*max(.1,gap))
        cap=.001 if v>=0 else .01
        cur=np.sign(v)*min(cap,abs(v)*g)
        self.snapshots.append(self.ids.copy())
        self.records.append(dict(phase=phase,cycle=cycle,point=point,voltage_V=v,current_A=cur,
            read_G_S=g,Ag_total=int(metal.sum()),Ag_bottom10nm=int(metal[-20:].sum()),
            Ag_bottom_interface=int(metal[-1].sum()),injected=self.injected,returned=self.returned,
            connected=bool(unit>0),added=counts[0],removed=counts[1],VN_adjacent_Ag=int((self.vn&(convolve(metal.astype(int),KERNEL,mode='constant')>0)).sum())))

    def scan(self,phase,cycle,scale,negative=False):
        memory=None
        for j,u in enumerate(V):
            v=-u if negative else u;added=removed=0
            if negative:
                n=int(self.rng.poisson(.16*np.expm1(u/.42)))
                for _ in range(n):
                    ok,memory=self.oxidize(v,phase,cycle,j,memory);removed+=ok
            elif phase=='forming' and j>100:
                for _ in range(self.rng.poisson(.7*np.exp(-u/.5))):
                    ok,memory=self.oxidize(u,phase,cycle,j,memory);removed+=ok
            else:
                _,unit=geom.current_field(self.metallic())
                veff=min(u,.001/(.04*unit)) if unit else u
                n=int(self.rng.poisson(scale*np.expm1(veff/(.3 if phase=='forming' else .11))))
                for _ in range(n):
                    if self.deposit(v,phase,cycle,j):added+=1
            self.checkpoint(v,phase,cycle,j,(added,removed))

    def save(self,path):
        path.mkdir(parents=True,exist_ok=True)
        records=self.records
        args={k:np.array([r[k] for r in records]) for k in records[0]}
        np.savez_compressed(path/'states.npz',metal_id=np.array(self.snapshots),VN=self.vn,**args)
        offsets=np.r_[0,np.cumsum([len(p) for p in self.paths])]
        np.savez_compressed(path/'transport.npz',offsets=offsets,nodes=np.concatenate(self.paths),
                            reservoir_offsets=np.r_[0,np.cumsum([len(p) for p in self.reservoir_paths])],reservoir_nodes=np.concatenate(self.reservoir_paths),
                            Ag_id=np.array([e['Ag_id'] for e in self.events]),kind=np.array([e['kind'] for e in self.events]))
        for filename,rows in [('events.csv',self.events),('iv.csv',records)]:
            with (path/filename).open('w',encoding='utf8',newline='') as f:
                w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
        ends=[r for r in records if r['point']==201]
        settings=dict(seed=self.seed,eta=self.eta,migration_saddle_reduction_eV=self.assist,
            pt_mode=self.pt_mode,pair_formation_pristine_eV=2.48,pair_formation_VN_eV=-.5,
            VN_fraction=self.vn_fraction,temperature_K=300,metal_packet_count='0.5nm cell effective packet, not atom',
            forming_supply_prefactor=self.form_scale,SET_supply_prefactor=.025,RESET_supply_prefactor=.16,
            source='Fixed top electrode connected to exposed active Ag surface through instant metallic-reservoir mass exchange; ions start at adjacent host sites',
            catalogue='Phase-conditioned deposition or oxidation. Local paired barriers are reversible, full phase simulation is driven coarse-graining.',
            interfacial_adhesion_eV=.08,common_field_redox_coefficient=.08,delivery_window_us=1000.,
            transfer='Shortest mean-edge-time dielectric path plus reaction time; success probability 1-exp(-1000us/tau). Failed packets explicitly return. Between events free ions integrated out.',
            phase_ends=ends)
        (path/'summary.json').write_text(json.dumps(settings,indent=2),encoding='utf8')
        checkpoint=dict(rng_state=self.rng.bit_generator.state,seed=self.seed,eta=self.eta,
            assist=self.assist,pt_mode=self.pt_mode,vn_fraction=self.vn_fraction,
            injected=self.injected,returned=self.returned,mu=self.mu,sigma=self.sigma,
            form_scale=self.form_scale,phase_end=ends[-1])
        (path/'restart.json').write_text(json.dumps(checkpoint,indent=2),encoding='utf8')
        np.savez_compressed(path/'restart.npz',metal_id=self.ids,VN=self.vn)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mode',choices=['assisted','binding_only','neutral','arrival_pt'],default='assisted')
    p.add_argument('--scale',type=float,default=2.5)
    p.add_argument('--cycles',type=int,default=3)
    p.add_argument('--vn-fraction',type=float,default=.15)
    p.add_argument('--resume',type=Path)
    a=p.parse_args()
    if not 0<=a.vn_fraction<=1:raise ValueError('VN fraction must lie between 0 and 1.')
    if a.resume:
        model=Pilot.restore(a.resume);first_cycle=model.resumed_cycle+1
        destination=OUT/(a.resume.name+f'_continued_{first_cycle}')
    else:
        model=Pilot(314,eta=0 if a.mode=='neutral' else .03,assist=.05 if a.mode in ['assisted','arrival_pt'] else 0,pt_mode='arrival_pt' if a.mode=='arrival_pt' else 'top_front',vn_fraction=a.vn_fraction)
        model.form_scale=a.scale;first_cycle=1
        model.scan('forming',0,a.scale);print('forming',model.records[-1],flush=True)
        destination=OUT/(a.mode if a.vn_fraction==.15 else f'{a.mode}_VN{a.vn_fraction:g}')
    for cycle in range(first_cycle,first_cycle+a.cycles):
        model.scan('SET',cycle,.025);print('SET',cycle,model.records[-1],flush=True)
        model.scan('RESET',cycle,0,negative=True);print('RESET',cycle,model.records[-1],flush=True)
    model.save(destination)


if __name__=='__main__':main()
