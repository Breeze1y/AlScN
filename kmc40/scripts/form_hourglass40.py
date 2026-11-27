"""Voltage-clocked stochastic interfacial Ag deposition on an initially empty film."""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.ndimage import label, convolve, distance_transform_edt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/hourglass40/forming"
NX,NY,DX=100,80,.5
V=.015*np.r_[np.arange(101),np.arange(100,-1,-1)]
KERNEL=np.array([[0,1,0],[1,0,1],[0,1,0]])
NEIGH=((-1,0),(1,0),(0,-1),(0,1))
PILOT=(2020,2021,2022,2023,2024)
SEEDS=tuple(range(310,322))


def span(mask):
    labs,_=label(mask)
    ids=np.intersect1d(labs[0][labs[0]>0],labs[-1][labs[-1]>0])
    return np.isin(labs,ids)


def arms(mask):
    labs,_=label(mask)
    return [np.isin(labs,np.unique(labs[r][labs[r]>0])) for r in (0,-1)]


def min_gap(mask):
    a,b=arms(mask)
    if not a.any() or not b.any(): return 40.
    if (a&b).any(): return 0.
    # Shortest separation of non-overlapping square Ag domains.
    from scipy.spatial import cKDTree
    aa=np.argwhere(a);bb=np.argwhere(b)
    _,ii=cKDTree(aa).query(bb,k=min(5,len(aa)))
    ii=np.asarray(ii).reshape(len(bb),-1)
    sep=np.maximum(abs(aa[ii]-bb[:,None,:])-1,0)
    return float(np.linalg.norm(sep,axis=2).min()*DX)


def current_field(mask):
    active=span(mask); q=np.zeros_like(mask,float)
    if not active.any(): return q,0.
    ix=np.full(mask.shape,-1,int); ix[active]=np.arange(active.sum())
    h=active[:,1:]&active[:,:-1];z=active[1:]&active[:-1]
    a=np.r_[ix[:,1:][h],ix[1:][z]]; b=np.r_[ix[:,:-1][h],ix[:-1][z]]
    tc=ix[0][active[0]];bc=ix[-1][active[-1]]
    n=int(active.sum()); ends=np.zeros(n);ends[np.r_[tc,bc]]=2.
    di=np.bincount(np.r_[a,b],minlength=n)+ends
    lap=coo_matrix((np.r_[di,-np.ones(2*len(a))],(np.r_[np.arange(n),a,b],np.r_[np.arange(n),b,a])),shape=(n,n)).tocsc()
    src=np.zeros(n);src[tc]=2.;phi=spsolve(lap,src)
    dw=(phi[a]-phi[b])**2
    power=np.bincount(np.r_[a,b],np.r_[dw,dw]*.5,minlength=n)
    power[tc]+=2*(phi[tc]-1)**2;power[bc]+=2*phi[bc]**2
    q[active]=power
    return q,float(power.sum())


class Deposition:
    def __init__(self,seed,top_fraction=.5):
        self.rng=np.random.default_rng(seed)
        self.top_fraction=top_fraction
        top=int(self.rng.integers(24,76))
        self.origins=[top,int(np.clip(round(top+self.rng.normal(0,5.5)),6,NX-7))]
        self.cells=np.zeros((NY+2,NX+2),np.int8)
        self.deposited_at=np.zeros_like(self.cells,dtype=np.int32)
        self.cells[[0,-1],:]=-1;self.cells[:,[0,-1]]=-1
        self.front=[{(1,top+1)},{(NY,self.origins[1]+1)}]
        self.bridge=False
        self.n=0

    def step(self):
        side=int(self.rng.random()>=self.top_fraction)
        c=self.cells
        sites=np.array(sorted(self.front[side]))
        y,x=sites.T;sgn=1 if side==0 else -1
        direct=(c[y-sgn,x]==side+1)
        left=c[y,x-1]==side+1;right=c[y,x+1]==side+1
        n=direct.astype(int)+left.astype(int)+right.astype(int)+(c[y+sgn,x]==side+1)
        # Forward-growth events have a larger rate than side deposition. Concave sites
        # sharing more metal bonds have greater incorporation probability.
        weight=(.95+3.0*direct)*np.exp(.38*np.maximum(n-1,0))
        boundary=(y==(1 if side==0 else NY))
        weight*=np.where(boundary,1.8,1.)
        j=int(np.searchsorted(np.cumsum(weight),self.rng.random()*weight.sum()))
        yy,xx=sites[j]
        c[yy,xx]=side+1
        self.deposited_at[yy,xx]=self.n+1
        for f in self.front: f.discard((yy,xx))
        for dy,dx in NEIGH:
            adj=(yy+dy,xx+dx);value=c[adj]
            if value==0:self.front[side].add(adj)
            elif value==2-side:self.bridge=True
        self.n+=1

    def metal(self): return self.cells[1:-1,1:-1]>0


def grow(seed):
    dep=Deposition(seed);frames=[dep.metal().copy()]
    while not dep.bridge and dep.n<5000:
        dep.step()
        if dep.n%10==0:frames.append(dep.metal().copy())
    if not dep.bridge:raise RuntimeError('No forming bridge')
    return dep,np.stack(frames)


def widths(m):
    s=span(m);profiles=s.sum(1)*DX
    mean=np.convolve(profiles,np.ones(7)/7,'valid')
    pos=int(np.argmin(mean))+3
    return dict(top_width_nm=float(profiles[:16].mean()),mid_width_nm=float(profiles[32:48].mean()),
                pt_width_nm=float(profiles[-16:].mean()),waist_3p5nm_window_depth_nm=(pos+.5)*DX,
                waist_3p5nm_window_width_nm=float(mean.min()),spanning_cells=int(s.sum()))


def run(seed,scale):
    dep=Deposition(seed);rng=np.random.default_rng(seed+80000)
    states=[];events=[];target=rng.exponential()
    f=None;formed=None;opened=False;lost=0;f_at=None
    # Separate zero-field relaxation clock for capillarity of the newly deposited tip.
    for t,v in enumerate(V):
        attached=0;released=0
        if t<=100:
            dose=scale*np.expm1(v/.30) if f is None else .32*v/1.5
            left=1.
            while dose>0 and target<=left*dose:
                left-=target/dose
                dep.step();attached+=1
                if f is None and dep.bridge:
                    f=dep.metal().copy();f_at=t
                if f is not None:dose=.32*v/1.5
                target=rng.exponential()
            if left>0 and dose>0:target-=left*dose
        if t==100:
            formed=dep.metal().copy()
            q,_=current_field(formed);q/=max(q.max(),1.e-12)
            memory=q.copy();m=formed.copy()
        if t>100:
            # Undersaturated young Ag releases packets at low forward bias. Before
            # disconnection current-crowding favours narrow carrying sites. Surface
            # coordination controls detachment energy at both terminal arms.
            for _ in range(rng.poisson(.7*np.exp(-v/.5))):
                nn=convolve(m.astype(int),KERNEL,mode='constant')
                nn[[0,-1]]+=1
                expose=m&(nn<4)&(nn>0)
                rates=expose*np.exp(-2.1*(nn-2))*(.00003+memory**1.8)
                if rates.sum()<=0:continue
                k=np.searchsorted(np.cumsum(rates),rng.random()*rates.sum())
                y,x=divmod(k,NX);m[y,x]=False;lost+=1;released+=1
                if not opened:
                    if not span(m).any():opened=True
                    else:
                        q,_=current_field(m); memory=q/max(q.max(),1.e-12)
            mask=m.copy()
        else:mask=dep.metal().copy()
        states.append(mask);events.append((attached,released))
    if f is None and formed is not None: raise RuntimeError('Selected scale did not form a bridge.')
    return dict(state=np.array(states),F=f,Lformed=formed,H0=states[-1],voltage_V=V,
                P=np.zeros((NY,NX),bool),growth_release_events=np.array(events),
                nuclei_x=np.array(dep.origins),F_index=f_at,
                **{'seed':seed})


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf8') as stream:
        w=csv.DictWriter(stream,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True,parents=True)
    pilot=[]
    for seed in PILOT:
        dep,_=grow(seed);pilot.append(dep.n)
    expected=np.median(pilot);scale=expected/np.expm1(V[:90]/.30).sum()
    results=[];rows=[]
    for seed in SEEDS:
        obj=run(seed,scale);w=widths(obj['Lformed'])
        r=dict(seed=seed,first_bridge_V=V[obj['F_index']],top_x_nm=(obj['nuclei_x'][0]+.5)*.5,
               pt_x_nm=(obj['nuclei_x'][1]+.5)*.5,H0_gap_nm=min_gap(obj['H0']),**w)
        r['hourglass']=bool(min(w['top_width_nm'],w['pt_width_nm'])>=1.5*w['mid_width_nm'] and 14<=w['waist_3p5nm_window_depth_nm']<=26)
        results.append(obj);rows.append(r)
        print(r,flush=True)
    # Select the median mid-band width among ensemble members that re-open and have bilateral flaring.
    eligible=[i for i,r in enumerate(rows) if r['hourglass'] and r['H0_gap_nm']>=.5]
    if not eligible:raise ValueError('No compatible hourglass/open return in this parameter set')
    mid=np.median([rows[i]['mid_width_nm'] for i in eligible])
    winner=min(eligible,key=lambda i:abs(rows[i]['mid_width_nm']-mid))
    obj=results[winner];r=rows[winner]
    np.savez_compressed(OUT/'selected.npz',**obj)
    np.savez_compressed(OUT/'ensemble.npz',seeds=np.array(SEEDS),F=np.stack([o['F'] for o in results]),
                        Lformed=np.stack([o['Lformed'] for o in results]),H0=np.stack([o['H0'] for o in results]))
    write_csv(OUT/'ensemble.csv',rows)
    (OUT/'summary.json').write_text(json.dumps(dict(selected=r,selection='median middle-band width among hourglass and open H0',n_eligible=len(eligible),pilot_events=pilot,scale=scale),indent=2,default=lambda x:x.item()),encoding='utf8')
    analyze()


def analyze():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    with np.load(OUT/'selected.npz') as d:
        frames=d['state']; f=d['F']; formed=d['Lformed']; h0=d['H0']; v=d['voltage_V']
        events=d['growth_release_events']; seed=int(d['seed']); first=int(d['F_index'])
    dep,_=grow(seed); owners=dep.cells[1:-1,1:-1].copy()
    np.testing.assert_array_equal(owners>0,f)
    assert len(np.unique(dep.deposited_at[dep.cells>0]))==f.sum()
    np.savez_compressed(OUT/'origins.npz',F_owner=owners,F_deposition_event=dep.deposited_at[1:-1,1:-1],
                        codes=np.array(['dielectric','Ag-side incorporation','Pt-side reduction']))
    amounts=frames.sum(axis=(1,2))
    np.testing.assert_array_equal(amounts,np.cumsum(events[:,0]-events[:,1]))
    exp_path=next((ROOT.parent/'baseline').glob('SET*forming*.csv'))
    with exp_path.open(encoding='utf-8-sig') as fp:
        raw=np.array([(float(row[1]),float(row[2])) for row in csv.reader(fp) if row and row[0]=='DataValue'])
    np.testing.assert_allclose(raw[:,0],v)
    p=json.loads((ROOT/'results/single_layer_20_40_run/simulation/AlScN40/L1.json').read_text())['params']
    assert not span(frames[-1]).any()
    data=[]; last_state=None
    for i,(voltage,m) in enumerate(zip(v,frames)):
        if last_state is None or not np.array_equal(m,last_state):
            _,unit=current_field(m)
            gap=min_gap(m) if unit==0 else 0.
            last_state=m
        z=(voltage/p['fresh_turn_V'])**2
        leak=p['fresh_G0']+p['fresh_delta_G']*z/(1+z)
        # Scale the metal graph with an effective sheet conductance; near-contact
        # tunnelling is weighted as a single shortest channel in this forming diagnostic.
        g=leak+.1*unit+(.001*np.exp(-10.25*max(gap,.15)) if unit==0 else 0.)
        i_model=min(.001,voltage*g/(1+10*g))
        data.append(dict(index=i,voltage_V=voltage,meas_forming_ref_A=raw[i,1],model_A=i_model,
                         connected=int(unit>0),min_edge_gap_nm=gap,metal_cells=int(m.sum()),
                         added=int(events[i,0]),released=int(events[i,1]),G_ohmic_unit=unit))
    write_csv(OUT/'forming_trace.csv',data)
    summary=json.loads((OUT/'summary.json').read_text())
    j=next(i for i in range(101,202) if data[i]['connected']==0)
    summary['selected_trace']=dict(bridge_point=first,open_return_point=j,open_return_V=v[j],
        F_metal_cells=int(f.sum()),positive_peak_metal_cells=int(formed.sum()),H0_metal_cells=int(h0.sum()),
        H0_direct_spans=0,return_released_cells=int(events[:,1].sum()),reference_500uA_onset_V=1.335,
        mass_balance_every_frame=True,source_reference=str(exp_path.relative_to(ROOT.parent)),
        seeds_with_hourglass=sum(r['hourglass']=='True' for r in csv.DictReader((OUT/'ensemble.csv').open())))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
    plt.rcParams.update({'font.family':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,
                         'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,4,figsize=(13.8,3.85),layout='constrained',gridspec_kw={'width_ratios':[1,1,1,1.26]})
    lead=np.argmin(abs(v[:first]-(v[first]-.12)))
    for ax,stage,title in zip(axs[:3],[frames[lead],f,h0],
             [f'接近贯通  +{v[lead]:.3f} V',f'首次贯通  +{v[first]:.3f} V','形成回零  ·  双端残余 Ag']):
        colors=np.where(stage,1,0)
        if stage is f:colors=owners
        ax.imshow(colors,extent=(0,50,40,0),interpolation='nearest',vmin=0,vmax=2,
                  cmap=ListedColormap(['#fafbfc','#294f6f','#2e8986']))
        ax.set(title=title,xlabel='横向位置（nm）',xticks=[0,10,20,30,40,50])
        ax.axhline(0,color='#8895a0',lw=2);ax.axhline(40,color='#8895a0',lw=2)
    axs[0].set_ylabel('距 Ag 电极（nm）')
    ax=axs[-1];pts=np.array([d['model_A'] for d in data])
    for sl,name,col in [(slice(1,101),'去程','#294f6f'),(slice(101,201),'回程','#bd612e')]:
        ax.scatter(v[sl],np.abs(raw[sl,1]),s=11,lw=.55,facecolors='none',edgecolors=col)
        ax.plot(v[sl],pts[sl],lw=1.7,color=col,label=name)
    ax.set(title='形成 I–V  ·  空心点为共用参考',xlabel='扫描电压（V）',ylabel='|I|（A）',yscale='log')
    ax.legend(frameon=False)
    fig.savefig(OUT/'forming_evolution.png',dpi=175);plt.close(fig)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analyze-only',action='store_true')
    if parser.parse_args().analyze_only:analyze()
    else:main()
