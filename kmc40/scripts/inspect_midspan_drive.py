"""Independent depth audit on the saved SET morphology. All drives are normalized."""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.ndimage import label, binary_erosion
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[1]
STATE = ROOT/"results/wide100_vacancy40/growth/AlScN40/set_states.npz"
OUT = ROOT/"results/midspan_reset40/drive"
CELL = .5


def edges(mask):
    ix = np.full(mask.shape, -1, dtype=int)
    ix[mask] = np.arange(mask.sum())
    h = mask[:,1:] & mask[:,:-1]
    v = mask[1:] & mask[:-1]
    a = np.r_[ix[:,:-1][h],ix[:-1][v]]
    b = np.r_[ix[:,1:][h],ix[1:][v]]
    return ix, a, b


def laplacian(n, a, b, g, end):
    return sparse.coo_matrix((np.r_[g,g,-g,-g,end], (np.r_[a,b,a,b,np.arange(n)],np.r_[a,b,b,a,np.arange(n)])),shape=(n,n)).tocsc()


def electrical(span):
    ix, a, b = edges(span)
    n = span.sum()
    contacts = [ix[y][span[y]] for y in (0,-1)]
    end = np.zeros(n)
    end[np.r_[*contacts]]=2.
    rhs = np.zeros(n)
    rhs[contacts[0]]=2.
    V = spsolve(laplacian(n,a,b,np.ones(len(a)),end),rhs)
    power = (V[a]-V[b])**2
    Q = np.bincount(np.r_[a,b],np.r_[power,power]/2,minlength=n)
    for voltage,index in zip((1.,0.),contacts):
        Q[index] += 2*(voltage-V[index])**2
    np.testing.assert_allclose(Q.sum(),2*np.sum(1-V[contacts[0]]), rtol=1e-9)
    voltage = np.zeros(span.shape)
    voltage[span] = V
    heat = np.zeros(span.shape)
    heat[span] = Q
    return voltage,heat


def temperature(metal, Q, ratio):
    """Finite-volume 2-D heat conduction; no lateral flux, isothermal electrode faces."""
    ix,a,b=edges(np.ones_like(metal))
    k = np.where(metal,ratio,1.).ravel()
    g = 2/(1/k[a]+1/k[b])
    bc=np.zeros(k.shape)
    bc[ix[[0,-1]]]=2*k[ix[[0,-1]]]
    t=spsolve(laplacian(k.size,a,b,g,bc),Q.ravel())
    np.testing.assert_allclose(bc@t,Q.sum(),rtol=2e-9)
    return t.reshape(Q.shape)/t.max()


def peak(array, mask):
    i = np.unravel_index(np.argmax(np.where(mask,array,-np.inf)),array.shape)
    return dict(x_nm=(i[1]+.5)*CELL,depth_nm=(i[0]+.5)*CELL,value=float(array[i]))


def connected_metal(metal):
    labs,_=label(metal)
    ids=np.intersect1d(labs[0][labs[0]>0],labs[-1][labs[-1]>0])
    return np.isin(labs,ids)


def recede(metal,ratio):
    # A deterministic high-selectivity limit. Ranking uses no depth coordinate.
    current=metal.copy()
    masks=[current.copy()]
    removed=[]
    drives=[]
    conductance=[]
    while np.any(span:=connected_metal(current)):
        _,Q=electrical(span)
        conductance.append(Q.sum())
        drive=temperature(current,Q,ratio) if ratio is not None else Q/Q.max()
        exposed=current & ~binary_erosion(current, border_value=1)
        eligible=exposed & span & (Q > 1e-8*Q.max())
        if not np.any(eligible):
            raise RuntimeError("No exposed conducting cell")
        where=np.unravel_index(np.argmax(np.where(eligible,drive,-np.inf)),current.shape)
        drives.append(drive)
        current[where]=False
        removed.append(where)
        masks.append(current.copy())
    return np.asarray(removed),np.stack(masks),np.stack(drives),np.r_[conductance,0.]/conductance[0]


def plot_summary(metal,span,heat,temp,removed,masks):
    plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"axes.unicode_minus":False})
    fig,axs=plt.subplots(1,3,figsize=(12,4.4),layout="constrained")
    extent=(0,50,40,0)
    for ax,image,cmap,title in [(axs[0],np.ma.masked_where(~span,heat/heat.max()),"inferno","贯通簇上的相对焦耳热源"),
                                (axs[1],temp[3],"inferno","电极等温时的相对温升")]:
        im=ax.imshow(image,extent=extent,vmin=0,vmax=1,interpolation="nearest",cmap=cmap)
        fig.colorbar(im,ax=ax,shrink=.75)
        ax.set_title(title,size=12)
    xv=(np.arange(metal.shape[1])+.5)*CELL;yv=(np.arange(metal.shape[0])+.5)*CELL
    axs[1].contour(xv,yv,span,levels=[.5],colors=["#94beca"],linewidths=.30)
    ax=axs[2]
    ax.imshow(masks[-1],extent=extent,cmap="Blues",vmin=0,vmax=1,interpolation="nearest")
    for step,(y,x) in enumerate(removed,1):
        cx,cy=(x+.5)*CELL,(y+.5)*CELL
        ax.add_patch(plt.Rectangle((cx-.25,cy-.25),.5,.5,facecolor="#e56d23"))
        delta=[(-22,-18),(-11,23),(14,-19),(18,16)][step-1]
        ax.annotate(str(step),(cx,cy),xytext=delta,textcoords="offset points",color="#992b06",
                    weight='bold',arrowprops={'arrowstyle':'-','color':'#883b23','lw':1.1})
    ax.set(xlim=(12.25,21),ylim=(23,15.5),title="四次局部退缩后的金属形貌")
    for ax in axs:
        ax.set_xlabel("横向位置（nm）")
        ax.set_ylabel("距 Ag 电极（nm）")
    fig.savefig(OUT/"electrothermal_selection.png",dpi=180)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True,parents=True)
    with np.load(STATE) as z:
        metal=z["L1"][1:-1]==2
    span=connected_metal(metal)
    voltage,heat=electrical(span)
    depth=(np.arange(metal.shape[0])+.5)*CELL
    mid = (depth>=18)&(depth<=22)
    ratios = [1,3,10,30,100]
    summary = dict(source=str(STATE.relative_to(ROOT)),number_spanning_cells=int(span.sum()),normalized_joule_peak=peak(heat/heat.max(),span),thermal_cases=[])
    rows=[]
    for ix,d in enumerate(np.arange(1,metal.shape[0])*CELL):
        vertical=span[ix]&span[ix+1]
        currents=np.abs((voltage[ix]-voltage[ix+1])[vertical])
        eff=currents.sum()**2/(currents@currents)
        rows.append(dict(face_depth_nm=d,metal_crossings=int(vertical.sum()),effective_parallel_bonds=eff,current_weighted_power=np.sum(currents**2)))
    temp={}
    for ratio in ratios:
        T=temperature(metal,heat,ratio)
        temp[ratio]=T
        center_hotspot=peak(T,span&mid[:,None])
        summary['thermal_cases'].append(dict(k_Ag_over_matrix=ratio,hottest_silver_site=peak(T,span),
                hottest_middle_silver_site=center_hotspot,center_to_hotspot=center_hotspot['value']/T[span].max()))
    for i,r in enumerate(rows):
        for ratio in ratios:
            v=span[i]&span[i+1]
            r[f"max_relative_temperature_k{ratio}"]=float(((temp[ratio][i]+temp[ratio][i+1])[v]/2).max())
    with (OUT/'depth_scan.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    with (OUT/'hotspots.json').open('w') as f:
        json.dump(summary,f,indent=2)
    np.savez_compressed(OUT/'field_maps.npz',spanning=span,metal=metal,voltage_on_span=voltage,
                        joule_source_normalized=heat/heat.max(),conductivity_ratios=ratios,
                        relative_temperature=np.stack([temp[k] for k in ratios]))
    released={}
    events=[]
    for key in [None,*ratios]:
        removed,masks,drives,rel_g=recede(metal,key)
        assert all(connected_metal(st).any() for st in masks[:-1]) and not connected_metal(masks[-1]).any()
        name="power" if key is None else f"thermal_k{key}"
        released[f"{name}_row_col"]=removed
        released[f"{name}_metal"]=masks
        released[f"{name}_drive_before"]=drives
        released[f"{name}_relative_direct_G"]=rel_g
        for step,(y,x) in enumerate(removed,1):
            events.append(dict(drive=name,step=step,depth_nm=(y+.5)*CELL,x_nm=(x+.5)*CELL,
                               relative_drive=drives[step-1,y,x],
                               direct_G_before_relative=rel_g[step-1],direct_G_after_relative=rel_g[step],
                               remains_spanning=bool(step<len(removed))))
        if key==3:
            plot_summary(metal,span,heat,temp,removed,masks)
    np.savez_compressed(OUT/"retreat_sequences.npz",**released)
    with (OUT/'retreat_events.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=events[0]);w.writeheader();w.writerows(events)
    print(json.dumps(summary))
    return summary


if __name__ == '__main__':
    main()
