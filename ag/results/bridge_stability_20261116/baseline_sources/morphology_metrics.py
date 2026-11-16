"""Explicit lattice morphology definitions; no inference of experimental atoms."""
from collections import deque
import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.csgraph import maximum_flow

CROSS=np.array([[0,1,0],[1,1,1],[0,1,0]])


def thin(mask):
    """Two-subiteration topology thinning using an eight-neighbor ring."""
    image=np.pad(mask.astype(bool),1)
    changed=True
    while changed:
        changed=False
        for phase in [0,1]:
            p=[image[:-2,1:-1],image[:-2,2:],image[1:-1,2:],image[2:,2:],
               image[2:,1:-1],image[2:,:-2],image[1:-1,:-2],image[:-2,:-2]]
            count=sum(x.astype(int) for x in p)
            transitions=sum((~p[i]&p[(i+1)%8]).astype(int) for i in range(8))
            if phase==0: constraint=~(p[0]&p[2]&p[4]) & ~(p[2]&p[4]&p[6])
            else: constraint=~(p[0]&p[2]&p[6]) & ~(p[0]&p[4]&p[6])
            remove=image[1:-1,1:-1]&(count>=2)&(count<=6)&(transitions==1)&constraint
            if remove.any(): image[1:-1,1:-1][remove]=False;changed=True
    return image[1:-1,1:-1]


def _neighbors(index,nz,nx):
    z,x=divmod(index,nx)
    for dz,dx in [(-1,0),(1,0),(0,-1),(0,1)]:
        zz,xx=z+dz,x+dx
        if 0<=zz<nz and 0<=xx<nx: yield zz*nx+xx


def missing_cells(metal):
    nz,nx=metal.shape;flat=metal.ravel();dist=np.full(flat.size,flat.size,int)
    q=deque(range(nx));dist[:nx]=0
    while q:
        i=q.popleft()
        for j in _neighbors(i,nz,nx):
            cost=int(not flat[j]);candidate=dist[i]+cost
            if candidate<dist[j]:
                dist[j]=candidate
                if cost: q.append(j)
                else: q.appendleft(j)
    return int(dist[-nx:].min())


def metal_path_edges(metal):
    nz,nx=metal.shape;flat=metal.ravel();dist=np.full(flat.size,-1,int)
    q=deque(range(nx));dist[:nx]=0
    while q:
        i=q.popleft()
        if i>=(nz-1)*nx: return int(dist[i])
        for j in _neighbors(i,nz,nx):
            if flat[j] and dist[j]<0: dist[j]=dist[i]+1;q.append(j)
    return None


def vertex_bottleneck(metal):
    nz,nx=metal.shape;ids=np.flatnonzero(metal.ravel())
    mapping={int(i):k for k,i in enumerate(ids)}
    source=2*len(ids);sink=source+1;capacity=len(ids)+1
    rows=[];cols=[];values=[]
    def edge(a,b,c): rows.append(a);cols.append(b);values.append(c)
    for i in ids:
        k=mapping[int(i)];z=i//nx
        edge(2*k,2*k+1,capacity if z in [0,nz-1] else 1)
        if z==0: edge(source,2*k,capacity)
        if z==nz-1: edge(2*k+1,sink,capacity)
        for j in _neighbors(int(i),nz,nx):
            if j in mapping: edge(2*k+1,2*mapping[j],capacity)
    graph=sparse.csr_matrix((values,(rows,cols)),shape=(sink+1,)*2,dtype=np.int64)
    return int(maximum_flow(graph,source,sink).flow_value)


def metrics(state, spacing_nm, layer_cells):
    metal=state==2;inside=metal[1:-1];depth=inside.shape[0]*spacing_nm
    labels,_=ndimage.label(metal,CROSS)
    cathode=labels==labels[-1,0];connected=bool(cathode[0].any())
    skeleton=thin(inside)
    degree=ndimage.convolve(skeleton.astype(int),np.ones((3,3),int),mode='constant')-skeleton
    _,junctions=ndimage.label(skeleton&(degree>=3),np.ones((3,3),int))
    path=metal_path_edges(metal) if connected else None
    return dict(metal=int(inside.sum()),ions=int(np.sum(state[1:-1]==1)),connected_4=connected,
                cathode_connected_metal=int(cathode[1:-1].sum()),
                minimum_missing_cells_4=missing_cells(metal),
                vertex_bottleneck_4=vertex_bottleneck(metal) if connected else 0,
                shortest_path_nm=(path-1)*spacing_nm if path is not None else None,
                path_tortuosity=(path-1)*spacing_nm/depth if path is not None else None,
                skeleton_pixels_8=int(skeleton.sum()),skeleton_endpoints_8=int(np.sum(skeleton&(degree==1))),
                skeleton_junction_clusters_8=int(junctions),
                width_by_depth_nm=(inside.sum(axis=1)*spacing_nm).tolist(),
                layer_metals=[int(inside[z*layer_cells:(z+1)*layer_cells].sum()) for z in range(3)],
                layer_ions=[int(np.sum(state[1+z*layer_cells:1+(z+1)*layer_cells]==1)) for z in range(3)])
