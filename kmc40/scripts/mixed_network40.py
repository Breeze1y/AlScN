"""Passive Ag/contact/tunnel/defect network, used on both sides of rupture.

Each occupied Ag cell is a sphere of diameter dx for distance purposes; only
orthogonal touching Ag cells have a metallic bond. A diagonal pair therefore
has a positive tunnelling distance. All Ag islands participate, including
islands without a metallic connection to either electrode. Defect nodes carry
an effective electronic conductance, not ionic drift current.

The circuit is linear at each frozen electronic activity. Its potentials are
conduction-node potentials. Host pixels without a node get a linear background
potential in the returned map; this is not a Poisson solution in the host.
"""
from dataclasses import dataclass
import time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve


@dataclass(frozen=True)
class NetworkParams:
    metal_bond_S: float = .04
    tunnel_prefactor_S: float = .001
    tunnel_beta_nm: float = 10.25
    trap_prefactor_S: float = .0001
    trap_beta_nm: float = 2.
    leakage_S: float = 1.e-6


def prepare(metal, vn, dx_nm=.5, radius_nm=1.5):
    """Build geometry once; numerical conductances can be changed in solve()."""
    metal = np.asarray(metal, dtype=bool)
    vn = np.asarray(vn, dtype=bool)
    if metal.ndim != 2 or metal.shape != vn.shape:
        raise ValueError('metal and vn must have the same two-dimensional shape')
    if dx_nm <= 0 or radius_nm < dx_nm:
        raise ValueError('positive dx and radius >= dx are required')
    coords = np.argwhere(metal | vn)
    n = len(coords)
    is_metal = metal[tuple(coords.T)] if n else np.empty(0, bool)
    pairs = cKDTree(coords * dx_nm).query_pairs(radius_nm + 1.e-10, output_type='ndarray') if n else np.empty((0, 2), int)
    a, b = pairs.T if len(pairs) else (np.empty(0, int), np.empty(0, int))
    distance = np.linalg.norm(coords[a]-coords[b], axis=1) * dx_nm
    both = is_metal[a] & is_metal[b]
    touch = both & (np.abs(distance-dx_nm) < 1.e-10)
    # kinds: 0 metallic contact; 1 direct tunnel; 2 defect-assisted edge.
    kind = np.where(touch, 0, np.where(both, 1, 2)).astype(np.int8)
    radii = .5 * dx_nm * is_metal
    gap = np.maximum(distance-radii[a]-radii[b], 0.)
    # Explicit electrode links for all nodes within the geometric cutoff.
    y = (coords[:, 0]+.5)*dx_nm if n else np.empty(0)
    top_gap = np.maximum(y-radii, 0.)
    bottom_gap = np.maximum(metal.shape[0]*dx_nm-y-radii, 0.)
    top = np.flatnonzero(y <= radius_nm+1.e-10)
    bottom = np.flatnonzero(metal.shape[0]*dx_nm-y <= radius_nm+1.e-10)
    def boundary_kind(nodes, gaps):
        return np.where(is_metal[nodes], np.where(gaps[nodes] < 1.e-10, 0, 1), 2).astype(np.int8)
    return dict(shape=metal.shape, dx_nm=float(dx_nm), radius_nm=float(radius_nm),
                coords=coords, is_metal=is_metal, a=a, b=b, kind=kind, gap=gap,
                top=top, bottom=bottom, top_gap=top_gap[top], bottom_gap=bottom_gap[bottom],
                top_kind=boundary_kind(top, top_gap), bottom_kind=boundary_kind(bottom, bottom_gap))


def _conductances(kind, gap, p, activity, boundary=False):
    values = np.empty(len(kind), float)
    m, t, d = kind == 0, kind == 1, kind == 2
    values[m] = p.metal_bond_S * (2. if boundary else 1.)
    values[t] = p.tunnel_prefactor_S*np.exp(-p.tunnel_beta_nm*gap[t])
    values[d] = p.trap_prefactor_S*activity*np.exp(-p.trap_beta_nm*gap[d])
    return values


def solve(topology, params=None, activity=1., fields=False):
    """Return G at one volt, including a separate terminal leakage branch.

    activity is the effective fraction/strength of conducting defect states.
    It scales every defect edge, before and after metallic disconnection.
    `unit_power` assigns each internal edge's heat half to each endpoint;
    electrode-interface heat is assigned to its material node. The sum equals
    the material-network G. Leakage is not included in this spatial map.
    """
    p = NetworkParams() if params is None else params
    if not 0 <= activity <= 1 or not np.isfinite(activity):
        raise ValueError('activity must be finite and in [0,1]')
    if min(p.metal_bond_S, p.tunnel_prefactor_S, p.trap_prefactor_S,
           p.leakage_S, p.tunnel_beta_nm, p.trap_beta_nm) < 0:
        raise ValueError('conductances and attenuation coefficients must be nonnegative')
    t = topology
    n = len(t['coords'])
    a, b = t['a'], t['b']
    g = _conductances(t['kind'], t['gap'], p, activity)
    gt = _conductances(t['top_kind'], t['top_gap'], p, activity, True)
    gb = _conductances(t['bottom_kind'], t['bottom_gap'], p, activity, True)
    top, bottom = t['top'], t['bottom']
    phi = np.full(n, .5)
    numerical_residual = 0.
    # Components attached to one terminal have exactly constant potential;
    # unattached components have arbitrary potential and carry no current.
    active = g > 0
    adj = coo_matrix((np.ones(2*active.sum(), bool),
                      (np.r_[a[active], b[active]], np.r_[b[active], a[active]])), shape=(n, n)).tocsr()
    count, labels = connected_components(adj, directed=False)
    has_top = np.zeros(count, bool); has_bottom = np.zeros(count, bool)
    has_top[labels[top[gt > 0]]] = True
    has_bottom[labels[bottom[gb > 0]]] = True
    bridge = has_top & has_bottom
    useful = bridge[labels]
    phi[has_top[labels] & ~has_bottom[labels]] = 1.
    phi[has_bottom[labels] & ~has_top[labels]] = 0.
    if useful.any():
        ix = np.full(n, -1, int); ix[useful] = np.arange(useful.sum())
        edges = useful[a] & useful[b] & active
        aa, bb, gg = ix[a[edges]], ix[b[edges]], g[edges]
        ti, bi = useful[top] & (gt > 0), useful[bottom] & (gb > 0)
        tt, bt = ix[top[ti]], ix[bottom[bi]]
        gtt, gbt = gt[ti], gb[bi]
        nn = useful.sum()
        diagonal = np.bincount(np.r_[aa, bb, tt, bt], np.r_[gg, gg, gtt, gbt], minlength=nn)
        source = np.bincount(tt, gtt, minlength=nn)
        scale = diagonal.max()
        lap = coo_matrix((np.r_[diagonal, -gg, -gg]/scale,
                          (np.r_[np.arange(nn), aa, bb], np.r_[np.arange(nn), bb, aa])), shape=(nn, nn)).tocsc()
        solution = spsolve(lap, source/scale, permc_spec='MMD_AT_PLUS_A')
        if not np.isfinite(solution).all():
            raise FloatingPointError('nonfinite network potential')
        phi[useful] = solution
        if fields:
            numerical_residual = float(np.max(np.abs(lap@solution-source/scale))*scale)
    power_edges = g*(phi[a]-phi[b])**2
    power_top = gt*(1.-phi[top])**2
    power_bottom = gb*phi[bottom]**2
    material_g = float(power_edges.sum()+power_top.sum()+power_bottom.sum())
    result = dict(conductance_S=material_g+p.leakage_S, material_conductance_S=material_g)
    if fields:
        shape = t['shape']; coords = t['coords']
        pmap = np.zeros(shape)
        power = np.bincount(np.r_[a,b,top,bottom], np.r_[.5*power_edges,.5*power_edges,power_top,power_bottom], minlength=n)
        pmap[tuple(coords.T)] = power
        phimap = np.broadcast_to(1.-(np.arange(shape[0])[:,None]+.5)/shape[0], shape).copy()
        phimap[tuple(coords.T)] = phi
        result.update(unit_phi=phimap, unit_power=pmap,
                      top_current_A=float(np.dot(gt,1.-phi[top]))+p.leakage_S,
                      bottom_current_A=float(np.dot(gb,phi[bottom]))+p.leakage_S,
                      max_kcl_residual_A=numerical_residual, nodes=n,
                      component_count=int(count), conducting_nodes=int(useful.sum()),
                      metal_edge_power_W=float(power_edges[t['kind']==0].sum()),
                      tunnel_edge_power_W=float(power_edges[t['kind']==1].sum()),
                      defect_edge_power_W=float(power_edges[t['kind']==2].sum()),
                      contact_power_W=float(power_top.sum()+power_bottom.sum()))
    return result


def read(metal, vn, params=None, activity=1., fields=False, **kwargs):
    return solve(prepare(metal,vn,**kwargs),params,activity,fields)
