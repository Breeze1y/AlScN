"""Kirchhoff network for local, fractional-volume retreat of existing Ag.

An Ag lattice cell is an effective volume: ``fraction`` is its remaining
metal volume, not a separately resolved atom count. Sphere radius scales as
f**(1/3), and an orthogonal contact's effective conductance scales as
min(f_i, f_j)**contact_exponent. Metal, tunnel and exposed-VN channels share
one local potential. Their conductances, rather than the resulting current,
change continuously as Ag retreats.

At every binary fraction map this model reproduces ``mixed_network40`` for
the surviving Ag mask. In particular, a full metallic contact carries no
additional parallel tunnel channel. No extra short-gap prefactor or activity
enhancement is introduced.
"""
from dataclasses import asdict

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

from mixed_network40 import NetworkParams


def prepare(initial_metal, vn, dx_nm=.5, radius_nm=1.5):
    """Retain possible Ag/VN nodes and links for a retreat-only trajectory."""
    initial_metal = np.asarray(initial_metal, dtype=bool)
    vn = np.asarray(vn, dtype=bool)
    if initial_metal.ndim != 2 or initial_metal.shape != vn.shape:
        raise ValueError('Ag and VN arrays must have the same 2-D shape')
    if not np.isfinite([dx_nm, radius_nm]).all() or dx_nm <= 0 or radius_nm < dx_nm:
        raise ValueError('Finite positive dx and radius >= dx are required')
    coords = np.argwhere(initial_metal | vn)
    n = len(coords)
    pairs = (cKDTree(coords * dx_nm).query_pairs(radius_nm + 1.e-10, output_type='ndarray')
             if n else np.empty((0, 2), int))
    a, b = pairs.T if len(pairs) else (np.empty(0, int), np.empty(0, int))
    distance = np.linalg.norm(coords[a] - coords[b], axis=1) * dx_nm
    y = (coords[:, 0] + .5) * dx_nm
    bottom_distance = initial_metal.shape[0] * dx_nm - y
    top = np.flatnonzero(y <= radius_nm + 1.e-10)
    bottom = np.flatnonzero(bottom_distance <= radius_nm + 1.e-10)
    return dict(
        shape=initial_metal.shape, dx_nm=float(dx_nm), radius_nm=float(radius_nm),
        coords=coords, initial_ag=initial_metal[tuple(coords.T)],
        vn=vn[tuple(coords.T)], a=a, b=b, distance=distance,
        touch=np.abs(distance - dx_nm) < 1.e-10,
        top=top, bottom=bottom, top_distance=y[top],
        bottom_distance=bottom_distance[bottom],
        top_touch=np.abs(y[top] - .5 * dx_nm) < 1.e-10,
        bottom_touch=np.abs(bottom_distance[bottom] - .5 * dx_nm) < 1.e-10)


def edge_channels(topology, fraction, params=None, activity=1., contact_exponent=2.,
                  defect_site_gain=None, extra_defect_weight=None):
    """Return (metal, tunnel, defect) arrays for internal/top/bottom links.

    ``VN * (1-f)`` is the exposed electronic defect fraction of a cell.
    The residual tunnel weight ``1-contact_weight`` removes double counting
    at a full metal contact and opens the tunnel contribution during retreat.
    Fraction values outside the original Ag mask do not create new Ag nodes.
    """
    t = topology
    p = NetworkParams() if params is None else params
    fraction = np.asarray(fraction, dtype=float)
    if (fraction.shape != t['shape'] or not np.isfinite(fraction).all()
            or np.any((fraction < 0.) | (fraction > 1.))):
        raise ValueError('Fraction must match the lattice shape and lie in [0,1]')
    if not np.isfinite(activity) or not 0 <= activity <= 1:
        raise ValueError('Activity must be finite and in [0,1]')
    if not np.isfinite(contact_exponent) or contact_exponent <= 0:
        raise ValueError('Contact exponent must be finite and positive')
    values = np.asarray(list(asdict(p).values()), dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError('Network parameters must be finite and nonnegative')
    f = fraction[tuple(t['coords'].T)] * t['initial_ag']
    d = t['vn'] * (1. - f)
    if extra_defect_weight is not None:
        extra = np.asarray(extra_defect_weight, dtype=float)
        if extra.shape != d.shape or not np.all((extra >= 0.) & (extra <= 1.)):
            raise ValueError('Interface state weights must lie in [0,1] for every network node')
        d = d + (~t['vn']) * (1.-f) * extra
    r = .5 * t['dx_nm'] * np.cbrt(f)
    a, b, dist = t['a'], t['b'], t['distance']
    gap = np.maximum(dist - r[a] - r[b], 0.)
    contact = t['touch'] * np.minimum(f[a], f[b]) ** contact_exponent
    metal = p.metal_bond_S * contact
    tunnel = (p.tunnel_prefactor_S * f[a] * f[b]
              * np.exp(-p.tunnel_beta_nm * gap) * (1. - contact))
    # A local change of hopping barrier scales only defect-assisted links.
    # A defect pair sees the mean barrier change of its two endpoints.
    d_contact = d
    d_pair = d
    if defect_site_gain is not None:
        gain = np.asarray(defect_site_gain, dtype=float)
        if gain.shape != d.shape or not np.isfinite(gain).all() or np.any(gain <= 0):
            raise ValueError('Defect gain must be positive and provided for every network node')
        d_contact = d * gain
        d_pair = d * np.sqrt(gain)
    defect = activity * p.trap_prefactor_S * (
        f[a] * d_contact[b] * np.exp(-p.trap_beta_nm * np.maximum(dist - r[a], 0.))
        + d_contact[a] * f[b] * np.exp(-p.trap_beta_nm * np.maximum(dist - r[b], 0.))
        + d_pair[a] * d_pair[b] * np.exp(-p.trap_beta_nm * dist))

    def boundary(nodes, distance, touch):
        contact = touch * f[nodes] ** contact_exponent
        gap = np.maximum(distance - r[nodes], 0.)
        metal = 2. * p.metal_bond_S * contact
        tunnel = (p.tunnel_prefactor_S * f[nodes]
                  * np.exp(-p.tunnel_beta_nm * gap) * (1. - contact))
        defect = activity * p.trap_prefactor_S * d_contact[nodes] * np.exp(-p.trap_beta_nm * distance)
        return np.asarray((metal, tunnel, defect))

    return (np.asarray((metal, tunnel, defect)),
            boundary(t['top'], t['top_distance'], t['top_touch']),
            boundary(t['bottom'], t['bottom_distance'], t['bottom_touch']))


def solve(topology, fraction, params=None, activity=1., contact_exponent=2., fields=False):
    """Solve the passive network at unit terminal bias.

    ``conductance_S`` includes the separate terminal leakage branch;
    ``material_conductance_S`` and all spatial/channel powers exclude it.
    With ``fields=True``, channel powers include their electrode links; the
    ``*_edge_power_W`` entries include only internal links, matching the old
    solver's convention. ``contact_power_W`` is the electrode subtotal.
    """
    p = NetworkParams() if params is None else params
    channels = edge_channels(topology, fraction, p, activity, contact_exponent)
    return solve_channels(topology, channels, p.leakage_S, fields)


def solve_channels(topology, link_channels, leakage_S=0., fields=False):
    """Unit-bias KCL and power readout for supplied nonnegative link conductances.

    The three link arrays (internal, top, bottom) each have shape (3, n)
    for metal, direct tunnelling and defect transport. Parallel links between
    the same nodes are allowed and dissipate power at the common voltage drop.
    """
    t = topology
    channels, ct, cb = link_channels
    g, gt, gb = channels.sum(axis=0), ct.sum(axis=0), cb.sum(axis=0)
    a, b, top, bottom = t['a'], t['b'], t['top'], t['bottom']
    n = len(t['coords'])
    phi = np.full(n, .5)
    numerical_residual = 0.
    active = g > 0.
    adj = coo_matrix((np.ones(2 * active.sum(), bool),
                      (np.r_[a[active], b[active]], np.r_[b[active], a[active]])),
                     shape=(n, n)).tocsr()
    count, labels = connected_components(adj, directed=False)
    ht, hb = np.zeros(count, bool), np.zeros(count, bool)
    ht[labels[top[gt > 0]]] = True
    hb[labels[bottom[gb > 0]]] = True
    useful = (ht & hb)[labels]
    phi[ht[labels] & ~hb[labels]] = 1.
    phi[hb[labels] & ~ht[labels]] = 0.
    if useful.any():
        ix = np.full(n, -1, int)
        ix[useful] = np.arange(useful.sum())
        edge = useful[a] & useful[b] & active
        aa, bb, gg = ix[a[edge]], ix[b[edge]], g[edge]
        ti, bi = useful[top] & (gt > 0), useful[bottom] & (gb > 0)
        tt, bt = ix[top[ti]], ix[bottom[bi]]
        gtt, gbt = gt[ti], gb[bi]
        nn = useful.sum()
        diagonal = np.bincount(np.r_[aa, bb, tt, bt], np.r_[gg, gg, gtt, gbt], minlength=nn)
        source = np.bincount(tt, gtt, minlength=nn)
        scale = diagonal.max()
        lap = coo_matrix((np.r_[diagonal, -gg, -gg] / scale,
                          (np.r_[np.arange(nn), aa, bb], np.r_[np.arange(nn), bb, aa])),
                         shape=(nn, nn)).tocsc()
        potential = spsolve(lap, source / scale, permc_spec='MMD_AT_PLUS_A')
        if not np.isfinite(potential).all():
            raise FloatingPointError('Nonfinite fractional-network potential')
        phi[useful] = potential
        if fields:
            numerical_residual = float(np.max(np.abs(lap @ potential - source / scale)) * scale)
    dv2 = (phi[a] - phi[b]) ** 2
    vt2, vb2 = (1. - phi[top]) ** 2, phi[bottom] ** 2
    pe, pt, pb = g * dv2, gt * vt2, gb * vb2
    material = float(pe.sum() + pt.sum() + pb.sum())
    result = dict(conductance_S=material + leakage_S, material_conductance_S=material)
    if fields:
        internal = channels @ dv2
        cp = internal + ct @ vt2 + cb @ vb2
        coords, shape = t['coords'], t['shape']
        pmap = np.zeros(shape)
        pmap[tuple(coords.T)] = np.bincount(
            np.r_[a, b, top, bottom], np.r_[.5 * pe, .5 * pe, pt, pb], minlength=n)
        phimap = np.broadcast_to(1. - (np.arange(shape[0])[:, None] + .5) / shape[0], shape).copy()
        phimap[tuple(coords.T)] = phi
        result.update(
            unit_phi=phimap, unit_power=pmap,
            top_current_A=float(gt @ (1. - phi[top])) + leakage_S,
            bottom_current_A=float(gb @ phi[bottom]) + leakage_S,
            max_kcl_residual_A=numerical_residual, nodes=n,
            component_count=int(count), conducting_nodes=int(useful.sum()),
            metal_power_W=float(cp[0]), tunnel_power_W=float(cp[1]),
            defect_power_W=float(cp[2]),
            metal_edge_power_W=float(internal[0]), tunnel_edge_power_W=float(internal[1]),
            defect_edge_power_W=float(internal[2]),
            nonmetal_power_W=float(cp[1:].sum()),
            nonmetal_power_fraction=float(cp[1:].sum() / material) if material else 1.,
            contact_power_W=float(pt.sum() + pb.sum()))
    return result
