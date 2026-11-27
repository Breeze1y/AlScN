"""Reconstruct the archived C68 actions and separate the frozen-state electrical effects."""
import csv
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
import local_fraction_network as network
from cycle_unified40 import CircuitParams, read_circuit
from mixed_network40 import NetworkParams

SOURCE = ROOT/"replay/cycles100/cycle_68/RESET"


def write_csv(name, rows):
    with (OUT/name).open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def main():
    summary = json.loads((SOURCE/"summary.json").read_text())
    circ = CircuitParams(**summary["clock"]["circuit"])
    params = NetworkParams(**summary["network"])
    a = np.load(SOURCE/"states.npz")
    with (SOURCE/"events.csv").open(encoding="utf-8-sig") as h:
        actions = list(csv.DictReader(h))
    n = int(a["subunits"])
    units = (a["initial_metal_id"] > 0).astype(int)*n
    topo = network.prepare(units > 0, a["VN"])

    def read(f, q, fields=False, topology=topo):
        mat = network.solve(topology, f, params, q, summary["parameters"]["contact_exponent"], fields)
        terminal = read_circuit(mat["conductance_S"], q, -1., circ)["read_G_S"]
        return terminal, mat

    events = []
    steps = []
    pointer = 0
    g0 = None
    for point, v in enumerate(a["voltage_V"]):
        prev_f = units/n
        while pointer < len(actions) and int(actions[pointer]["point"]) == point:
            e = actions[pointer]
            y, x = int(e["y"]), int(e["x"])
            q = float(e["q"])
            before = int(units[y,x])
            g1, _ = read(units/n, q)
            units[y,x] -= int(e["removed_units"])
            g2, _ = read(units/n, q)
            g_without_vn = None
            if a["VN"][y,x] and units[y,x] == 0:
                site = (topo["coords"][:,0] == y) & (topo["coords"][:,1] == x)
                muted = dict(topo, vn=topo["vn"] & ~site)
                g_without_vn, _ = read(units/n, q, topology=muted)
            assert units[y,x] == int(e["remaining_units"])
            events.append(dict(point=point, branch=a["branch"][point], voltage_V=v,
                               action=pointer+1, event=int(e["event_index"]), time=float(e["time"]),
                               y=y, x=x, kind=e["kind"], units_before=before,
                               units_after=units[y,x], q_at_event=q, route_nodes=e["route_nodes"],
                               connected_before=e["connected_before"], connected_after=e["connected_after"],
                               pre_G_uS=g1*1e6, post_G_uS=g2*1e6, loss_dec_at_fixed_q=np.log10(g1/g2),
                               VN_at_event_site=bool(a["VN"][y,x]),
                               post_G_without_revealed_VN_uS=1e6*g_without_vn if g_without_vn is not None else None))
            pointer += 1
        np.testing.assert_array_equal(units, a["Ag_units"][point])
        q = a["activity"][point]
        f = units/n
        g, material = read(f,q,True)
        if point:
            gpq, _ = read(prev_f, q)
            gfp, _ = read(f, a["activity"][point-1])
            # Log loss for changed geometry at both endpoint activities.
            # Averaging the two orders splits the nonlinear interaction.
            geo = .5 * (np.log10(g0/gfp) + np.log10(gpq/g))
            electronic = .5 * (np.log10(g0/gpq) + np.log10(gfp/g))
            loss = np.log10(g0/g)
            assert abs(loss-geo-electronic) < 1.e-12
        else:
            geo = electronic = loss = 0.
        g0 = g
        total = material["material_conductance_S"]
        steps.append(dict(point=point, branch=a["branch"][point], voltage_V=v,
                          archived_dwell_current_uA=1e6*a["current_A"][point],
                          archived_endpoint_G_uS=1e6*a["read_G_S"][point],
                          exact_endpoint_G_uS=1e6*g, q=q, Ag_in_film=a["Ag_total"][point],
                          metal_connected=a["metal_connected"][point], microactions=a["microevents"][point],
                          geometry_loss_dec=geo, activity_loss_dec=electronic, total_loss_dec=loss,
                          metal_power_fraction=material["metal_power_W"]/total,
                          tunnel_power_fraction=material["tunnel_power_W"]/total,
                          defect_power_fraction=material["defect_power_W"]/total))
    assert pointer == len(actions)
    write_csv("actual_event_effects.csv", events)
    write_csv("endpoint_attribution.csv", steps)


if __name__ == "__main__":
    main()
