"""Fit local sub-grid retreat of the narrowest neck in the existing AlScN40 KMC L state."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.ndimage import label
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "baseline/diff/fit_baselines/AlScN_40nm_file074.csv"
GEOMETRY = ROOT / "results/refined_lattice_20_40/spatial/AlScN40/states.npz"
OUTPUT = ROOT / "results/reset40_gradual/spatial"
CELL = 0.5
BETA = 10.25   # 1/nm, rectangular 1 eV barrier, free-electron mass
AG_LAYER = .408626 / np.sqrt(3)
GAP_STEP = 2 * AG_LAYER    # high-field retreat scale, two FCC(111) spacings
PARAMETERS = ["onset_V", "dose_scale_V", "feedback_length_nm", "depin_scale_V",
              "depin_slope_V", "parallel_floor_uS"]


def source_rows():
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    resets = [r for r in rows if r["phase"] == "RESET"]
    v = np.array([float(r["voltage_V"]) for r in resets])
    i = np.array([float(r["current_A"]) for r in resets])
    read = [r for r in rows if r["phase"] == "SET" and r["branch"] == "return"
            and .014 < float(r["voltage_V"]) < .061]
    g = np.median([float(r["current_A"]) / float(r["voltage_V"]) for r in read])
    return v, i, resets, g


def find_neck(initial):
    metal = initial[1:-1] == 2
    regions, _ = label(metal)
    spanning = np.intersect1d(regions[0][regions[0] > 0], regions[-1][regions[-1] > 0])
    assert len(spanning) == 1
    cluster = np.isin(regions, spanning)
    row = np.argmin(cluster.sum(axis=1)) + 1
    contacts = np.flatnonzero(cluster[row-1])
    local_coord = np.array([np.count_nonzero(initial[row-2:row+3,x-2:x+3] == 2)
                            for x in contacts])
    susceptibility = local_coord.mean() / local_coord
    return int(row), contacts, susceptibility, cluster


def interface(p, voltage, step, susceptibility, high_field):
    onset, scale, feedback, depin, slope, _ = p
    magnitude = np.maximum(-voltage[:, None], 0)
    # Equal dwell per acquisition tick (absolute time unknown). step fixes the dose unit.
    exposure = step * np.maximum(magnitude-onset, 0) * susceptibility / scale**2
    slow = feedback * np.log1p(np.cumsum(exposure, axis=0))
    # Difference of exponentials enforces zero reverse reaction at exactly zero voltage.
    hazard = np.exp(np.clip((magnitude-depin)/slope, -700, 60))
    hazard *= -np.expm1(-magnitude/slope)
    hazard = step/slope * hazard * susceptibility
    fast = GAP_STEP * (-np.expm1(-np.cumsum(hazard, axis=0))) if high_field else np.zeros_like(slow)
    return slow+fast, slow, fast


def conductance(gap, floor_uS, g0):
    return floor_uS * 1e-6 + (g0-floor_uS*1e-6) * np.exp(-BETA*gap).mean(axis=1)


def calibrate(v, i, mask, step, susceptibility, g0, high_field):
    lo = np.array([0, .02, .025, 1.12, .008, .03])
    hi = np.array([.17, .7, .7, 1.70, .22, 3.5])
    active = np.arange(6) if high_field else np.array([0,1,2,5])
    def unpack(x):
        p = np.array([.05, .12, .16, 1.4, .05, .85])
        p[active] = x
        return p
    def residual(x):
        p = unpack(x)
        gap, _, _ = interface(p, v, step, susceptibility, high_field)
        pred = v * conductance(gap, p[-1], g0)
        return np.r_[np.log10(np.abs(pred[mask]/i[mask])), 2*(pred[mask]-i[mask])/.001]
    fits = []
    for guess in [.02, .055, .08, .11]:
        p0 = np.array([guess,.12,.16,1.42,.04,.85])
        fits.append(least_squares(residual, p0[active], bounds=(lo[active],hi[active])))
    return unpack(min(fits, key=lambda r:r.cost).x)


def fractional_geometry(initial, row, contacts, gap):
    # Strip metal from the Pt-facing end of this Ag neck into its Ag-side reservoir.
    # Include the wider rim in layers above the initial narrowest row.
    face = row*CELL
    columns = np.arange(initial.shape[1])
    lane = np.argmin(np.abs(columns[:,None]-contacts), axis=1)
    distance = gap[:,lane]
    top = (np.arange(len(initial)) - 1)*CELL
    above = (top+CELL <= face) & (np.arange(len(initial)) > 0)
    # A partially penetrated cell loses only the fraction between the plane and its bottom.
    full_distance_from_face = (face-top-CELL)[None,:,None]
    removed = np.clip((distance[:,None,:] - full_distance_from_face)/CELL,0,1) * above[None,:,None]
    return (initial == 2)[None,:,:] * (1 - removed)


def score(fit, measure, masks):
    out = {}
    for name, m in masks.items():
        r = np.log10(np.abs(fit[m]/measure[m]))
        out[name] = dict(n=int(m.sum()), log10_current_RMSE=float(np.mean(r*r)**.5),
                         rms_uA=float(1e6*np.mean((fit[m]-measure[m])**2)**.5))
    return out


def plots(v, i, fit_i, gaps, full, partial, times):
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig, (ax, logax, gx)=plt.subplots(1,3,figsize=(13.8,3.9),layout="constrained")
    for sl, name, co in [(slice(0,101),"Outward","#b74721"),(slice(101,202),"Return","#225e89")]:
        for panel, scale, unit in [(ax,1000,"mA"),(logax,1e6,"µA")]:
            valid = v[sl] != 0 if panel is logax else np.ones(v[sl].size,dtype=bool)
            panel.scatter(-v[sl][valid],np.abs(i[sl][valid])*scale,s=10,facecolors="none",edgecolors=co,alpha=.65)
            panel.plot(-v[sl][valid],np.abs(fit_i[sl][valid])*scale,color=co,label=name,linewidth=1.8)
            panel.set(xlabel="Reverse bias magnitude (V)",ylabel=f"|I| ({unit})")
        gx.plot(-v[sl],gaps[sl].mean(axis=1),color=co,label=name,linewidth=1.8)
    ax.set(title="A   Current peak and decay")
    logax.set(yscale="log",ylim=(.025,1500),title="B   Return branch and tail")
    gx.set(xlabel="Reverse bias magnitude (V)",ylabel="Mean local retreat (nm)",title="C   State carried through reversal")
    ax.legend(frameon=False);logax.text(.06,.49,"Circles  measurement\nLines     model",transform=logax.transAxes,fontsize=9)
    fig.savefig(OUTPUT/"reset_iv.png", dpi=200);plt.close(fig)
    fig = plt.figure(figsize=(13,5.4),layout="constrained")
    gs=fig.add_gridspec(2,4,width_ratios=[.9,1.3,1.3,1.3]); ax=fig.add_subplot(gs[:,0])
    ax.imshow(full[1:-1]==2,extent=[0,16,40,0],cmap="Blues",vmin=0,vmax=1,interpolation="nearest")
    ax.add_patch(Rectangle((5,3.25),7,3.75,fill=False,edgecolor="#d45223",lw=1.8))
    ax.set(xlabel="x (nm)",ylabel="Depth from Ag (nm)",title="Initial Ag network")
    for n,t in enumerate(times):
        sub=fig.add_subplot(gs[n//3,1+n%3])
        patches=[]
        for y,x in zip(*np.nonzero(partial[t,1:-1,:]>0)):
            fill=partial[t,y+1,x]
            if 2.75 <= y*CELL <= 8 and 5<= x*CELL<=12:
                # Pt-facing surfaces recede, so remaining metal occupies the upper side.
                patches.append(Rectangle((x*CELL,y*CELL), CELL, CELL*fill))
        sub.add_collection(PatchCollection(patches,facecolor="#204c69",edgecolor="#d6e2e7",linewidth=.18))
        sub.set(xlim=(5,12),ylim=(7,3.25),aspect="equal",xlabel="x (nm)",
                title=f"{'Return' if t>=101 else 'Out'}  {v[t]:.3f} V    gap {gaps[t].mean():.2f} nm")
        if n%3==0:sub.set_ylabel("Depth (nm)")
    fig.suptitle("AlScN 40 nm   |   Sub-grid Ag retreat at the narrowest neck")
    fig.savefig(OUTPUT/"neck_states.png",dpi=200);plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True,exist_ok=True)
    v,i,rows,g0=source_rows()
    initial=np.load(GEOMETRY)["L"]
    row, contacts, susceptibility, cluster=find_neck(initial)
    step=float(np.median(np.abs(np.diff(v[:101]))))
    all_points=v != 0
    train=all_points & (np.arange(len(v))%2==0)
    masks={"all_nonzero":all_points,"calibration_interleaved":train,
           "interleaved_withheld_same_curve":all_points & ~train,
           "outward":all_points&(np.arange(len(v))<101),"return":all_points&(np.arange(len(v))>=101)}
    summary={}
    for field in (False,True):
        p=calibrate(v,i,train,step,susceptibility,g0,high_field=field)
        gaps,slow,fast=interface(p,v,step,susceptibility,high_field=field)
        G=conductance(gaps,p[-1],g0); pred=v*G
        summary["with_high_field" if field else "continuous_only"] = dict(
            parameters={k:val for k,val in zip(PARAMETERS,p) if field or not k.startswith("depin")},
            scores=score(pred,i,masks))
    occupancy=fractional_geometry(initial,row,contacts,gaps)
    # Independent gap reading along the original two contiguous Ag columns.
    removed=(initial==2)[None,:,:]-occupancy
    rebuilt_gap=removed[:,:,contacts].sum(axis=1)*CELL
    assert np.allclose(rebuilt_gap,gaps)
    assert np.allclose(conductance(rebuilt_gap,p[-1],g0)*v,pred)
    assert np.allclose(occupancy[0],initial==2) and np.all(np.diff(occupancy,axis=0)<=1e-12)
    assert not np.any(removed[:,[0,-1],:])
    final_labels, _ = label(occupancy[-1,1:-1,:]>0)
    final_spanning = np.intersect1d(final_labels[0][final_labels[0]>0],
                                    final_labels[-1][final_labels[-1]>0])
    assert final_spanning.size == 0
    inventory = occupancy[:,1:-1,:].sum(axis=(1,2))
    atoms_per_cell = CELL**2 * 1.0 * 4/.408626**3
    times=[0,9,20,60,92,201]
    summary.update(dict(input=str(SOURCE),input_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        geometry_source=str(GEOMETRY),geometry_sha256=hashlib.sha256(GEOMETRY.read_bytes()).hexdigest(),
        geometry_state="L, independent prior 0.5 nm KMC post-positive-pulse surrogate",
        cell_nm=CELL, row_index_with_top_electrode=row, contact_columns=contacts.tolist(),
        susceptibility=susceptibility.tolist(), paired_SET_return_G0_S=g0,
        measured_reset_initial_G_S=float(np.median((i/v.clip(max=-1e-15))[1:5])),
        H1_read_G_S=G[-1], H1_read_R_ohm=1/G[-1], H1_mean_retreat_nm=gaps[-1].mean(),
        beta_per_nm=BETA, high_field_saturation_depth_nm=GAP_STEP,
        voltage_step_V=step, curve_length=len(v), calibration_points=int(train.sum()),
        measured_out_peak_V=v[np.argmax(abs(i[:101]))], fitted_out_peak_V=v[np.argmax(abs(pred[:101]))],
        measured_out_peak_mA=max(abs(i[:101]))*1e3, fitted_out_peak_mA=max(abs(pred[:101]))*1e3,
        max_reconstructed_gap_error_nm=float(np.max(np.abs(rebuilt_gap-gaps))),
        H1_metal_loss_equivalent_cells=removed[-1].sum(), snapshots=times,
        Ag_inventory_at_depth_1_nm=dict(L_metal_cells=float(inventory[0]),
            L_spanning_metal_cells=int(cluster.sum()), L_metal_atoms=inventory[0]*atoms_per_cell,
            L_spanning_atoms=cluster.sum()*atoms_per_cell, H_metal_cell_equivalent=inventory[-1],
            H_metal_atoms=inventory[-1]*atoms_per_cell, H_spanning_metal_cells=0,
            H_spanning_clusters=int(final_spanning.size), H_spanning_atoms=0,
            metal_depleted_atoms=(inventory[0]-inventory[-1])*atoms_per_cell)))
    (OUTPUT/"fit.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    with (OUTPUT/"reset_trajectory.csv").open("w",newline="",encoding="utf-8") as fp:
        out=csv.writer(fp);out.writerow(["index0","branch","V","I_exp_A","I_sim_A","G_read_S",
               "left_gap_nm","right_gap_nm","mean_continuous_nm","mean_high_field_nm",
               "removed_metal_cell_equiv","used_for_fit","source_line"])
        for j,r in enumerate(rows):
            out.writerow([j,r["branch"],v[j],i[j],pred[j],G[j],*gaps[j],slow[j].mean(),
                          fast[j].mean(),removed[j].sum(),int(train[j]),r["source_line"]])
    np.savez_compressed(OUTPUT/"neck_states.npz",voltage_V=v,initial_occupancy_codes=initial,
        metal_fraction=occupancy,contact_columns=contacts,contact_row=row,positive_pulse_L_state=initial,
        gap_nm=gaps,slow_dose=np.expm1(slow/p[2]),high_field_coverage=fast/GAP_STEP,
        read_G_S=G,parameters=p,G0_S=g0,susceptibility=susceptibility,step_V=step)
    plots(v,i,pred,gaps,initial,occupancy,times)
    print(json.dumps({k: summary[k] for k in ("with_high_field","H1_read_R_ohm","H1_mean_retreat_nm",
          "contact_columns","fitted_out_peak_mA","H1_metal_loss_equivalent_cells")},indent=2))

if __name__ == "__main__":
    main()

