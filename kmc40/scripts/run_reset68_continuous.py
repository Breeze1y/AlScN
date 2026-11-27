"""Continue the shared-residual KMC chain from the original INIT_RESET endpoint.

Only archived spatial catalogues with identical Ag-ID/VN arrays, cumulative
counts, seed and input RNG may be reused. Catalogue event generation is
independent of electronic state (erosion field at q=1). SET/RESET clocks and
readout always start from the live electronic checkpoint.
"""
import argparse
import csv
import shutil
from pathlib import Path

import numpy as np

from field_growth40 import FieldGrowthParams, checkpoint_at_event as growth_checkpoint
from rate_reset40 import checkpoint_at_event as erosion_checkpoint
from mixed_network40 import NetworkParams
import run_full_chain as full
from run_ag_constrained_cycles import run_reset, run_set
from run_local_retreat import Params
from run_ten_cycles import source_info

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "results/cycle_reset_fit40/mid_clock_set"
OUT = ROOT / "replay/reset68/continuous"
ORIGINAL = ROOT / "replay/cycles100"
COLS = ("point","branch","voltage_V","current_A","endpoint_current_A","global_elapsed_units","read_G_S")


def archived_spatial_inputs_match(folder, checkpoint, phase, cfg):
    meta = full.js(folder/"restart.json")
    cp = full.js(checkpoint/"restart.json")
    if (meta["initial_rng_state"] != cp["rng_state"] or meta["seed"] != cp["seed"]
            or any(meta["initial_source"][k] != cp[k] for k in ("injected","returned"))):
        return False
    if phase == "SET":
        if meta["params"] != cfg["growth"]:
            return False
    else:
        if (meta["field_solver"].get("activity") != 1.
                or meta["field_solver"].get("params") != cfg["network"]
                or meta["params"]["background_weight"] != cfg["negative"]["background_oxidation_weight"]):
            return False
    with np.load(folder/"sequence.npz") as seq, np.load(checkpoint/"restart.npz") as old:
        return (np.array_equal(seq["metal_id"][0],old["metal_id"])
                and np.array_equal(seq["VN"],old["VN"]))


def prepare_catalogue(archived, dest, checkpoint, phase, cfg):
    if (dest/"sequence.npz").exists():
        original = full.js(archived/"restart.json")
        existing = full.js(dest/"restart.json")
        return ("reused_exact_spatial_inputs" if existing["initial_rng_state"]==original["initial_rng_state"]
                and archived_spatial_inputs_match(archived,checkpoint,phase,cfg) else "existing")
    if not archived_spatial_inputs_match(archived,checkpoint,phase,cfg):
        return "generated"
    shutil.copytree(archived,dest,dirs_exist_ok=True)
    meta = full.js(dest/"restart.json")
    cp = full.js(checkpoint/"restart.json")
    provenance = meta["initial_source"]
    provenance["path"] = str(checkpoint.relative_to(ROOT))
    hash_key = "restart_sha256" if phase=="SET" else "sha256"
    provenance[hash_key] = full.sha(checkpoint/"restart.npz")
    if phase=="SET":
        meta["electronic_state"] = cp["electronic_state"]
    else:
        provenance["electronic_state"] = cp["electronic_state"]
    full.dump(dest/"restart.json",meta)
    return "reused_exact_spatial_inputs"


def checkpoint_writer(archived, phase, reuse):
    normal = growth_checkpoint if phase=="SET" else erosion_checkpoint
    def write(folder, index, dest, observation, *electronics):
        read = full.js(archived/"read_state.json") if reuse else {}
        if read.get("event_index")!=index:
            return normal(folder,index,dest,observation,*electronics)
        # Same spatial event prefix gives the original material IDs and RNG at
        # this event. Its surrounding waveform and electronic state can differ.
        with np.load(folder/"sequence.npz") as seq, np.load(archived/"restart.npz") as cp:
            np.testing.assert_array_equal(seq["metal_id"][index],cp["metal_id"])
            np.testing.assert_array_equal(seq["VN"],cp["VN"])
            clock=seq["cumulative_hazard_clock"]
            assert read["last_event_exposure"]==clock[index] and observation>=clock[index]
            assert index+1>=len(clock) or observation<clock[index+1]
        shutil.copytree(archived,dest,dirs_exist_ok=True)
        data=full.js(dest/"restart.json")
        data["initial_source"]=full.js(folder/"restart.json")["initial_source"]
        if electronics:
            data["electronic_state"]=electronics[0]
        full.dump(dest/"restart.json",data)
        read.update(source_sequence=str(folder.resolve()),observation_exposure=observation,
                    original_spatial_prefix=str(archived.relative_to(ROOT)))
        full.dump(dest/"read_state.json",read)
    return write


def export_iv(total):
    folder = OUT/"iv"
    folder.mkdir(exist_ok=True)
    temp = folder/"iv_all.tmp"
    with temp.open("w",encoding="utf-8-sig",newline="") as h:
        writer = csv.writer(h)
        writer.writerow(("cycle","phase",*COLS))
        for n in range(1,total+1):
            for phase in ("SET","RESET"):
                with (OUT/f"cycle_{n:02d}"/phase/"iv.csv").open(encoding="utf-8-sig",newline="") as src:
                    rows = list(csv.DictReader(src))
                if len(rows)!=202:
                    raise ValueError(f"Incomplete saved phase {n} {phase}")
                writer.writerows((n, phase, *(row[key] for key in COLS)) for row in rows)
    temp.replace(folder/"iv_all.csv")


def summarize_chain(total):
    rows = []
    parent = REF/"INIT_RESET/checkpoint"
    for n in range(1,total+1):
        original=(REF if n<=5 else ORIGINAL)/f"cycle_{n:02d}"
        r={"cycle":n}
        for phase in ("SET","RESET"):
            d = OUT/f"cycle_{n:02d}"/phase
            cat = d.parent/"catalogues"/phase
            before, s = full.js(parent/"restart.json"), full.js(d/"summary.json")
            checkpoint=d/"checkpoint"
            after = full.js(checkpoint/"restart.json")
            assert s["initial_electronic_state"]==before["electronic_state"]
            assert s["final_electronic_state"]==after["electronic_state"]
            assert full.js(cat/"restart.json")["initial_rng_state"]==before["rng_state"]
            with np.load(parent/"restart.npz") as cp, np.load(cat/"sequence.npz") as seq:
                np.testing.assert_array_equal(seq["metal_id"][0],cp["metal_id"])
                np.testing.assert_array_equal(seq["VN"],cp["VN"])
                for field in ("injected","returned"):
                    assert seq["initial_"+field]==before[field]
            with np.load(d/"states.npz") as st, np.load(checkpoint/"restart.npz") as cp:
                np.testing.assert_array_equal(cp["metal_id"],st["metal_id"][-1])
                np.testing.assert_array_equal(cp["VN"],st["VN"])
                np.testing.assert_array_equal(st["Ag_total"]+st["returned"],st["injected"])
                assert len(st["point"])==202
                r[f"{phase}_start_elapsed"]=s["initial_electronic_state"]["elapsed_units"]
                r[f"{phase}_end_elapsed"]=s["final_electronic_state"]["elapsed_units"]
                r[f"{phase}_R_ohm"]=1/st["read_G_S"][-1]
                r[f"{phase}_Ag"]=st["Ag_total"][-1]
            if phase=="RESET":
                assert not s["pending_microsteps"] and s["retreat_mode"]=="shared_residual"
                r["shared_sites_at_break"]=(s.get("residual_cohort") or {}).get("sites",0)
            r[f"{phase}_catalogue"]=s["catalogue_reuse"]
            r[f"{phase}_event"]=s["final_event"]
            r[f"{phase}_waveform_reused"]=bool(s.get("phase_reuse"))
            with np.load(original/phase/"checkpoint/restart.npz") as baseline, np.load(checkpoint/"restart.npz") as cp:
                r[f"{phase}_geometry_equals_original"]=np.array_equal(cp["metal_id"],baseline["metal_id"])
            r[f"{phase}_RNG_equals_original"]=after["rng_state"]==full.js(original/phase/"checkpoint/restart.json")["rng_state"]
            r[f"{phase}_end_electronics_equals_original"]=(
                s["final_electronic_state"]==full.js(original/phase/"summary.json")["final_electronic_state"])
            parent=checkpoint
        rows.append(r)
    with (OUT/"chain_summary.csv").open("w",encoding="utf-8-sig",newline="") as h:
        writer=csv.DictWriter(h,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    return rows


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--cycles",type=int,default=100,choices=range(1,101),metavar="1..100")
    total = cli.parse_args().cycles
    cfg = full.js(REF/"parameters.json")
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/"progress.json").exists() and full.js(OUT/"progress.json").get("completed",0)>total:
        cli.error("The saved trajectory already contains more cycles than requested")
    initial = REF/"INIT_RESET/checkpoint"
    configuration = dict(model="shared_residual",initial_checkpoint=str(initial.relative_to(ROOT)),
                         integration_substeps=cfg["substeps"],positive=cfg["positive"],negative=cfg["negative"],
                         growth=cfg["growth"],network=cfg["network"],local_retreat=cfg["local_retreat"],
                         logq_grid=cfg["logq_grid"])
    params = OUT/"parameters.json"
    if params.exists() and full.js(params)!=configuration:
        raise ValueError("The saved continuous trajectory has different parameters")
    full.dump(params,configuration)
    checkpoint = initial
    growth, net, retreat = FieldGrowthParams(**cfg["growth"]), NetworkParams(**cfg["network"]), Params(**cfg["local_retreat"])
    grid = np.asarray(cfg["logq_grid"])
    for n in range(1,total+1):
        current = OUT/f"cycle_{n:02d}"
        original = (REF if n<=5 else ORIGINAL)/f"cycle_{n:02d}"
        for phase in ("SET","RESET"):
            dest = current/phase
            folder = current/"catalogues"/phase
            done = dest/"checkpoint/restart.npz"
            info = dest/"summary.json"
            if done.exists() and info.exists() and full.js(info).get("phase")==phase:
                saved = full.js(info)
                if (saved["source_checkpoint_sha256"] != full.sha(checkpoint/"restart.npz")
                        or saved["inherited_rng_state"]!=full.js(checkpoint/"restart.json")["rng_state"]
                        or saved["initial_electronic_state"]!=full.js(checkpoint/"restart.json")["electronic_state"]):
                    raise ValueError(f"Checkpoint mismatch on resume at {n} {phase}")
                checkpoint = done.parent
                continue
            # save_phase writes an intermediate summary before the final checkpoint.
            if info.exists():
                info.replace(dest/"intermediate_summary.json")
            full.dump(OUT/"progress.json",dict(status="running",cycle=n,phase=phase,completed=n-1))
            reused = prepare_catalogue(original/"catalogues"/phase,folder,checkpoint,phase,cfg)
            write_cp = checkpoint_writer(original/phase/"checkpoint",phase,
                                         reused=="reused_exact_spatial_inputs")
            parent = checkpoint
            cp_meta = full.js(parent/"restart.json")
            if phase=="SET":
                template = full.js(original/phase/"summary.json")
                same_electrons = cp_meta["electronic_state"]==template["initial_electronic_state"]
                if reused=="reused_exact_spatial_inputs" and same_electrons:
                    shutil.copytree(original/phase,dest,dirs_exist_ok=True)
                    checkpoint = dest/"checkpoint"
                    cs = full.js(checkpoint/"restart.json")
                    cs["initial_source"] = full.js(folder/"restart.json")["initial_source"]
                    full.dump(checkpoint/"restart.json",cs)
                    selected=full.js(checkpoint/"read_state.json")
                    selected.update(source_sequence=str(folder.resolve()),
                                    original_spatial_prefix=str((original/phase/"checkpoint").relative_to(ROOT)))
                    full.dump(checkpoint/"read_state.json",selected)
                    saved = full.js(info)
                    saved.update(source_info(parent,folder,n,phase),source_chain=str(OUT.relative_to(ROOT)))
                    saved["phase_reuse"]="identical starting geometry, RNG, Ag balance, electronic state and SET parameters"
                    full.dump(info,saved)
                else:
                    checkpoint=run_set(OUT,n,parent,cfg["positive"],net,grid,cfg["substeps"],growth,
                                       checkpoint_writer=write_cp)
            else:
                checkpoint=run_reset(OUT,n,parent,cfg["negative"],net,grid,retreat,cfg["substeps"],
                                     wait_for_return_route=True,share_residual_retreat=True,checkpoint_writer=write_cp)
            saved = full.js(info)
            assert saved["initial_electronic_state"]==cp_meta["electronic_state"]
            assert saved["inherited_rng_state"]==cp_meta["rng_state"]
            saved["catalogue_reuse"] = reused
            full.dump(info,saved)
        print(f"Continuous cycle {n}/{total} finished",flush=True)
    export_iv(total)
    summarize_chain(total)
    full.dump(OUT/"progress.json",dict(status="complete",completed=total,iv_rows=404*total))


if __name__=="__main__":
    main()
