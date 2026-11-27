"""Update the three static IV deliverables using the full shared-flux trajectory."""
import csv
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "replay/reset68/continuous/iv/iv_all.csv"
OUT = ROOT / "replay/cycles100/report"
PHASES = ("SET", "RESET")
COLORS = {"SET": "#126775", "RESET": "#ab4c27"}
BACKGROUND = "#73808b"
FOCUS = 68


def read_curves():
    curves = {}
    with SOURCE.open(encoding="utf-8-sig", newline="") as h:
        reader = csv.DictReader(h)
        for cycle in range(1, 101):
            for phase in PHASES:
                rows = [next(reader) for _ in range(202)]
                for idx, r in enumerate(rows):
                    if (int(r["cycle"]), r["phase"], int(r["point"]), r["branch"]) != (
                            cycle, phase, idx, "outward" if idx <= 100 else "return"):
                        raise ValueError(f"Unexpected scan order in cycle {cycle}, {phase}, point {idx}")
                v = np.array([float(r["voltage_V"]) for r in rows])
                amps = np.array([float(r["current_A"]) for r in rows])
                signed_grid = .015*np.r_[np.arange(101), np.arange(100, -1, -1)]
                np.testing.assert_allclose(v, signed_grid if phase=="SET" else -signed_grid, atol=1e-9)
                if not np.all(np.isfinite(amps)) or np.any(amps[1:-1] == 0):
                    raise ValueError(f"Invalid biased current in cycle {cycle}, {phase}")
                curves[cycle,phase] = v, amps, 1/float(rows[-1]["read_G_S"])
        if next(reader,None) is not None:
            raise ValueError("Unexpected records after cycle 100")
    return curves


def draw(curves, output=OUT):
    mpl.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial","DejaVu Sans"],
                         "font.size":8.5, "axes.titlesize":10, "pdf.fonttype":42, "svg.fonttype":"none",
                         "axes.spines.top":False, "axes.spines.right":False, "axes.linewidth":.65,
                         "axes.labelcolor":"#25323b", "text.color":"#25323b", "legend.frameon":False})
    fig, (iv, rax) = plt.subplots(1,2, figsize=(8.35,4.05), gridspec_kw={"width_ratios":(1.4,1)})
    fig.subplots_adjust(left=.10,right=.985,bottom=.28,top=.77,wspace=.30)
    fig.text(.10,.94,"100 continuous KMC cycles", fontsize=15, fontweight="bold", va="center")
    fig.text(.10,.87,"Shared residual flux during RESET", fontsize=9, va="center")
    for ax in (iv,rax):
        ax.set_yscale("log")
        ax.grid(axis="y",which="major",lw=.55,color="#dee3e6")
        ax.set_axisbelow(True)
        ax.tick_params(which="major",width=.65,length=3)
    # Each trace has two separate paths. Exclude precisely the two 0 V/0 A
    # endpoints on this logarithmic plot; the source CSV retains both.
    branches = (("outward",slice(1,101),"-"),("return",slice(101,201),"--"))
    for cycle in [n for n in range(1,101) if n!=FOCUS]+[FOCUS]:
        focus = cycle==FOCUS
        for phase in PHASES:
            volts, amps, _ = curves[cycle,phase]
            for name, sub, style in branches:
                line, = iv.plot(volts[sub], np.abs(amps[sub])*1e6, style,
                                color=COLORS[phase] if focus else BACKGROUND,
                                lw=1.65 if focus else .52, alpha=1. if focus else .23,
                                zorder=4 if focus else 2)
                line.set_gid(f"cycle-{cycle:03d}-{phase}-{name}")
    iv.set(xlim=(-1.55,1.55), xticks=np.arange(-1.5,1.6,.5),
           xlabel="Programmed voltage (V)", ylabel="Mean |current| (µA)",
           title="SET and RESET sweeps")
    cycle_ids = np.arange(1,101)
    for phase in PHASES:
        ohms = [curves[n,phase][2] for n in cycle_ids]
        line, = rax.plot(cycle_ids,ohms,color=COLORS[phase],lw=1.15)
        line.set_gid(f"endpoint-resistance-{phase}")
        rax.plot(FOCUS,ohms[FOCUS-1],"o",ms=5.3,markerfacecolor="white",
                 markeredgecolor=COLORS[phase],markeredgewidth=1.5)
    rax.set(xlim=(1,100),xticks=(1,25,50,68,100),xlabel="Cycle",
            ylabel="Read resistance (Ω)",title="Resistance at each phase endpoint")
    handles = [Line2D([],[],color=BACKGROUND,lw=1.5,label="Other 99 cycles"),
               Line2D([],[],color=COLORS["SET"],lw=1.8,label="Cycle 68 SET"),
               Line2D([],[],color=COLORS["RESET"],lw=1.8,label="Cycle 68 RESET"),
               Line2D([],[],color="#35444c",lw=1.1,label="Outward"),
               Line2D([],[],color="#35444c",lw=1.1,ls="--",label="Return")]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.53,.063),
               ncols=5,fontsize=8.4,handlelength=2.1,columnspacing=1.3)
    set_r, reset_r = (curves[FOCUS,phase][2] for phase in PHASES)
    fig.text(.10,.029,f"The cycle 68 endpoints (open circles) are {set_r:.1f} Ω after SET and {reset_r/1e6:.3f} MΩ after RESET.",
             fontsize=8, color="#4b565e")
    source = "Simulation source: replay/reset68/continuous/iv/iv_all.csv. 100 continuous SET–RESET cycles."
    output.mkdir(parents=True,exist_ok=True)
    fig.savefig(output/"iv_100cycles.png", dpi=300, metadata={"Description":source})
    fig.savefig(output/"iv_100cycles.pdf", metadata={"Title":"100 continuous KMC cycles","Subject":source})
    fig.savefig(output/"iv_100cycles.svg", metadata={"Title":"100 continuous KMC cycles","Description":source})
    plt.close(fig)


if __name__ == "__main__":
    draw(read_curves())
