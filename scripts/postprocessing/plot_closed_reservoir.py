#!/usr/bin/env python3
"""Closed-reservoir figure: emergent core temperature vs water temperature for the
defended and energy-limited variants. Shows that a generous sustained-scope ceiling
(4x Kleiber BMR) holds the core everywhere, while the model's own low FMR budget
cannot sustain the prescribed reservoir in cold water (core collapses)."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "results" / "closed_reservoir" / "closed_reservoir_summary.csv"
OUT = ROOT / "paper" / "figures" / "closed_reservoir.pdf"

rows = list(csv.DictReader(open(CSV)))
species = ["humpback", "bowhead", "minke"]
titles = {"humpback": "Humpback (thin blubber)", "bowhead": "Bowhead (thick blubber)",
          "minke": "Minke (calibration)"}

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
INK = "#1b2a4a"; RES = "#c44e34"; HOLD = "#2e7d5b"; GUARD = "#9aa0a6"
for ax, sp in zip(axes, species):
    r = [x for x in rows if x["species"] == sp]
    tw = sorted({float(x["T_water_C"]) for x in r})
    defended = {float(x["T_water_C"]): float(x["core_defended_C"]) for x in r}
    fmr = {float(x["T_water_C"]): float(x["core_emergent_C"])
           for x in r if x["ceiling_label"] == "fmr_nominal"}
    setp = float(r[0]["setpoint_C"])
    d = [defended[t] for t in tw]
    e = [fmr[t] for t in tw]
    # defended / kleiber x4 (identical - holds everywhere)
    ax.plot(tw, d, "-o", color=HOLD, ms=4, lw=2,
            label="Defended & $4\\times$BMR ceiling")
    # energy-limited under the tight FMR budget
    ax.plot(tw, e, "-s", color=RES, ms=4, lw=2,
            label="Energy-limited (baseline FMR)")
    ax.plot(tw, tw, ":", color=GUARD, lw=1, label="Water temperature")
    ax.axhline(28, ls="--", color=GUARD, lw=0.8)
    ax.set_title(titles[sp], fontsize=9)
    ax.set_xlabel("Water temperature (°C)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_xlim(-4, 32)
axes[0].set_ylabel("Core temperature (°C)", fontsize=8)
axes[0].set_ylim(-4, 40)
axes[-1].legend(fontsize=6.5, loc="lower right", framealpha=0.9)
fig.tight_layout()
fig.savefig(OUT, bbox_inches="tight")
print(f"wrote {OUT}")
# console summary
for sp in species:
    r = [x for x in rows if x["species"] == sp and x["ceiling_label"] == "fmr_nominal"]
    falls = [x for x in r if x["status"] == "core_falls"]
    print(f"{sp}: FMR-budget core falls in {len(falls)}/{len(r)} temps; "
          f"coldest emergent core = {min(float(x['core_emergent_C']) for x in r):.1f} C")
