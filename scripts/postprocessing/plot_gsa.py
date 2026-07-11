#!/usr/bin/env python3
"""Global sensitivity analysis figure for the manuscript.

Panel (a): Morris elementary-effects screening -- mean absolute effect mu* of each
of the 8 uncertain parameters on log(Q_surface) at Tw = -2 C, for the three species
spanning the strategy space (humpback/fin rorquals, bowhead balaenid).
Panel (b): LHS ranking-stability -- mean mass-specific surface cost per species with
the fraction of the 2000-sample ensemble in which the balaenid stays below both
rorquals annotated.

Reads results/gsa/{gsa_morris.csv, gsa_ranking_stability.json}.
"""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
GSA = ROOT / "results" / "gsa"
OUT = ROOT / "paper" / "figures" / "gsa.pdf"

# pretty parameter labels
PLABEL = {
    "k_blub": r"$k_{\mathrm{blub}}$", "BCI": "BCI", "h": r"$h$",
    "omega_core": r"$\omega_{\mathrm{core}}$", "cche": "CCHE",
    "ava": "AVA", "Ttn": r"$T_{\mathrm{tn}}$", "fmr": "FMR",
}
SP_LABEL = {"humpback": "Humpback", "fin_whale": "Fin", "bowhead": "Bowhead"}
SP_COLOR = {"humpback": "#c44e34", "fin_whale": "#d9a441", "bowhead": "#1b6f8c"}
SPECIES = ["humpback", "fin_whale", "bowhead"]

morris = list(csv.DictReader(open(GSA / "gsa_morris.csv")))
stab = json.load(open(GSA / "gsa_ranking_stability.json"))

# --- order parameters by mean mu* across species (most influential at top)
params = list(dict.fromkeys(r["param"] for r in morris))
mu_by = {(r["species"], r["param"]): float(r["mu_star"]) for r in morris}
order = sorted(params, key=lambda p: -np.mean([mu_by.get((s, p), 0) for s in SPECIES]))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 3.8),
                               gridspec_kw={"width_ratios": [1.55, 1]})

# ---- Panel (a): Morris mu* grouped bars
y = np.arange(len(order))
bw = 0.26
for k, sp in enumerate(SPECIES):
    vals = [mu_by.get((sp, p), 0.0) for p in order]
    axA.barh(y + (k - 1) * bw, vals, height=bw, color=SP_COLOR[sp],
             label=SP_LABEL[sp], edgecolor="white", lw=0.4)
axA.set_yticks(y)
axA.set_yticklabels([PLABEL.get(p, p) for p in order], fontsize=9)
axA.invert_yaxis()
axA.set_xlabel(r"Morris $\mu^{*}$  (mean $|$effect$|$ on $\log Q_{\mathrm{surface}}$)",
               fontsize=9)
axA.tick_params(labelsize=8)
axA.legend(fontsize=8, loc="lower right", framealpha=0.9)
axA.grid(axis="x", alpha=0.25, lw=0.5)
axA.set_title("(a) Parameter influence (Morris screening)", fontsize=9, loc="left")

# ---- Panel (b): ranking stability
st = stab["stability"]
msc = st["mean_msc_W_per_kg"]
sp_present = [s for s in SPECIES if s in msc]
vals = [msc[s] for s in sp_present]
xb = np.arange(len(sp_present))
bars = axB.bar(xb, vals, color=[SP_COLOR[s] for s in sp_present],
               edgecolor="white", lw=0.6, width=0.62)
axB.set_xticks(xb)
axB.set_xticklabels([SP_LABEL[s] for s in sp_present], fontsize=9)
axB.set_ylabel(r"Mean mass-specific surface cost ($\mathrm{W\,kg^{-1}}$)", fontsize=9)
axB.tick_params(labelsize=8)
axB.grid(axis="y", alpha=0.25, lw=0.5)
for b, v in zip(bars, vals):
    axB.text(b.get_x() + b.get_width() / 2, v, f"{v:.3f}",
             ha="center", va="bottom", fontsize=7.5)
pct = 100 * st["frac_bowhead_below_both"]
axB.set_title("(b) Family ranking stability", fontsize=9, loc="left")
axB.annotate(f"Bowhead $<$ both rorquals\nin {pct:.0f}% of "
             f"{st['n_samples']}-sample LHS ensemble",
             xy=(0.5, 0.97), xycoords="axes fraction", ha="center", va="top",
             fontsize=8, bbox=dict(boxstyle="round,pad=0.35", fc="#f2f2f2",
                                   ec="#9aa0a6", lw=0.6))
axB.set_ylim(0, max(vals) * 1.28)

fig.tight_layout()
fig.savefig(OUT, bbox_inches="tight")
print(f"wrote {OUT}")
print(f"ranking stability: bowhead<both = {pct:.1f}%  "
      f"(< humpback {100*st['frac_bowhead_below_humpback']:.1f}%, "
      f"< fin {100*st['frac_bowhead_below_fin']:.1f}%)")
print("Morris mu* order:", ", ".join(order[:4]))
