#!/usr/bin/env python3
"""Generate the GSA appendix table (Sobol total-effect + Morris mu*) from the
results/gsa CSVs. Emits paper/generated/gsa_table.tex as a standalone-safe table
(plain-text citations, no \\citet). Mirrors the other build_*_tables.py generators.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GSA = ROOT / "results" / "gsa"
OUT = ROOT / "paper" / "generated" / "gsa_table.tex"

PLABEL = {
    "k_blub": r"$k_{\mathrm{blub}}$", "BCI": "BCI", "h": r"$h$",
    "omega_core": r"$\omega_{\mathrm{core}}$", "cche": "CCHE",
    "ava": "AVA", "Ttn": r"$T_{\mathrm{tn}}$", "fmr": "FMR",
}
SP_LABEL = {"humpback": "Humpback", "fin_whale": "Fin", "bowhead": "Bowhead"}
SPECIES = ["humpback", "fin_whale", "bowhead"]

sobol = list(csv.DictReader(open(GSA / "gsa_sobol.csv")))
morris = list(csv.DictReader(open(GSA / "gsa_morris.csv")))
stab = json.load(open(GSA / "gsa_ranking_stability.json"))
params = list(dict.fromkeys(r["param"] for r in sobol))

ST = {(r["species"], r["param"]): float(r["ST"]) for r in sobol}
MU = {(r["species"], r["param"]): float(r["mu_star"]) for r in morris}
# order by mean ST across species
order = sorted(params, key=lambda p: -sum(ST.get((s, p), 0) for s in SPECIES))
r2 = stab.get("emulator_R2", {})

L = []
L.append(r"\begin{table}[htbp]")
L.append(r"\centering")
st = stab["stability"]
pct = 100 * st["frac_bowhead_below_both"]
L.append(
    r"\caption{Global sensitivity analysis of surface heat loss $\log Q_{\mathrm{surface}}$ "
    r"at $T_w=-2\,^\circ$C. For each of the eight uncertain inputs and three species we "
    r"report the Morris mean-absolute elementary effect $\mu^{*}$ (screening) and the "
    r"Sobol total-effect index $S_T$ from a polynomial emulator; parameters are ordered by "
    r"mean $S_T$. In a "
    + f"{st['n_samples']}"
    + r"-sample Latin-hypercube ensemble over all eight parameters the balaenid (bowhead) "
      r"retained a lower mass-specific surface cost than both rorquals in "
    + f"{pct:.0f}"
    + r"\% of draws, confirming the family ranking is robust to joint parameter uncertainty. "
      r"Emulator held-out $R^2$: "
    + ", ".join(f"{SP_LABEL[s]} {r2[s]['r2_test']:.2f}" for s in SPECIES if s in r2)
    + r".}"
)
L.append(r"\label{tab:gsa}")
L.append(r"\begin{tabular}{l" + "rr" * len(SPECIES) + "}")
L.append(r"\hline")
L.append(r"& " + " & ".join(
    r"\multicolumn{2}{c}{" + SP_LABEL[s] + "}" for s in SPECIES) + r" \\")
L.append("Parameter & " + " & ".join(r"$\mu^{*}$ & $S_T$" for _ in SPECIES) + r" \\")
L.append(r"\hline")
for p in order:
    cells = []
    for s in SPECIES:
        cells.append(f"{MU.get((s, p), 0):.2f}")
        cells.append(f"{ST.get((s, p), 0):.2f}")
    L.append(PLABEL.get(p, p) + " & " + " & ".join(cells) + r" \\")
L.append(r"\hline")
L.append(r"\end{tabular}")
L.append(r"\end{table}")

OUT.write_text("\n".join(L) + "\n")
print(f"wrote {OUT}")
print("param order (by mean ST):", ", ".join(order))
