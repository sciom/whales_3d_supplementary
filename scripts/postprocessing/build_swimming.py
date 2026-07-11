#!/usr/bin/env python3
"""Swimming (moving-mode) analysis for the manuscript. Compares the resting active
baseline against the swimming (moving) mode at the baseline food/condition scenario
across water temperature, for the species with full moving-mode coverage.

Two messages, both data-backed:
 (a) the defended core is nearly unchanged by swimming (<0.2 C) in all species except
     the smallest (minke, 7 t), where added muscular heat raises core by ~0.7-1.0 C;
 (b) swimming raises total surface heat loss by 4-28 %, and the fractional increase
     grows with water temperature -- swimming loads both sides of the heat budget and
     is not a thermoregulatory subsidy.

Emits paper/figures/swimming.pdf and paper/generated/swimming_table.tex.
"""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUTFIG = ROOT / "paper" / "figures" / "swimming.pdf"
OUTTAB = ROOT / "paper" / "generated" / "swimming_table.tex"

SPECIES = ["humpback", "right_whale", "blue_whale", "bowhead", "fin_whale", "minke"]
LABEL = {"humpback": "Humpback", "right_whale": "Right whale", "blue_whale": "Blue whale",
         "bowhead": "Bowhead", "fin_whale": "Fin whale", "minke": "Minke"}
COLOR = {"humpback": "#c44e34", "right_whale": "#7b4fa3", "blue_whale": "#1b6f8c",
         "bowhead": "#2e7d5b", "fin_whale": "#d9a441", "minke": "#b23b6a"}
MARK = {"humpback": "o", "right_whale": "s", "blue_whale": "^", "bowhead": "D",
        "fin_whale": "v", "minke": "*"}


def load(sp):
    d = json.load(open(ROOT / "results" / sp / "multicase_comparison_data.json"))
    cases = d["cases"]

    def base(mode):
        return {c["T_water_C"]: c for c in cases if c.get("mode") == mode
                and c.get("food_availability", 1) == 1.0
                and c.get("body_condition_index", 1) == 1.0}
    a, m = base("active"), base("moving")
    tw = sorted(set(a) & set(m))
    rows = []
    for t in tw:
        ca, cm = a[t], m[t]
        qa = ca.get("total_heat_loss_W", 0.0)
        qm = cm.get("total_heat_loss_W", 0.0)
        rows.append(dict(tw=t, core_a=ca["core_temp_C"], core_m=cm["core_temp_C"],
                         dcore=cm["core_temp_C"] - ca["core_temp_C"],
                         qa=qa, qm=qm, ratio=(qm / qa if qa else float("nan"))))
    return rows


data = {sp: load(sp) for sp in SPECIES}

# ---------------- figure ----------------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.8, 3.9))
for sp in SPECIES:
    r = data[sp]
    tw = [x["tw"] for x in r]
    ms = 8 if sp == "minke" else 4.5
    axA.plot(tw, [x["dcore"] for x in r], "-", marker=MARK[sp], ms=ms, lw=1.6,
             color=COLOR[sp], label=LABEL[sp])
    axB.plot(tw, [100 * (x["ratio"] - 1) for x in r], "-", marker=MARK[sp], ms=ms,
             lw=1.6, color=COLOR[sp], label=LABEL[sp])

axA.axhline(0, color="#9aa0a6", lw=0.8, ls=":")
axA.set_xlabel("Water temperature (°C)", fontsize=9)
axA.set_ylabel("Swimming $-$ resting core temperature (°C)", fontsize=9)
axA.set_title("(a) Core temperature is defended", fontsize=9, loc="left")
axA.tick_params(labelsize=8)
axA.grid(alpha=0.25, lw=0.5)
axA.annotate("minke (7 t): added heat\nnot fully shed", xy=(10, 0.76), xytext=(13, 0.45),
             fontsize=7, color=COLOR["minke"],
             arrowprops=dict(arrowstyle="->", color=COLOR["minke"], lw=0.8))

axB.set_xlabel("Water temperature (°C)", fontsize=9)
axB.set_ylabel("Increase in surface heat loss with swimming (%)", fontsize=9)
axB.set_title("(b) But surface heat loss rises", fontsize=9, loc="left")
axB.tick_params(labelsize=8)
axB.grid(alpha=0.25, lw=0.5)
axB.legend(fontsize=7, loc="upper left", ncol=2, framealpha=0.9)

fig.tight_layout()
fig.savefig(OUTFIG, bbox_inches="tight")
print(f"wrote {OUTFIG}")

# ---------------- table (compact: -2 C and 30 C endpoints) ----------------
L = []
L.append(r"\begin{table}[htbp]")
L.append(r"\centering")
L.append(
    r"\caption{Swimming (moving mode) versus resting (active mode) at the baseline "
    r"food/condition scenario, at the cold and warm ends of the range. $\Delta T_{\mathrm{core}}$ "
    r"is the swimming minus resting defended core temperature; $\dot{Q}$ is total surface heat "
    r"loss. Swimming leaves the defended core almost unchanged in all species except the "
    r"smallest (minke), while raising surface heat loss by 4--28\,\%, confirming that added "
    r"muscular heat is shed rather than banked --- swimming is not a thermoregulatory subsidy.}")
L.append(r"\label{tab:swimming}")
L.append(r"\begin{tabular}{lrrrr}")
L.append(r"\hline")
L.append(r"& \multicolumn{2}{c}{$T_w=-2\,^\circ$C} & \multicolumn{2}{c}{$T_w=30\,^\circ$C} \\")
L.append(r"Species & $\Delta T_{\mathrm{core}}$ (°C) & $\Delta\dot{Q}$ (\%) & $\Delta T_{\mathrm{core}}$ (°C) & $\Delta\dot{Q}$ (\%) \\")
L.append(r"\hline")
for sp in SPECIES:
    r = {x["tw"]: x for x in data[sp]}
    cold = r.get(-2.0) or r.get(min(r))
    warm = r.get(30.0) or r.get(max(r))
    L.append(f"{LABEL[sp]} & {cold['dcore']:+.2f} & {100*(cold['ratio']-1):+.1f} & "
             f"{warm['dcore']:+.2f} & {100*(warm['ratio']-1):+.1f} \\\\")
L.append(r"\hline")
L.append(r"\end{tabular}")
L.append(r"\end{table}")
OUTTAB.write_text("\n".join(L) + "\n")
print(f"wrote {OUTTAB}")
for sp in SPECIES:
    r = data[sp]
    print(f"  {sp:12s} dcore range {min(x['dcore'] for x in r):+.2f}..{max(x['dcore'] for x in r):+.2f} C, "
          f"Q +{100*(min(x['ratio'] for x in r)-1):.0f}..{100*(max(x['ratio'] for x in r)-1):.0f}%")
