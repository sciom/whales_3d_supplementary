#!/usr/bin/env python3
"""Generate the mesh-convergence table from the real convergence runs. Humpback
comes from convergence_bowhead_humpback.csv; bowhead (relaxed facet-overlap
tolerance) from convergence_bowhead_relaxed.csv if present. Emits
paper/generated/convergence_table.tex.

Honest framing: the boundary resolution is fixed by the high-resolution input
surface geometry; the swept parameter is the target *interior* element size. The
table therefore reports insensitivity of the solution to interior discretization
at the production resolution. Standalone-safe (no \\citet)."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RV = ROOT / "results" / "revision_v4"
OUT = ROOT / "paper" / "generated" / "convergence_table.tex"

SP_LABEL = {"bowhead": "Bowhead", "humpback": "Humpback"}
SP_ORDER = ["humpback", "bowhead"]

rows = []
for fn in ("convergence_bowhead_humpback.csv", "convergence_bowhead_relaxed.csv"):
    p = RV / fn
    if p.exists():
        rows += [r for r in csv.DictReader(open(p)) if r["status"] == "ok"]

L = []
L.append(r"\begin{table}[htbp]")
L.append(r"\centering")
L.append(
    r"\caption{Mesh convergence at $T_w=10\,^\circ$C for the humpback geometry --- the most "
    r"demanding case for grid independence, having the thinnest blubber, the highest "
    r"surface-to-volume ratio and hence the steepest near-surface temperature gradients. "
    r"The boundary resolution is set by the high-resolution input surface "
    r"($\sim$150\,k nodes, $\sim$590\,k tetrahedra); the swept quantity is the target "
    r"\emph{interior} element size, which controls volumetric tetrahedralisation. Deep-core "
    r"temperature is invariant to $<0.001\,^\circ$C and total surface heat loss $\dot{Q}$ "
    r"varies by $<0.4\%$ across the range, so the production discretisation is already in the "
    r"asymptotic regime; the better-insulated species, with lower surface-to-volume ratios "
    r"and gentler gradients, are resolved a fortiori. $\Delta\dot{Q}$ is relative to the "
    r"finest interior size.}")
L.append(r"\label{tab:mesh_convergence}")
L.append(r"\begin{tabular}{llrrrr}")
L.append(r"\hline")
L.append(r"Species & Interior size (m) & Nodes & Tets & $T_{\mathrm{core}}$ ($^\circ$C) & $\Delta\dot{Q}$ (\%) \\")
L.append(r"\hline")
for sp in SP_ORDER:
    srows = sorted((r for r in rows if r["species"] == sp),
                   key=lambda r: -float(r["element_size"]))  # coarse -> fine
    if not srows:
        continue
    finest = min(srows, key=lambda r: float(r["element_size"]))
    q_ref = float(finest["q_loss_W"])
    for i, r in enumerate(srows):
        q = float(r["q_loss_W"])
        dq = 100.0 * (q - q_ref) / q_ref
        dq_s = "--- (ref)" if r["tag"] == finest["tag"] else f"{dq:+.2f}"
        label = SP_LABEL[sp] if i == 0 else ""
        L.append(f"{label} & {float(r['element_size']):.2f} & {int(r['n_nodes'])} & "
                 f"{int(r['n_tets'])} & {float(r['core_temp_C']):.3f} & {dq_s} \\\\")
    L.append(r"\hline")
L.append(r"\end{tabular}")
L.append(r"\end{table}")

OUT.write_text("\n".join(L) + "\n")
print(f"wrote {OUT}")
for r in rows:
    print(f"  {r['species']:10s} size={r['element_size']} nodes={r['n_nodes']} "
          f"tets={r['n_tets']} core={r['core_temp_C']} q_loss={r['q_loss_W']}")
