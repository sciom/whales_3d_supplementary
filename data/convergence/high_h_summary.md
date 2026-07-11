# Deliverable 2 — High effective-conductance warm-water test

**Question (reviewer):** Is the "warm water is fine" result an artifact of the low
effective convection coefficient `h` used in the config (Body 35 … Fluke 60 W m⁻² K⁻¹)?
If `h` were realistically high, would warm water cause core overheating?

**Test:** Humpback (thin blubber) and bowhead (thick blubber), active mode, FA = BCI = 1,
at warm water T_w = 28 and 30 °C. The whole convection map was scaled up so that the
body coefficient reached ~200 W m⁻² K⁻¹ (×5.7) and ~500 W m⁻² K⁻¹ (×14.3), compared with
the baseline (body 35). Reported: volume-weighted core temp, deep-body temp, surface
heat loss `q_loss`, and whether the core rises materially above the defended set-point.

## Result

| Species | T_w (°C) | h level (body h) | core (°C) | deep body (°C) | q_loss (W) | core − set-point (°C) | overheating |
|---|---|---|---|---|---|---|---|
| Humpback | 28 | baseline (35) | 36.008 | 36.069 | 9 727 | +0.008 | no |
| Humpback | 28 | h200 | 35.998 | 36.067 | 12 555 | −0.002 | no |
| Humpback | 28 | h500 | 35.997 | 36.067 | 13 070 | −0.003 | no |
| Humpback | 30 | baseline (35) | 36.019 | 36.073 | 7 661 | +0.019 | no |
| Humpback | 30 | h200 | 36.011 | 36.071 | 10 042 | +0.011 | no |
| Humpback | 30 | h500 | 36.010 | 36.071 | 10 485 | +0.010 | no |
| Bowhead | 28 | baseline (35) | 33.949 | 34.058 | 4 852 | +0.149 | no |
| Bowhead | 28 | h200 | 33.947 | 34.058 | 6 474 | +0.147 | no |
| Bowhead | 28 | h500 | 33.947 | 34.058 | 6 791 | +0.147 | no |
| Bowhead | 30 | baseline (35) | 33.950 | 34.058 | 3 189 | +0.150 | no |
| Bowhead | 30 | h200 | 33.949 | 34.058 | 4 254 | +0.150 | no |
| Bowhead | 30 | h500 | 33.949 | 34.058 | 4 462 | +0.150 | no |

Set-points: humpback core 36.0 °C, bowhead 33.8 °C; lethal-hyperthermia guard 40.0 °C.

## Finding

**No — realistic high `h` does NOT produce warm-water overheating; the core stays defended.**
Raising the body coefficient from 35 to 500 W m⁻² K⁻¹ (a ~14× increase) leaves the core
temperature essentially unchanged: Δcore < 0.02 °C for the humpback and < 0.002 °C for the
bowhead across the entire `h` range, and the deep-body temperature moves by ≲0.005 °C. If
anything the higher `h` cools the animal *slightly* (core edges downward, heat loss rises by
25–40 %), because at T_w = 28–30 °C the water is still well below core temperature, so a
larger surface conductance can only extract *more* heat, never add it. Overheating would
require the water to exceed the core set-point (≥34–36 °C), which does not occur in any
oceanographically realistic scenario for these species.

The "warm water is fine" conclusion is therefore **not** an artifact of the chosen `h`: it is
robust to a >10× increase in effective surface conductance. Heat loss is the quantity that
responds to `h` (energy cost), not core temperature, consistent with the manuscript's thesis
that active thermoregulation manages heat *flux* rather than defending against warm-water
hyperthermia. The `q_loss` values (≈3–13 kW at these small water-to-core gradients) are
physically sane, confirming the solver operates in a correct physical (not normalized) frame.

Data: `results/revision_v4/high_h_test.csv`
