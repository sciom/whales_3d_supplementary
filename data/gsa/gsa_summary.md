# Global Sensitivity & Uncertainty Analysis (GSA)

Water temperature: **Tw = -2 C** (cold, costly end). Response of interest: **log(Q_surface)** = log(energy['q_loss_W']).

Method: Morris elementary-effects screening (r trajectories), a polynomial emulator (main + quadratic + 2-way) of log(Q_surface), Sobol variance decomposition on the emulator, and a 2000-sample LHS ensemble for family-ranking stability.


## Emulator fit (log Q_surface)

| species | R2 train | R2 5-fold CV | R2 held-out LHS |
|---|---|---|---|
| humpback | 1.000 | 1.000 | 1.000 |
| fin_whale | 0.999 | 0.999 | 0.998 |
| bowhead | 0.998 | 0.995 | 0.993 |

## Top-3 most influential parameters (Sobol total-effect ST)

| species | 1st | 2nd | 3rd |
|---|---|---|---|
| humpback | ava (ST=0.68) | omega_core (ST=0.13) | h (ST=0.09) |
| fin_whale | ava (ST=0.56) | k_blub (ST=0.12) | omega_core (ST=0.09) |
| bowhead | Ttn (ST=0.63) | ava (ST=0.30) | k_blub (ST=0.03) |

## Morris mu* ranking (screening; log Q_surface)

- **humpback**: ava (0.33), omega_core (0.15), h (0.12), k_blub (0.09)
- **fin_whale**: ava (0.30), k_blub (0.13), omega_core (0.12), h (0.12)
- **bowhead**: Ttn (0.64), ava (0.44), k_blub (0.14), h (0.13)

## Full Sobol indices

| species | param | S1 | ST |
|---|---|---|---|
| humpback | ava | 0.677 | 0.681 |
| humpback | omega_core | 0.129 | 0.132 |
| humpback | h | 0.085 | 0.086 |
| humpback | k_blub | 0.061 | 0.061 |
| humpback | BCI | 0.041 | 0.041 |
| humpback | cche | 0.003 | 0.003 |
| humpback | fmr | 0.000 | 0.000 |
| humpback | Ttn | 0.000 | 0.000 |
| fin_whale | ava | 0.552 | 0.556 |
| fin_whale | k_blub | 0.124 | 0.124 |
| fin_whale | omega_core | 0.089 | 0.091 |
| fin_whale | BCI | 0.082 | 0.082 |
| fin_whale | h | 0.081 | 0.082 |
| fin_whale | Ttn | 0.068 | 0.069 |
| fin_whale | cche | 0.000 | 0.000 |
| fin_whale | fmr | 0.000 | 0.000 |
| bowhead | Ttn | 0.627 | 0.631 |
| bowhead | ava | 0.296 | 0.299 |
| bowhead | k_blub | 0.028 | 0.031 |
| bowhead | h | 0.023 | 0.024 |
| bowhead | BCI | 0.018 | 0.019 |
| bowhead | cche | 0.001 | 0.001 |
| bowhead | omega_core | 0.001 | 0.001 |
| bowhead | fmr | 0.000 | 0.000 |

## Ranking stability (headline)

Across a **2000-sample** LHS ensemble over all 8 uncertain parameters, the balaenid (bowhead) had **lower mass-specific surface cost than BOTH rorquals (humpback & fin)** in **100.0%** of samples at Tw=-2 C.

- bowhead < humpback: 100.0%
- bowhead < fin: 100.0%
- rank-reversal frequency (bowhead not below both): 0.0%

Mean mass-specific surface cost (W/kg): humpback=0.9250, fin_whale=0.6103, bowhead=0.2833

## Validation: monotonic response checks (rank correlation, log Q)

Expect POSITIVE for both (higher h or k_blub -> more heat loss).

| species | rho(h) | rho(k_blub) |
|---|---|---|
| humpback | +0.149 | +0.501 |
| fin_whale | +0.193 | +0.569 |
| bowhead | +0.143 | +0.249 |

## Runs

Total FEM solves: **576** (humpback=192, fin_whale=192, bowhead=192).
