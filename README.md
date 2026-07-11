# Supplementary materials for the mysticete thermoregulation FEM study

This repository contains curated supplementary material for the manuscript:

**Thermal maintenance under a defended arterial set-point: a comparative 3D
finite-element bioheat analysis of mysticete cold-water thermoregulation**

The repository is intended for journal submission and peer review.  It contains
processed simulation outputs, the global sensitivity and closed-reservoir
analyses added during revision, model-support and parameter-provenance tables,
one complete raw worked example, selected validation figures, and scripts needed
to regenerate figure and table outputs from the processed data.

## Contents

- `data/processed/` — processed multicase result summaries for each whale model
  (all modes, including `moving` for the swimming analysis).
- `data/gsa/` — **global sensitivity analysis** (revision): Morris screening,
  Sobol indices, emulator ranking-stability (2000-sample LHS) and the raw 576
  per-solve records. See `gsa_summary.md`.
- `data/closed_reservoir/` — **closed-reservoir variant** (revision): emergent
  core temperature when arterial temperature is energy-limited rather than fixed.
- `data/convergence/` — mesh-convergence, high-conductance robustness test, and
  quantitative geometry/mass-consistency check (revision).
- `data/derived_tables/` — machine-readable versions of the manuscript tables
  (energy partition, provenance, geometry, GSA, swimming, convergence, `T_tn`).
- `data/sensitivity/` — one-at-a-time sensitivity-audit summaries and configs.
- `data/model_support/` — literature-derived model-support tables.
- `data/species/` — species parameter files used by the simulations.
- `data/csv_tables/` — structured parameter tables used during model assembly.
- `data/raw_example/` — **one complete raw volumetric example** (right whale:
  input mesh + solved field + config + regeneration instructions).
- `figures/` — exported figure files and mesh-validation images, including the
  new GSA, closed-reservoir and swimming figures.
- `scripts/` — postprocessing utilities to regenerate figure and table outputs.
- `docs/` — reviewer/maintainer notes and the file-to-manuscript manifest.

Large raw solver files (`*.vtu`, full transient fields, and duplicate historical
metadata exports) are otherwise excluded; a single complete example is provided
in `data/raw_example/`, and any case can be regenerated from the main modelling
repository.

## Reproducing the revision analyses

In the `humpback_fem` conda environment, from the main modelling repository the
new analyses are produced by:

- Global sensitivity: `python revision_v4/gsa.py {generate,run,analyze}`
- Closed reservoir: `python revision_v4/closed_reservoir.py`
- Mesh convergence / high-h / geometry: `python revision_v4/convergence.py`

The postprocessing scripts in `scripts/postprocessing/` (e.g. `plot_gsa.py`,
`build_gsa_tables.py`, `plot_closed_reservoir.py`, `build_swimming.py`,
`build_convergence_table.py`, `build_energy_partition.py`,
`build_provenance_tables.py`) are archived here as the exact provenance of each
manuscript figure and table. They reference the main modelling repository's
directory layout (`results/`, `paper/`); the processed inputs they consume are
mirrored in this repository under `data/` (see `docs/manuscript_manifest.md` for
the file-to-figure/table mapping).

## Reproducibility Scope

The files here support reproduction of the reported tables, figures and
postprocessed biomarkers from processed model outputs, plus one end-to-end raw
example. They are not a full archive of the modelling source tree or every raw
FEM field, which are regenerable from the main repository.

## Citation

If this repository is archived through Zenodo, Figshare or another repository
service, keep `CITATION.cff` in sync with the assigned DOI before submission.
