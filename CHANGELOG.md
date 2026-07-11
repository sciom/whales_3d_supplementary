# Changelog

All notable changes to this supplementary-materials package are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — 2026-07-11

First archived release, accompanying the peer-review submission of the manuscript
*"Thermal maintenance under a defended arterial set-point: a comparative 3D
finite-element bioheat analysis of mysticete cold-water thermoregulation"*
(Hackenberger, Djerdj & Hackenberger). Archived at DOI
[10.5281/zenodo.21308022](https://doi.org/10.5281/zenodo.21308022) under CC-BY-4.0.

### Added
- Processed multicase result summaries for all eight whale models, all thermal
  modes including `moving` (swimming) — `data/processed/`.
- Global sensitivity analysis: Morris screening, Sobol indices, emulator
  ranking-stability (2000-sample LHS), and the raw 576 per-solve records —
  `data/gsa/`.
- Closed-reservoir variant: emergent core temperature under an energy-limited
  arterial temperature — `data/closed_reservoir/`.
- Mesh-convergence, high-conductance robustness, and geometry/mass-consistency
  checks — `data/convergence/`.
- Machine-readable derived tables mirroring the manuscript tables (energy
  partition, provenance, geometry, GSA, swimming, convergence, `T_tn`) —
  `data/derived_tables/`.
- One complete raw worked example (right whale: input mesh, solved field, config,
  regeneration instructions) — `data/raw_example/`.
- One-at-a-time sensitivity audits (`data/sensitivity/`), literature-derived
  model-support tables (`data/model_support/`), species parameter files
  (`data/species/`), and structured assembly tables (`data/csv_tables/`).
- Exported figures and mesh-validation images, including the GSA,
  closed-reservoir, and swimming figures — `figures/`.
- Postprocessing utilities archived as the provenance of each manuscript figure
  and table — `scripts/`.
- Reviewer/maintainer notes and the file-to-manuscript manifest — `docs/`.
- `CITATION.cff` with the assigned Zenodo DOI, and this changelog.

### Notes
- Scope: these files reproduce the reported tables, figures, and postprocessed
  biomarkers from processed model outputs, plus one end-to-end raw example. They
  are not a full archive of the modelling source tree or every raw FEM field,
  which are regenerable from the main modelling repository.
- Large raw solver files (`*.vtu`, full transient fields, duplicate historical
  metadata exports) are excluded except for the single provided example.

[1.0.0]: https://doi.org/10.5281/zenodo.21308022
