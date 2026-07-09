# Supplementary materials for the mysticete thermoregulation FEM study

This repository contains curated supplementary material for the manuscript:

**Thermal maintenance, not hypothermia, governs mysticete cold-water thermoregulation: a comparative 3D finite-element bioheat analysis**

The repository is intended for journal submission and peer review.  It contains
processed simulation outputs, sensitivity summaries, model-support tables,
selected validation figures, and scripts needed to regenerate the manuscript
figures from the processed outputs.

## Contents

- `manuscript/` - current manuscript PDF snapshot.
- `data/processed/` - processed multicase result summaries for each whale model.
- `data/sensitivity/` - sensitivity-analysis summaries and configuration files.
- `data/model_support/` - literature-derived model-support tables.
- `data/species/` - species parameter files used by the simulations.
- `data/csv_tables/` - structured parameter tables used during model assembly.
- `figures/` - exported manuscript figures and mesh-validation images.
- `scripts/` - postprocessing utilities required to regenerate figure outputs.
- `docs/` - notes for reviewers and maintainers.

Large raw solver files (`*.vtu`, full transient fields, and duplicate historical
metadata exports) are intentionally excluded.  They can be regenerated from the
main modelling repository if needed.

## Reproducibility Scope

The files here support reproduction of the reported tables, figures and
postprocessed biomarkers from processed model outputs.  They are not a full
archive of the modelling source tree or every raw FEM field.

## Citation

If this repository is archived through Zenodo, Figshare or another repository
service, update `CITATION.cff` with the assigned DOI before submission.
