# Model Support Data (Curated)

This folder contains curated parameter-support data used to justify or bound
inferred parameters in the whale thermal FEM model.

## Files

- `data/literature/model_support_data.csv`
  - Flat table of parameter values/ranges and source links.
  - Includes both primary measurements (e.g., digestibility) and engineering
    correlations used to justify scaling exponents.

## Notes

- Some quantities are explicitly derived (see `parameter` fields ending in `_derived`).
- Energy-density conversions for krill are stored as an envelope because both dry
  energy content and moisture fraction vary seasonally.
- This is not a feeding/intake model; it supports translating modeled thermal power
  (W) into an implied energy requirement (kcal/day) for scenario comparisons.

