#!/usr/bin/env python3
"""
Results Inventory
=================
Scans all species result directories and reports which
(species, water_temp_C, mode) combinations have VTU files + metadata JSONs.
Highlights missing combinations relative to the expected 4-mode × 9-temp grid.

Usage:
    conda activate whales
    python scripts/utils/results_inventory.py
    python scripts/utils/results_inventory.py --species humpback
"""

import argparse
import json
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_BASE = PROJECT_ROOT / "results"

EXPECTED_TEMPS = [-2.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 28.0, 30.0]
EXPECTED_MODES = ['active', 'no_control', 'moving', 'moving_no_control']

# Pattern: thermal_T<temp>C_<mode>_<timestamp>.vtu
_VTU_RE = re.compile(r"thermal_T(-?\d+(?:\.\d+)?)C_([\w]+)_\d{8}_\d{6}\.vtu$")


def scan_species(species_name):
    results_dir = RESULTS_BASE / species_name
    if not results_dir.exists():
        return {'species': species_name, 'error': 'results directory not found'}

    found = {}  # (temp, mode) -> {'vtu': bool, 'meta': bool}

    for vtu in results_dir.glob("thermal_T*C_*.vtu"):
        m = _VTU_RE.match(vtu.name)
        if not m:
            continue
        temp = float(m.group(1))
        mode = m.group(2)
        key = (temp, mode)
        if key not in found:
            found[key] = {'vtu': False, 'meta': False}
        found[key]['vtu'] = True
        # Check paired metadata JSON
        meta_path = vtu.with_name(vtu.stem + "_metadata.json")
        if meta_path.exists():
            found[key]['meta'] = True

    # Check for gaps against expected grid
    gaps = []
    for t in EXPECTED_TEMPS:
        for mode in EXPECTED_MODES:
            key = (t, mode)
            if key not in found:
                gaps.append(key)
            elif not found[key]['meta']:
                gaps.append(('META_MISSING', t, mode))

    return {
        'species': species_name,
        'found': found,
        'gaps': gaps,
        'n_found': len(found),
        'n_expected': len(EXPECTED_TEMPS) * len(EXPECTED_MODES),
    }


def print_report(report):
    if 'error' in report:
        print(f"\n  {report['species']}: {report['error']}")
        return

    n = report['n_found']
    exp = report['n_expected']
    status = "COMPLETE" if n >= exp and not report['gaps'] else f"{n}/{exp}"
    print(f"\n{'='*60}")
    print(f"  {report['species']:20s}  [{status}]")

    if report['gaps']:
        print(f"  Missing combinations:")
        for g in sorted(report['gaps']):
            if g[0] == 'META_MISSING':
                print(f"    META MISSING   T={g[1]:5.1f}°C  {g[2]}")
            else:
                print(f"    NOT FOUND      T={g[0]:5.1f}°C  {g[1]}")
    else:
        print(f"  All {exp} expected cases present.")


def main():
    parser = argparse.ArgumentParser(
        description="Inventory simulation result files for all species")
    parser.add_argument('--species', '-s', nargs='+', default=None,
                        help='Species to scan (default: all with result dirs)')
    args = parser.parse_args()

    if args.species:
        species_list = args.species
    else:
        import sys
        sys.path.insert(0, str(PROJECT_ROOT))
        from scripts.utils.species_config import SPECIES_DIR
        species_list = [p.stem for p in sorted(SPECIES_DIR.glob('*.yaml'))]

    print("RESULTS INVENTORY")
    print(f"Expected: {len(EXPECTED_TEMPS)} temps × {len(EXPECTED_MODES)} modes = "
          f"{len(EXPECTED_TEMPS)*len(EXPECTED_MODES)} cases per species")

    for sp in species_list:
        report = scan_species(sp)
        print_report(report)

    print()


if __name__ == "__main__":
    main()
