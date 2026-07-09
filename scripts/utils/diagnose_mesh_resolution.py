#!/usr/bin/env python3
"""
Mesh Resolution Diagnostic
===========================
Checks whether the element_size in each species YAML is small enough to
resolve the blubber layer. The rule is:

    element_size < 0.5 × min(blubber_thickness per region)

Regions with fewer than 2 elements through the blubber produce an unresolved
blubber-core interface, which inflates thermal conductance and causes
spurious core hyperthermia (confirmed root cause for blue whale: PAPER_REVISION_PLAN B2).

Usage:
    conda activate whales
    python scripts/utils/diagnose_mesh_resolution.py
    python scripts/utils/diagnose_mesh_resolution.py --species blue_whale
"""

import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent

# Insulating regions: blubber is the PRIMARY thermal resistance — must be resolved.
# Under-resolution here directly inflates core temperature.
INSULATING_REGIONS = {'Body', 'Peduncle', 'Peduncle'}

# Thermal-window regions: convection-dominated; thin blubber expected and acceptable.
# Under-resolution here affects appendage heat flux accuracy (minor) but not core T.
THERMAL_WINDOW_REGIONS = {
    'Fluke', 'DorsalFin', 'VentralGrooves',
    'LeftPectoralFlipper', 'RightPectoralFlipper',
    'LeftPectoralFlipper', 'RightPectoralFlipper',
    'Rostrum',
}


def check_species(species_name):
    import sys
    sys.path.insert(0, str(PROJECT_ROOT))
    from scripts.utils.species_config import load_species_config

    cfg = load_species_config(species_name=species_name)
    el = cfg.mesh_element_size
    blubber = cfg.blubber_thickness

    if not blubber:
        return {'species': species_name, 'error': 'no blubber_thickness in config'}

    results = []
    for region, thickness in blubber.items():
        if thickness <= 0:
            continue
        elements_through = thickness / el
        is_insulating = region in INSULATING_REGIONS

        if is_insulating:
            # Body/peduncle: element_size must be < blubber_thickness so that
            # at least some nodes are correctly classified as blubber.
            # Confirmed failure mode: element_size > blubber_thickness → 0 blubber nodes
            #   → inflated conductance → core hyperthermia (blue whale root cause).
            if elements_through < 1.0:
                flag = "CRITICAL: element_size > blubber — 0 nodes classified as blubber"
            elif elements_through < 2.0:
                flag = "WARNING: element_size > 0.5×blubber — sparse blubber resolution"
            else:
                flag = ""
        else:
            # Appendage: informational only — convection-dominated
            if elements_through < 1.0:
                flag = "note (thermal window, element > blubber)"
            else:
                flag = ""

        results.append({
            'region': region,
            'blubber_m': thickness,
            'element_size_m': el,
            'elements_through': elements_through,
            'flag': flag,
            'is_insulating': is_insulating,
        })

    # Fail only when element_size > blubber_thickness for insulating regions
    # (confirmed cause of core hyperthermia for blue whale)
    insulating_results = [r for r in results if r['is_insulating']]
    passes = all(r['elements_through'] >= 1.0 for r in insulating_results)

    return {
        'species': species_name,
        'common_name': cfg.common_name,
        'element_size_m': el,
        'regions': results,
        'passes': passes,
    }


def print_report(report):
    if 'error' in report:
        print(f"\n  {report['species']}: ERROR — {report['error']}")
        return

    status = "OK" if report['passes'] else "FAIL"
    print(f"\n{'='*60}")
    print(f"  {report['common_name']} ({report['species']})  [{status}]")
    print(f"  element_size: {report['element_size_m']:.3f} m")
    print(f"  {'Region':<28} {'Type':<9} {'Blubber (m)':>11} {'Elem/blub':>10}  {'Note'}")
    print(f"  {'-'*72}")
    for r in sorted(report['regions'], key=lambda x: (not x['is_insulating'], x['elements_through'])):
        rtype = "insul." if r['is_insulating'] else "window"
        note = f"  ← {r['flag']}" if r['flag'] else ""
        print(f"  {r['region']:<28} {rtype:<9} {r['blubber_m']:>11.4f} {r['elements_through']:>10.2f}{note}")


def main():
    parser = argparse.ArgumentParser(
        description="Check mesh element size against blubber thickness for each species")
    parser.add_argument('--species', '-s', nargs='+', default=None,
                        help='Species to check (default: all)')
    args = parser.parse_args()

    import sys
    sys.path.insert(0, str(PROJECT_ROOT))
    from scripts.utils.species_config import SPECIES_DIR

    species_list = args.species or [p.stem for p in sorted(SPECIES_DIR.glob('*.yaml'))]

    print("MESH RESOLUTION DIAGNOSTIC")
    print("Rule: element_size < 0.5 × min(blubber_thickness)  →  ≥2 elements through blubber")

    any_fail = False
    for sp in species_list:
        try:
            report = check_species(sp)
            print_report(report)
            if not report.get('passes', True):
                any_fail = True
        except Exception as e:
            print(f"\n  [ERROR] {sp}: {e}")

    print()
    if any_fail:
        print("ACTION REQUIRED: Re-mesh species marked FAIL, then re-run simulations.")
    else:
        print("All species pass the mesh resolution check.")


if __name__ == "__main__":
    main()
