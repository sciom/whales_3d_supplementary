"""
Result discovery utilities for postprocessing.
"""

import json
from pathlib import Path
import pyvista as pv


def load_metadata(meta_file: Path):
    with open(meta_file, "r", encoding="utf-8") as f:
        return json.load(f)


def load_latest_solution(results_dir: Path):
    """Load the most recent single-solution result."""
    vtu_files = sorted(results_dir.glob("thermal_solution_*.vtu"))
    if not vtu_files:
        raise FileNotFoundError("No simulation results found!")

    latest = vtu_files[-1]
    mesh = pv.read(latest)
    metadata = load_metadata(latest.with_name(latest.stem + "_metadata.json"))
    return latest, mesh, metadata


def iter_case_metadata(results_dir: Path, pattern: str = "thermal_T*C_*.vtu"):
    """Yield (vtu_file, metadata) for all matching cases."""
    for vtu_file in results_dir.glob(pattern):
        meta_file = vtu_file.with_name(vtu_file.stem + "_metadata.json")
        if not meta_file.exists():
            continue
        yield vtu_file, load_metadata(meta_file)


def find_latest_case(results_dir: Path, T_water_C=None, mode=None):
    """Find the latest multicase result matching water temp and mode."""
    candidates = []
    for vtu_file, metadata in iter_case_metadata(results_dir):
        if T_water_C is not None and metadata.get("T_water_C") != T_water_C:
            continue
        if mode is not None and metadata.get("mode") != mode:
            continue
        candidates.append((metadata.get("timestamp", ""), vtu_file, metadata))

    if not candidates:
        raise FileNotFoundError("No matching multicase results found.")

    candidates.sort()
    _, vtu_file, metadata = candidates[-1]
    mesh = pv.read(vtu_file)
    return vtu_file, mesh, metadata


def load_heat_flux_analysis(results_dir: Path, timestamp: str):
    """Load heat flux analysis JSON for a given timestamp."""
    flux_file = results_dir / f"heat_flux_analysis_{timestamp}.json"
    with open(flux_file, "r", encoding="utf-8") as f:
        return json.load(f)
