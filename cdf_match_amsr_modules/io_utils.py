from __future__ import annotations

import re
from pathlib import Path

import numpy as np


SEC_PER_DAY = 24 * 60 * 60
SLOT_PATTERN = re.compile(r"^soil_moisture_([AD])_(\d{2})$")


def list_input_files(path: Path) -> list[Path]:
    if path.is_file():
        if path.suffix.lower() != ".nc":
            raise ValueError(f"Input path must be .nc: {path}")
        return [path]
    if path.suffix.lower() == ".nc":
        raise FileNotFoundError(f"Input NetCDF file not found: {path}")
    files = sorted(path.glob("*.nc"))
    if not files:
        raise FileNotFoundError(f"No NetCDF files found under: {path}")
    return files


def collect_slot_pairs(ds) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for name in sorted(ds.variables):
        if not SLOT_PATTERN.match(name):
            continue
        suffix = name.replace("soil_moisture_", "", 1)
        tm_name = f"observation_time_min_{suffix}"
        if tm_name in ds.variables:
            out.append((name, tm_name))
    if out:
        return out
    if "soil_moisture" in ds.variables and "observation_time_min" in ds.variables:
        return [("soil_moisture", "observation_time_min")]
    return out


def normalize_lon(values: np.ndarray, target_axis: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    axis = np.asarray(target_axis, dtype=np.float64)
    if np.nanmin(axis) >= -180.0 and np.nanmax(axis) <= 180.0:
        return ((values + 180.0) % 360.0) - 180.0
    return values % 360.0


def nearest_index_1d(axis: np.ndarray, values: np.ndarray) -> np.ndarray:
    axis = np.asarray(axis, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if axis.ndim != 1 or axis.size == 0:
        raise ValueError("axis must be a non-empty 1D array")

    ascending = axis[0] <= axis[-1]
    work_axis = axis if ascending else axis[::-1]
    pos = np.searchsorted(work_axis, values, side="left")
    left = np.clip(pos - 1, 0, work_axis.size - 1)
    right = np.clip(pos, 0, work_axis.size - 1)
    choose_left = (values - work_axis[left]) <= (work_axis[right] - values)
    idx = np.where(choose_left, left, right)
    if ascending:
        return idx.astype(np.int64)
    return (work_axis.size - 1 - idx).astype(np.int64)


def input_to_m3m3(values: np.ndarray, units: str | None = None) -> np.ndarray:
    """Convert AMSR soil moisture to m3/m3.

    Existing AMSR products in this workspace are percent-like. If units are
    missing, values above 2 are also treated as percent for safety.
    """
    arr = np.asarray(values, dtype=np.float32)
    unit_text = (units or "").lower()
    if "%" in unit_text or "percent" in unit_text:
        return arr / np.float32(100.0)
    finite = arr[np.isfinite(arr)]
    if finite.size and float(np.nanmax(finite)) > 2.0:
        return arr / np.float32(100.0)
    return arr
