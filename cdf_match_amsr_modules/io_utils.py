from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable


SLOT_PATTERN = re.compile(r"^soil_moisture_([AD])_(\d{2})$")
SEC_PER_DAY = 24 * 60 * 60
EPOCH = 0


def list_input_files(path: Path) -> list[Path]:
    if path.is_file():
        if path.suffix.lower() != ".nc":
            raise ValueError(f"Input path must be .nc: {path}")
        return [path]

    files = sorted(path.glob("*.nc"))
    if not files:
        raise FileNotFoundError(f"No NetCDF files found under: {path}")
    return files


def collect_slot_pairs(ds) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for name in sorted(ds.variables.keys()):
        if SLOT_PATTERN.match(name):
            suffix = name.replace("soil_moisture_", "", 1)
            tm_name = f"observation_time_min_{suffix}"
            if tm_name in ds.variables:
                out.append((name, tm_name))
    return out


def normalize_reference_lon(q_lon: Iterable[float], source_lon: Iterable[float]) -> list[float]:
    q = list(float(x) for x in q_lon)
    src = list(float(x) for x in source_lon)
    if not src:
        return q

    src_min = min(src)
    src_max = max(src)
    if src_min >= -180.0 and src_max <= 180.0:
        # Reference axis is centered at 0 (e.g., -179.75..179.75).
        return [((x + 180.0) % 360.0) - 180.0 for x in q]

    # Reference axis is likely 0..360.
    return [x % 360.0 for x in q]


def nearest_index_1d(axis: list[float] | tuple[float, ...] | object, values):
    axis_arr = list(float(v) for v in axis)
    if not axis_arr:
        raise ValueError("axis must not be empty")
    vals = list(float(v) for v in values)
    asc = axis_arr if axis_arr[0] < axis_arr[-1] else list(reversed(axis_arr))

    import numpy as np

    a = np.array(asc, dtype=np.float64)
    v = np.array(vals, dtype=np.float64)
    pos = np.searchsorted(a, v, side="left")
    left = np.clip(pos - 1, 0, len(a) - 1)
    right = np.clip(pos, 0, len(a) - 1)
    choose_left = (v - a[left]) <= (a[right] - v)
    idx_asc = np.where(choose_left, left, right)
    if axis_arr[0] < axis_arr[-1]:
        return idx_asc.astype(np.int64)
    return (len(a) - 1 - idx_asc).astype(np.int64)
