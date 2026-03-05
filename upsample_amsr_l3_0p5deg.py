#!/usr/bin/env python3
"""
Upscale merged AMSR-L3 NetCDF observations to a coarser regular grid.

Input format is expected to match merge_amsr_l3_daily.py:
  - dimensions: (time, lat, lon)
  - per-direction slot variables:
      soil_moisture_A_01, soil_moisture_A_02, ...
      soil_moisture_D_01, soil_moisture_D_02, ...
      with matching observation_time_min_* variables

This script maps each observed grid cell to the target 0.5° cell by nearest-cell
assignment and keeps duplicate observations in additional slot variables without
averaging.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Tuple

import numpy as np
from netCDF4 import Dataset

SLOT_PATTERN = re.compile(r"^soil_moisture_([AD])_(\d{2})$")


@dataclass(frozen=True)
class SourceChunk:
    file_path: Path
    time_index: int


class UpscaledNetCDFWriter:
    def __init__(
        self,
        output_path: Path,
        nlat: int,
        nlon: int,
        lat: np.ndarray,
        lon: np.ndarray,
        overwrite: bool,
    ) -> None:
        if output_path.exists() and not overwrite:
            raise FileExistsError(f"{output_path} already exists. Use --overwrite.")
        if output_path.exists():
            output_path.unlink()

        self.nlat = nlat
        self.nlon = nlon
        self._max_slot_by_direction: DefaultDict[str, int] = defaultdict(int)
        self._slots: Dict[tuple[str, int], Tuple[object, object]] = {}

        output_path.parent.mkdir(parents=True, exist_ok=True)
        self.ds = Dataset(output_path, "w", format="NETCDF4")
        self.ds.createDimension("time", None)
        self.ds.createDimension("lat", nlat)
        self.ds.createDimension("lon", nlon)

        self.var_time = self.ds.createVariable("time", "i4", ("time",))
        var_lat = self.ds.createVariable("lat", "f4", ("lat",))
        var_lon = self.ds.createVariable("lon", "f4", ("lon",))
        var_lat[:] = lat.astype(np.float32)
        var_lon[:] = lon.astype(np.float32)

        self.var_time.units = "days since 1970-01-01"
        var_lat.units = "degrees_north"
        var_lon.units = "degrees_east"
        self.time_index = 0

    def ensure_slot(self, direction: str, slot_id: int) -> tuple[object, object]:
        key = (direction, slot_id)
        if key in self._slots:
            return self._slots[key]

        var_sm_name = f"soil_moisture_{direction}_{slot_id:02d}"
        var_tm_name = f"observation_time_min_{direction}_{slot_id:02d}"

        var_sm = self.ds.createVariable(
            var_sm_name,
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, self.nlat), min(360, self.nlon)),
            fill_value=np.float32(np.nan),
        )
        var_tm = self.ds.createVariable(
            var_tm_name,
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, self.nlat), min(360, self.nlon)),
            fill_value=np.float32(np.nan),
        )

        var_sm.units = "%"
        var_tm.units = "minutes since 00:00 UTC of corresponding time day"
        var_sm.long_name = f"Soil Moisture Content ({direction})"
        var_tm.long_name = f"Observation Time ({direction})"
        var_sm.comment = (
            "All observations are preserved by assigning duplicates to higher slot variables."
        )

        self._slots[key] = (var_sm, var_tm)
        self._max_slot_by_direction[direction] = max(self._max_slot_by_direction[direction], slot_id)
        return var_sm, var_tm

    def write_time(self, day_value: int) -> None:
        self.var_time[self.time_index] = int(day_value)

    def next_time(self) -> int:
        current = self.time_index
        self.time_index += 1
        return current

    def close(self) -> None:
        self.ds.setncattr("max_slots_A", int(self._max_slot_by_direction.get("A", 0)))
        self.ds.setncattr("max_slots_D", int(self._max_slot_by_direction.get("D", 0)))
        self.ds.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Upsample merged AMSR L3 observations into a 0.5-degree regular grid "
            "while preserving duplicates using A/D slot variables."
        )
    )
    parser.add_argument(
        "input_path",
        type=Path,
        nargs="?",
        default=Path("processed/l3_daily"),
        help=(
            "Input merged NetCDF file or directory containing merged NetCDF files. "
            "(default: processed/l3_daily)"
        ),
    )
    parser.add_argument(
        "output_path",
        type=Path,
        nargs="?",
        default=Path("processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg.nc"),
        help=(
            "Output NetCDF path for 0.5-degree grid. "
            "(default: processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg.nc)"
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output file if it already exists.",
    )
    return parser.parse_args()


def list_input_files(input_path: Path) -> List[Path]:
    if input_path.is_file():
        if input_path.suffix.lower() != ".nc":
            raise ValueError(f"Input file must be .nc: {input_path}")
        return [input_path]

    files = sorted(input_path.glob("*.nc"))
    if not files:
        raise FileNotFoundError(f"No .nc files found under: {input_path}")
    return files


def collect_slot_pairs(ds: Dataset) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    for name in sorted(ds.variables):
        m = SLOT_PATTERN.match(name)
        if not m:
            continue
        direction = m.group(1)
        suffix = f"{direction}_{m.group(2)}"
        tm_name = f"observation_time_min_{suffix}"
        if tm_name in ds.variables:
            out.append((name, tm_name))
    return out


def build_day_index(
    files: List[Path],
) -> tuple[
    dict[int, List[SourceChunk]],
    Dict[Path, List[Tuple[str, str]]],
    Dict[Path, tuple[np.ndarray, np.ndarray]],
    Dict[str, np.ndarray],
    np.ndarray,
    int,
]:
    day_index: dict[int, List[SourceChunk]] = defaultdict(list)
    slot_pairs_by_file: Dict[Path, List[Tuple[str, str]]] = {}
    map_by_file: Dict[Path, tuple[np.ndarray, np.ndarray]] = {}
    file_lat_lon: Dict[str, np.ndarray] = {}

    target_res = 0.5
    nlat = int(round(180.0 / target_res))
    nlon = int(round(360.0 / target_res))

    first_time_seen = None
    for fp in files:
        with Dataset(fp) as ds:
            required = {"time", "lat", "lon"}
            if not required.issubset(ds.variables):
                continue

            slot_pairs = collect_slot_pairs(ds)
            if not slot_pairs:
                continue

            times = np.array(ds.variables["time"][:], dtype=np.int64)
            if times.ndim != 1:
                raise ValueError(f"{fp} time variable is not 1-D.")

            slot_pairs_by_file[fp] = slot_pairs
            if first_time_seen is None and times.size > 0:
                first_time_seen = times[0]

            src_lat = np.array(ds.variables["lat"][:], dtype=np.float32)
            src_lon = np.array(ds.variables["lon"][:], dtype=np.float32)
            file_lat_lon[str(fp)] = np.array(src_lat, dtype=np.float32)
            map_by_file[fp] = build_grid_index_maps(src_lat, src_lon, target_res, nlat, nlon)

            for ti, day in enumerate(times):
                day_index[int(day)].append(SourceChunk(file_path=fp, time_index=int(ti)))

    if not day_index:
        raise ValueError("No merged AMSR observation variables found in input files.")

    return day_index, slot_pairs_by_file, map_by_file, file_lat_lon, np.array([], dtype=np.float32), nlat


def build_target_grid(resolution: float) -> tuple[np.ndarray, np.ndarray]:
    nlat = int(round(180.0 / resolution))
    nlon = int(round(360.0 / resolution))
    lat = (90.0 - resolution / 2.0) - np.arange(nlat, dtype=np.float32) * resolution
    lon = (-180.0 + resolution / 2.0) + np.arange(nlon, dtype=np.float32) * resolution
    return lat, lon


def build_grid_index_maps(
    src_lat: np.ndarray,
    src_lon: np.ndarray,
    target_resolution: float,
    target_nlat: int,
    target_nlon: int,
) -> tuple[np.ndarray, np.ndarray]:
    # Source lat in merge_amsr_l3_daily.py is north->south.
    lat_idx = np.floor((90.0 - src_lat - target_resolution / 2.0) / target_resolution).astype(np.int64)
    lon_wrapped = np.mod(src_lon + 180.0, 360.0) - 180.0
    lon_idx = np.floor((lon_wrapped - (-180.0 + target_resolution / 2.0)) / target_resolution).astype(
        np.int64
    )

    lat_idx = np.clip(lat_idx, 0, target_nlat - 1).astype(np.int32)
    lon_idx = np.clip(lon_idx, 0, target_nlon - 1).astype(np.int32)
    return lat_idx, lon_idx


def allocate_slots(counts: np.ndarray, target_cells: np.ndarray) -> np.ndarray:
    if target_cells.size == 0:
        return np.empty(0, dtype=np.int32)

    order = np.argsort(target_cells, kind="mergesort")
    sorted_cells = target_cells[order]
    starts = np.zeros_like(sorted_cells, dtype=np.bool_)
    starts[0] = True
    starts[1:] = sorted_cells[1:] != sorted_cells[:-1]
    group_start = np.where(starts, np.arange(len(sorted_cells), dtype=np.int64), -1)
    last_start = np.maximum.accumulate(group_start)
    rank_in_cell = np.arange(len(sorted_cells), dtype=np.int64) - last_start

    base = counts[sorted_cells].astype(np.int64)
    slot_ids_sorted = base + 1 + rank_in_cell

    cell_unique, n_per_cell = np.unique(sorted_cells, return_counts=True)
    counts[cell_unique] += n_per_cell

    slot_ids = np.empty_like(slot_ids_sorted, dtype=np.int32)
    slot_ids[order] = slot_ids_sorted.astype(np.int32)
    return slot_ids


def parse_direction_from_sm_name(name: str) -> str:
    m = SLOT_PATTERN.match(name)
    if not m:
        raise ValueError(f"Unsupported variable name: {name}")
    return m.group(1)


def main() -> None:
    args = parse_args()

    input_files = list_input_files(args.input_path)
    day_index, slot_pairs_by_file, index_maps, file_lat_lon, _, target_nlat = build_day_index(input_files)

    target_resolution = 0.5
    target_nlon = int(round(360.0 / target_resolution))
    target_lat, target_lon = build_target_grid(target_resolution)

    writer = UpscaledNetCDFWriter(
        output_path=args.output_path,
        nlat=target_nlat,
        nlon=target_nlon,
        lat=target_lat,
        lon=target_lon,
        overwrite=args.overwrite,
    )
    writer.ds.setncattr("source_files", ",".join(str(p) for p in input_files))
    writer.ds.setncattr("source_resolution", "original merged AMSR grid")
    writer.ds.setncattr("target_resolution_deg", target_resolution)
    writer.ds.setncattr("interpolation", "nearest-neighbor assignment to target cell center")
    writer.ds.setncattr("slot_strategy", "duplicates preserved in additional slot variables")

    target_cell_count = target_nlat * target_nlon

    for day in sorted(day_index.keys()):
        writer.write_time(day)
        out_t = writer.next_time()
        counts_by_direction: DefaultDict[str, np.ndarray] = defaultdict(
            lambda: np.zeros(target_cell_count, dtype=np.int32)
        )

        for chunk in day_index[day]:
            with Dataset(chunk.file_path) as ds:
                slot_pairs = slot_pairs_by_file.get(chunk.file_path, [])
                if not slot_pairs:
                    continue
                lat_idx_src, lon_idx_src = index_maps[chunk.file_path]
                for sm_name, tm_name in slot_pairs:
                    direction = parse_direction_from_sm_name(sm_name)
                    sm = np.array(ds.variables[sm_name][chunk.time_index, :, :], dtype=np.float32)
                    tm = np.array(ds.variables[tm_name][chunk.time_index, :, :], dtype=np.float32)

                    valid = np.isfinite(sm) & np.isfinite(tm)
                    if not np.any(valid):
                        continue

                    lat_idx2d, lon_idx2d = np.nonzero(valid)
                    t_lat = lat_idx_src[lat_idx2d]
                    t_lon = lon_idx_src[lon_idx2d]
                    target_cells = (t_lat.astype(np.int64) * target_nlon) + t_lon.astype(np.int64)

                    sm_vals = sm[valid]
                    tm_vals = tm[valid]
                    slot_ids = allocate_slots(counts_by_direction[direction], target_cells)

                    for slot in np.unique(slot_ids):
                        mask = slot_ids == slot
                        var_sm, var_tm = writer.ensure_slot(direction, int(slot))
                        var_sm[out_t, t_lat[mask], t_lon[mask]] = sm_vals[mask]
                        var_tm[out_t, t_lat[mask], t_lon[mask]] = tm_vals[mask]

        if not counts_by_direction["A"].any() and not counts_by_direction["D"].any():
            pass

    writer.close()


if __name__ == "__main__":
    main()
