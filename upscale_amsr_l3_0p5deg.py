#!/usr/bin/env python3
"""
Upscale merged AMSR-L3 NetCDF observations to a 0.5-degree grid by averaging.

For each day and each 0.5-degree cell, all observations that map into the cell are
averaged. No averaging is performed across day boundaries.
"""

from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, List, Tuple

import numpy as np
from netCDF4 import Dataset

SLOT_PATTERN = re.compile(r"^soil_moisture_([AD])_(\d{2})$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Upscale merged AMSR-L3 NetCDF observations to 0.5-degree grid using daily "
            "cell averages."
        )
    )
    parser.add_argument(
        "input_path",
        type=Path,
        nargs="?",
        default=Path("processed/l3_daily"),
        help=(
            "Input merged NetCDF file or directory containing merged NetCDF files "
            "(default: processed/l3_daily)"
        ),
    )
    parser.add_argument(
        "output_path",
        type=Path,
        nargs="?",
        default=Path("processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc"),
        help=(
            "Output NetCDF path for 0.5-degree averaged grid "
            "(default: processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc)"
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output file if it already exists.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=min(4, os.cpu_count() or 1),
        help="Number of parallel workers for per-day averaging.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable tqdm progress bar.",
    )
    return parser.parse_args()


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

        output_path.parent.mkdir(parents=True, exist_ok=True)
        self.ds = Dataset(output_path, "w", format="NETCDF4")
        self.ds.createDimension("time", None)
        self.ds.createDimension("lat", nlat)
        self.ds.createDimension("lon", nlon)

        self.var_time = self.ds.createVariable("time", "i4", ("time",))
        var_lat = self.ds.createVariable("lat", "f4", ("lat",))
        var_lon = self.ds.createVariable("lon", "f4", ("lon",))

        self.var_soil = self.ds.createVariable(
            "soil_moisture",
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, nlat), min(360, nlon)),
            fill_value=np.float32(np.nan),
        )
        self.var_count = self.ds.createVariable(
            "observation_count",
            "i4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, nlat), min(360, nlon)),
            fill_value=np.int32(0),
        )
        self.var_obs_time_min = self.ds.createVariable(
            "observation_time_min",
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, nlat), min(360, nlon)),
            fill_value=np.float32(np.nan),
        )

        var_lat[:] = lat.astype(np.float32)
        var_lon[:] = lon.astype(np.float32)

        self.var_time.units = "days since 1970-01-01"
        self.var_time.long_name = "day timestamp"
        var_lat.units = "degrees_north"
        var_lon.units = "degrees_east"
        self.var_soil.units = "%"
        self.var_soil.long_name = "daily mean soil moisture (0.5 degree)"
        self.var_soil.comment = "Mean over all observations that map to the 0.5-degree cell for the day."
        self.var_count.units = "1"
        self.var_count.long_name = "number of merged observations used in daily average"
        self.var_count.comment = "Duplicates preserved by raw accumulation before averaging."
        self.var_obs_time_min.units = "minutes since 00:00 UTC of corresponding time day"
        self.var_obs_time_min.long_name = "daily mean observation time (0.5 degree)"
        self.var_obs_time_min.comment = (
            "Mean observation_time_min over all observations used in the daily cell average."
        )

        self.time_index = 0

    def write_time(self, day_value: int) -> None:
        self.var_time[self.time_index] = int(day_value)

    def next_time(self) -> int:
        cur = self.time_index
        self.time_index += 1
        return cur

    def write_day(
        self,
        day_index: int,
        soil_sum_flat: np.ndarray,
        time_min_sum_flat: np.ndarray,
        count_flat: np.ndarray,
        nlat: int,
        nlon: int,
    ) -> None:
        self.write_time(day_index)
        t = self.next_time()

        soil_mean_grid = np.full((nlat * nlon,), np.nan, dtype=np.float32)
        time_mean_grid = np.full((nlat * nlon,), np.nan, dtype=np.float32)
        idx = count_flat > 0
        if np.any(idx):
            soil_mean_grid[idx] = (soil_sum_flat[idx] / count_flat[idx]).astype(np.float32)
            time_mean_grid[idx] = (time_min_sum_flat[idx] / count_flat[idx]).astype(np.float32)

        self.var_soil[t, :, :] = soil_mean_grid.reshape((nlat, nlon))
        self.var_count[t, :, :] = count_flat.reshape((nlat, nlon)).astype(np.int32)
        self.var_obs_time_min[t, :, :] = time_mean_grid.reshape((nlat, nlon))

    def close(self) -> None:
        self.ds.close()


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
        suffix = m.group(1) + "_" + m.group(2)
        tm_name = f"observation_time_min_{suffix}"
        if tm_name in ds.variables:
            out.append((name, tm_name))
    return out


def build_day_index(
    files: List[Path],
) -> tuple[
    Dict[int, List[SourceChunk]],
    Dict[Path, List[Tuple[str, str]]],
    int,
    int,
]:
    day_index: Dict[int, List[SourceChunk]] = defaultdict(list)
    slot_pairs_by_file: Dict[Path, List[Tuple[str, str]]] = {}

    target_res = 0.5
    nlat = int(round(180.0 / target_res))
    nlon = int(round(360.0 / target_res))

    for fp in files:
        with Dataset(fp) as ds:
            required = {"time", "lat", "lon"}
            if not required.issubset(ds.variables):
                continue

            slot_pairs = collect_slot_pairs(ds)
            if not slot_pairs:
                continue
            slot_pairs_by_file[fp] = slot_pairs

            times = np.array(ds.variables["time"][:], dtype=np.int64)
            if times.ndim != 1:
                raise ValueError(f"{fp} time variable is not 1-D.")

            for ti, day in enumerate(times):
                day_index[int(day)].append(SourceChunk(file_path=fp, time_index=int(ti)))

    if not day_index:
        raise ValueError("No merged AMSR observation variables found in input files.")

    return day_index, slot_pairs_by_file, nlat, nlon


def build_target_grid(resolution: float, nlat: int, nlon: int) -> tuple[np.ndarray, np.ndarray]:
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
    lat_idx = np.floor((90.0 - src_lat - target_resolution / 2.0) / target_resolution).astype(np.int64)
    lon_wrapped = np.mod(src_lon + 180.0, 360.0) - 180.0
    lon_idx = np.floor((lon_wrapped - (-180.0 + target_resolution / 2.0)) / target_resolution).astype(np.int64)

    lat_idx = np.clip(lat_idx, 0, target_nlat - 1).astype(np.int32)
    lon_idx = np.clip(lon_idx, 0, target_nlon - 1).astype(np.int32)
    return lat_idx, lon_idx


def process_day_payload(args: tuple[int, List[SourceChunk], Dict[Path, List[Tuple[str, str]]], int, int]) -> tuple[
    int,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """
    Aggregate one day's observations into sum and count arrays on 0.5-degree grid.
    """
    day, chunks, slot_pairs_by_file, target_nlat, target_nlon = args
    ncell = target_nlat * target_nlon
    target_soil_sum = np.zeros(ncell, dtype=np.float64)
    target_time_min_sum = np.zeros(ncell, dtype=np.float64)
    target_count = np.zeros(ncell, dtype=np.int32)

    chunks_by_file: Dict[Path, List[int]] = defaultdict(list)
    for chunk in chunks:
        chunks_by_file[chunk.file_path].append(chunk.time_index)

    cache: Dict[Path, tuple[np.ndarray, np.ndarray]] = {}
    for fp, time_indices in chunks_by_file.items():
        with Dataset(fp) as ds:
            slot_pairs = slot_pairs_by_file.get(fp, [])
            if not slot_pairs:
                continue

            # Source->target index maps are source-grid dependent.
            if fp in cache:
                lat_idx_src, lon_idx_src = cache[fp]
            else:
                src_lat = np.array(ds.variables["lat"][:], dtype=np.float32)
                src_lon = np.array(ds.variables["lon"][:], dtype=np.float32)
                lat_idx_src, lon_idx_src = build_grid_index_maps(src_lat, src_lon, 0.5, target_nlat, target_nlon)
                cache[fp] = (lat_idx_src, lon_idx_src)

            for t_idx in time_indices:
                for sm_name, tm_name in slot_pairs:
                    sm_2d = np.array(ds.variables[sm_name][t_idx, :, :], dtype=np.float32)
                    tm_2d = np.array(ds.variables[tm_name][t_idx, :, :], dtype=np.float32)
                    valid = np.isfinite(sm_2d) & np.isfinite(tm_2d)
                    if not np.any(valid):
                        continue

                    lat_idx2d, lon_idx2d = np.nonzero(valid)
                    target_lat = lat_idx_src[lat_idx2d]
                    target_lon = lon_idx_src[lon_idx2d]
                    target_cells = target_lat.astype(np.int64) * target_nlon + target_lon.astype(np.int64)

                    sm_vals = sm_2d[valid].astype(np.float64, copy=False)
                    tm_vals = tm_2d[valid].astype(np.float64, copy=False)
                    target_soil_sum += np.bincount(target_cells, weights=sm_vals, minlength=ncell)
                    target_time_min_sum += np.bincount(target_cells, weights=tm_vals, minlength=ncell)
                    target_count += np.bincount(target_cells, minlength=ncell)

    return day, target_soil_sum, target_time_min_sum, target_count


def _iterable_progress(iterable, total: int, show: bool, label: str):
    if not show:
        for item in iterable:
            yield item
        return

    try:
        from tqdm.auto import tqdm

        yield from tqdm(iterable, total=total, desc=label, unit="day")
        return
    except Exception:
        pass

    width = 30
    if total <= 0:
        for i, item in enumerate(iterable, start=1):
            yield item
        return

    for i, item in enumerate(iterable, start=1):
        ratio = i / total
        filled = int(ratio * width)
        bar = "#" * filled + "-" * (width - filled)
        percent = int(ratio * 100)
        print(f"\r{label}: [{bar}] {percent:3d}% ({i}/{total})", end="", flush=True)
        yield item
    print()


def main() -> None:
    args = parse_args()

    workers = max(1, args.workers)
    input_files = list_input_files(args.input_path)
    day_index, slot_pairs_by_file, target_nlat, target_nlon = build_day_index(input_files)
    days_sorted = sorted(day_index.keys())

    target_resolution = 0.5
    target_lat, target_lon = build_target_grid(target_resolution, target_nlat, target_nlon)

    tasks: list[tuple[int, List[SourceChunk], Dict[Path, List[Tuple[str, str]]], int, int]] = []
    for day in days_sorted:
        tasks.append((day, day_index[day], slot_pairs_by_file, target_nlat, target_nlon))

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
    writer.ds.setncattr("aggregation", "daily mean after mapping to nearest 0.5-degree cell")

    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for day, target_soil_sum, target_time_min_sum, target_count in _iterable_progress(
                ex.map(process_day_payload, tasks),
                total=len(tasks),
                show=not args.no_progress,
                label=f"Upscaling (parallel workers={workers})",
            ):
                writer.write_day(
                    day,
                    target_soil_sum,
                    target_time_min_sum,
                    target_count,
                    target_nlat,
                    target_nlon,
                )
    else:
        for day, target_soil_sum, target_time_min_sum, target_count in _iterable_progress(
            (process_day_payload(task) for task in tasks),
            total=len(tasks),
            show=not args.no_progress,
            label="Upscaling (single)",
        ):
            writer.write_day(
                day,
                target_soil_sum,
                target_time_min_sum,
                target_count,
                target_nlat,
                target_nlon,
            )

    writer.close()


if __name__ == "__main__":
    main()
