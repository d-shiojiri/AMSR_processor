from __future__ import annotations

import datetime as dt
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np
from netCDF4 import Dataset

from query_amsr_l3_0p5deg import AmsrUpsampledL3ObservationReader

from .io_utils import collect_slot_pairs, list_input_files, SEC_PER_DAY
from .reference import load_reference_paths, ReferenceReader
from .writer import AveragedProcessedWriter, YearlyProcessedWriter


def _to_datetime(value: dt.datetime | str) -> dt.datetime:
    if isinstance(value, dt.datetime):
        return value
    return dt.datetime.fromisoformat(value)


def _to_epoch_seconds(value: dt.datetime) -> int:
    return int((value - dt.datetime(1970, 1, 1)).total_seconds())


def _build_rank_mapper(amsr_values: np.ndarray, ref_values: np.ndarray):
    if amsr_values.size == 0:
        return None
    amsr_sorted = np.sort(np.array(amsr_values, dtype=np.float32))
    ref_sorted = np.sort(np.array(ref_values, dtype=np.float32))
    n = amsr_sorted.size
    if n == 0:
        return None

    cdf = (np.arange(1, n + 1, dtype=np.float64) - 0.5) / n
    if n == 1:
        single_val = float(ref_sorted[0])

        def mapper(values: np.ndarray) -> np.ndarray:
            return np.full(np.asarray(values, dtype=np.float32).shape, single_val, dtype=np.float32)

        return mapper

    def mapper(values: np.ndarray) -> np.ndarray:
        v = np.array(values, dtype=np.float32)
        p = np.interp(v.astype(np.float64), amsr_sorted.astype(np.float64), cdf, left=0.0, right=1.0)
        return np.interp(p, cdf, ref_sorted.astype(np.float64), left=float(ref_sorted[0]), right=float(ref_sorted[-1])).astype(
            np.float32
        )

    return mapper


def _auto_range_days(reader: AmsrUpsampledL3ObservationReader) -> tuple[dt.datetime, dt.datetime]:
    if not reader.day_index:
        raise ValueError("No AMSR records found in input.")
    day_keys = sorted(reader.day_index.keys())
    start_sec = int(day_keys[0]) * SEC_PER_DAY
    end_sec = int(day_keys[-1] + 1) * SEC_PER_DAY - 1
    return (
        dt.datetime.utcfromtimestamp(start_sec),
        dt.datetime.utcfromtimestamp(end_sec),
    )


def _is_averaged_input(ds: Dataset) -> bool:
    return (
        "soil_moisture" in ds.variables
        and "observation_count" in ds.variables
        and "time" in ds.variables
        and not collect_slot_pairs(ds)
    )


@dataclass
class CDFMatchConfig:
    input_path: Path
    reference: list[str] | None
    reference_dict: str | None
    reference_time: str | None
    reference_lat: str | None
    reference_lon: str | None
    reference_var: str | None
    start: str | None
    end: str | None
    step_hours: float
    window_hours: float
    output_dir: Path
    output_mode: str
    output_file: Path
    overwrite: bool
    reference_cache_size: int
    reference_depth: int = 0


def _within_time_range(day_value: int, start_sec: int | None, end_sec: int | None) -> bool:
    if start_sec is None and end_sec is None:
        return True
    day_start = int(day_value) * SEC_PER_DAY
    day_end = day_start + SEC_PER_DAY - 1
    if start_sec is not None and day_end < start_sec:
        return False
    if end_sec is not None and day_start > end_sec:
        return False
    return True


def run_cdf_matching(config: CDFMatchConfig) -> None:
    input_path = config.input_path
    if not input_path.exists():
        raise FileNotFoundError(f"Input path not found: {input_path}")

    reference_paths = load_reference_paths(config.reference, config.reference_dict)
    for p in reference_paths:
        if not p.exists():
            raise FileNotFoundError(f"Reference file not found: {p}")

    input_files = list_input_files(input_path)
    reader = AmsrUpsampledL3ObservationReader(
        input_path,
        interval_hours=config.step_hours,
        window_hours=config.window_hours,
    )

    start = _to_datetime(config.start) if config.start else None
    end = _to_datetime(config.end) if config.end else None
    if start is None or end is None:
        auto_start, auto_end = _auto_range_days(reader)
        if start is None:
            start = auto_start
        if end is None:
            end = auto_end
    if end < start:
        raise ValueError("end must be >= start")

    start_sec = _to_epoch_seconds(start)
    end_sec = _to_epoch_seconds(end)
    step = dt.timedelta(hours=config.step_hours)
    window = dt.timedelta(hours=config.window_hours)

    ref_reader = ReferenceReader(
        paths=reference_paths,
        time_name=config.reference_time,
        lat_name=config.reference_lat,
        lon_name=config.reference_lon,
        ref_var_name=config.reference_var,
        reference_depth=config.reference_depth,
        cache_size=config.reference_cache_size,
    )

    # pass 1: build global CDF mapping from paired AMSR/REF values
    amsr_values: list[float] = []
    ref_values: list[float] = []
    for _, _, sm, obs_time, lat, lon in reader.iterate_windows(start, end, step=step, window=window):
        if sm.size == 0:
            continue
        obs_sec = obs_time.astype("datetime64[s]").astype(np.int64)
        ref_at_obs = ref_reader.sample(obs_sec, lat, lon, require_within_reference=True)
        valid = np.isfinite(sm) & np.isfinite(ref_at_obs)
        if not np.any(valid):
            continue
        amsr_values.append(np.asarray(sm[valid], dtype=np.float32).tolist())
        ref_values.append(np.asarray(ref_at_obs[valid], dtype=np.float32).tolist())

    amsr_values_arr = np.array(np.concatenate(amsr_values), dtype=np.float32) if amsr_values else np.empty(0, dtype=np.float32)
    ref_values_arr = np.array(np.concatenate(ref_values), dtype=np.float32) if ref_values else np.empty(0, dtype=np.float32)
    if amsr_values_arr.size == 0:
        ref_reader.close()
        raise ValueError("No valid AMSR-reference pairs in the given range.")

    map_fn = _build_rank_mapper(amsr_values_arr, ref_values_arr)
    if map_fn is None:
        ref_reader.close()
        raise ValueError("CDF map creation failed.")

    # output writer (same shape/variables as input layout)
    with Dataset(input_files[0], mode="r") as probe:
        lat = np.array(probe.variables["lat"][:], dtype=np.float32)
        lon = np.array(probe.variables["lon"][:], dtype=np.float32)
        time_units = getattr(probe.variables["time"], "units", "days since 1970-01-01")
        has_slot_pairs = bool(collect_slot_pairs(probe))
        is_averaged = _is_averaged_input(probe)
        if not has_slot_pairs and not is_averaged:
            raise ValueError("Input NetCDF does not match supported formats (slot-based or 0.5-degree averaged).")

    if is_averaged:
        writer = AveragedProcessedWriter(
            output_dir=config.output_dir,
            output_file=config.output_file,
            mode=config.output_mode,
            overwrite=config.overwrite,
            nlat=lat.size,
            nlon=lon.size,
            lat=lat,
            lon=lon,
            time_units=time_units,
        )
    else:
        writer = YearlyProcessedWriter(
            output_dir=config.output_dir,
            output_file=config.output_file,
            mode=config.output_mode,
            overwrite=config.overwrite,
            nlat=lat.size,
            nlon=lon.size,
            lat=lat,
            lon=lon,
            time_units=time_units,
        )

    # pass 2: re-read input and overwrite AMSR values by CDF ranking
    for path in input_files:
        with Dataset(path, mode="r") as ds:
            times = np.array(ds.variables["time"][:], dtype=np.int64)
            if times.ndim != 1:
                raise ValueError(f"{path} time variable must be 1D.")
            lat_src = np.array(ds.variables["lat"][:], dtype=np.float32)
            lon_src = np.array(ds.variables["lon"][:], dtype=np.float32)
            if is_averaged:
                if "soil_moisture" not in ds.variables or "observation_count" not in ds.variables:
                    continue
                sm_var = ds.variables["soil_moisture"]
                count_var = ds.variables["observation_count"]

                for ti, day in enumerate(times):
                    if not _within_time_range(int(day), start_sec, end_sec):
                        continue
                    day_val = int(day)
                    day_origin = dt.datetime(1970, 1, 1) + dt.timedelta(seconds=day_val * SEC_PER_DAY)
                    year = day_origin.year

                    sm_2d = np.array(sm_var[ti, :, :], dtype=np.float32)
                    count_2d = np.array(count_var[ti, :, :], dtype=np.int32)
                    sm_out = np.array(sm_2d, copy=True)

                    valid = np.isfinite(sm_out) & (count_2d > 0)
                    if np.any(valid):
                        lat_idx, lon_idx = np.nonzero(valid)
                        sm_out[lat_idx, lon_idx] = map_fn(sm_out[lat_idx, lon_idx])

                    writer.write_day(day_val, year, sm_out.astype(np.float32), count_2d.astype(np.int32))
            else:
                slot_pairs = collect_slot_pairs(ds)
                if not slot_pairs:
                    continue

                for ti, day in enumerate(times):
                    if not _within_time_range(int(day), start_sec, end_sec):
                        continue
                    day_val = int(day)
                    day_origin = dt.datetime(1970, 1, 1) + dt.timedelta(seconds=day_val * SEC_PER_DAY)
                    year = day_origin.year

                    payload: list[Tuple[str, str, np.ndarray, np.ndarray]] = []
                    for sm_name, tm_name in slot_pairs:
                        sm_2d = np.array(ds.variables[sm_name][ti, :, :], dtype=np.float32)
                        tm_2d = np.array(ds.variables[tm_name][ti, :, :], dtype=np.float32)
                        sm_out = np.array(sm_2d, copy=True)

                        valid = np.isfinite(sm_out) & np.isfinite(tm_2d)
                        if np.any(valid):
                            lat_idx, lon_idx = np.nonzero(valid)
                            sm_out[lat_idx, lon_idx] = map_fn(sm_out[lat_idx, lon_idx])
                        payload.append((sm_name, tm_name, sm_out.astype(np.float32), tm_2d.astype(np.float32)))

                    writer.write_day(day_val, year, payload, ds)

    writer.close()
    ref_reader.close()
