from __future__ import annotations

import datetime as dt
import warnings
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from tqdm.auto import tqdm

from .io_utils import SEC_PER_DAY, collect_slot_pairs, input_to_m3m3, list_input_files
from .reference import ReferenceReader, load_reference_paths
from .writer import create_output_like


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
    output_dir: Path
    output_file: Path
    overwrite: bool
    reference_time_mode: str
    reference_depth: int
    block_rows: int
    map_columns: int
    workers: int


def _to_datetime(value: str | None) -> dt.datetime | None:
    if value is None:
        return None
    return dt.datetime.fromisoformat(value)


def _day_to_epoch_seconds(day_value: np.ndarray | int) -> np.ndarray:
    return np.asarray(day_value, dtype=np.int64) * np.int64(SEC_PER_DAY)


def _within_range(days: np.ndarray, start: dt.datetime | None, end: dt.datetime | None) -> np.ndarray:
    mask = np.ones(days.shape, dtype=bool)
    if start is not None:
        start_day = int((start - dt.datetime(1970, 1, 1)).total_seconds()) // SEC_PER_DAY
        mask &= days >= start_day
    if end is not None:
        end_day = int((end - dt.datetime(1970, 1, 1)).total_seconds()) // SEC_PER_DAY
        mask &= days <= end_day
    return mask


def _seconds_to_datetime(seconds: int) -> dt.datetime:
    return dt.datetime(1970, 1, 1) + dt.timedelta(seconds=int(seconds))


def _rank_map_one(raw: np.ndarray, ref: np.ndarray) -> np.ndarray:
    valid = np.isfinite(raw) & np.isfinite(ref)
    out = np.full(raw.shape, np.float32(np.nan), dtype=np.float32)
    if valid.sum() == 0:
        return out
    if valid.sum() == 1:
        out[np.isfinite(raw)] = np.float32(ref[valid][0])
        return out

    raw_valid = raw[valid].astype(np.float64)
    ref_valid = ref[valid].astype(np.float64)
    order_raw = np.argsort(raw_valid)
    order_ref = np.argsort(ref_valid)
    raw_sorted = raw_valid[order_raw]
    ref_sorted = ref_valid[order_ref]
    n = raw_sorted.size
    cdf = (np.arange(1, n + 1, dtype=np.float64) - 0.5) / n

    valid_pos = np.flatnonzero(valid)
    out[valid_pos[order_raw]] = ref_sorted.astype(np.float32)

    # If a raw sample is finite but has no paired reference value, keep it
    # usable by interpolating from the empirical paired CDF.
    unpaired_raw = np.isfinite(raw) & ~valid
    if np.any(unpaired_raw):
        p = np.interp(raw[unpaired_raw].astype(np.float64), raw_sorted, cdf, left=0.0, right=1.0)
        out[unpaired_raw] = np.interp(p, cdf, ref_sorted, left=ref_sorted[0], right=ref_sorted[-1]).astype(np.float32)
    return out


def _map_flat_range(raw_flat: np.ndarray, ref_flat: np.ndarray, start: int, stop: int) -> tuple[int, np.ndarray]:
    mapped = np.empty((stop - start, raw_flat.shape[1]), dtype=np.float32)
    for local, grid in enumerate(range(start, stop)):
        mapped[local] = _rank_map_one(raw_flat[grid], ref_flat[grid])
    return start, mapped


def _gridwise_map_columns(raw_cols: np.ndarray, ref_cols: np.ndarray, workers: int) -> np.ndarray:
    out_cols = np.full(raw_cols.shape, np.float32(np.nan), dtype=np.float32)

    valid = np.isfinite(raw_cols) & np.isfinite(ref_cols)
    full_valid = valid.all(axis=0)
    if np.any(full_valid):
        raw_full = raw_cols[:, full_valid]
        ref_full = ref_cols[:, full_valid]
        raw_order = np.argsort(raw_full, axis=0)
        ref_sorted = np.sort(ref_full, axis=0).astype(np.float32, copy=False)
        mapped_full = np.empty(raw_full.shape, dtype=np.float32)
        np.put_along_axis(mapped_full, raw_order, ref_sorted, axis=0)
        out_cols[:, full_valid] = mapped_full

    partial_valid = valid.any(axis=0) & ~full_valid
    if not np.any(partial_valid):
        return out_cols

    partial_indices = np.flatnonzero(partial_valid)
    raw_flat = raw_cols[:, partial_indices].T
    ref_flat = ref_cols[:, partial_indices].T
    out_flat = np.empty_like(raw_flat, dtype=np.float32)
    ngrid = raw_flat.shape[0]
    workers = max(1, int(workers))
    if workers == 1 or ngrid < 256:
        out_flat[:] = _map_flat_range(raw_flat, ref_flat, 0, ngrid)[1]
    else:
        chunk = max(1, (ngrid + workers - 1) // workers)
        ranges = [(i, min(i + chunk, ngrid)) for i in range(0, ngrid, chunk)]
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for start, mapped in ex.map(lambda r: _map_flat_range(raw_flat, ref_flat, r[0], r[1]), ranges):
                out_flat[start:start + mapped.shape[0]] = mapped
    out_cols[:, partial_indices] = out_flat.T
    return out_cols


def _gridwise_map(raw_stack: np.ndarray, ref_stack: np.ndarray, workers: int, map_columns: int) -> np.ndarray:
    # stacks are (sample, row, lon). Map independent grid cells.
    sample, nrow, nlon = raw_stack.shape
    raw_cols = raw_stack.reshape(sample, nrow * nlon)
    ref_cols = ref_stack.reshape(sample, nrow * nlon)
    out_cols = np.empty(raw_cols.shape, dtype=np.float32)

    map_columns = max(1, int(map_columns))
    for col0 in range(0, raw_cols.shape[1], map_columns):
        col1 = min(col0 + map_columns, raw_cols.shape[1])
        out_cols[:, col0:col1] = _gridwise_map_columns(raw_cols[:, col0:col1], ref_cols[:, col0:col1], workers)
    return out_cols.reshape(sample, nrow, nlon)


def _reference_for_samples(
    ref_reader: ReferenceReader,
    input_lat: np.ndarray,
    input_lon: np.ndarray,
    row_slice: slice,
    day_values: np.ndarray,
    obs_min: np.ndarray,
) -> np.ndarray:
    valid = np.isfinite(obs_min)
    obs_seconds = np.zeros(obs_min.shape, dtype=np.int64)
    if valid.any():
        day_seconds = _day_to_epoch_seconds(day_values)[:, None, None]
        obs_seconds[valid] = (day_seconds + np.rint(np.where(valid, obs_min, 0.0) * 60.0).astype(np.int64))[valid]
    flat_obs = obs_seconds[valid]
    out = np.full(obs_min.shape, np.float32(np.nan), dtype=np.float32)
    if flat_obs.size == 0:
        return out

    source_ids, local_ids = ref_reader.nearest_lookup(flat_obs)
    flat_out = np.full(flat_obs.shape, np.float32(np.nan), dtype=np.float32)
    flat_pos = np.flatnonzero(valid.ravel())
    block_shape = obs_min.shape

    for source_id in np.unique(source_ids):
        src = ref_reader.sources[int(source_id)]
        pos_source = np.nonzero(source_ids == source_id)[0]
        unique_local_ids = np.unique(local_ids[pos_source])
        grids = src.read_block_times(unique_local_ids, input_lat, input_lon, row_slice)
        for grid_index, local_id in enumerate(unique_local_ids):
            pos = pos_source[local_ids[pos_source] == local_id]
            sample_idx, row_idx, lon_idx = np.unravel_index(flat_pos[pos], block_shape)
            del sample_idx
            flat_out[pos] = grids[grid_index, row_idx, lon_idx]
    out[valid] = flat_out
    return out


def _reference_for_slot_median_samples(
    ref_reader: ReferenceReader,
    input_lat: np.ndarray,
    input_lon: np.ndarray,
    row_slice: slice,
    day_values: np.ndarray,
    obs_min: np.ndarray,
) -> np.ndarray:
    valid = np.isfinite(obs_min)
    out = np.full(obs_min.shape, np.float32(np.nan), dtype=np.float32)
    sample_has_obs = valid.reshape(valid.shape[0], -1).any(axis=1)
    if not np.any(sample_has_obs):
        return out

    representative_min = np.full(obs_min.shape[0], np.float32(np.nan), dtype=np.float32)
    for sample_idx in np.flatnonzero(sample_has_obs):
        representative_min[sample_idx] = np.nanmedian(obs_min[sample_idx])

    sample_indices = np.flatnonzero(np.isfinite(representative_min))
    obs_seconds = (
        _day_to_epoch_seconds(day_values[sample_indices])
        + np.rint(representative_min[sample_indices].astype(np.float64) * 60.0).astype(np.int64)
    )
    source_ids, local_ids = ref_reader.nearest_lookup(obs_seconds)

    for source_id in np.unique(source_ids):
        src = ref_reader.sources[int(source_id)]
        pos_source = np.nonzero(source_ids == source_id)[0]
        unique_local_ids = np.unique(local_ids[pos_source])
        grids = src.read_block_times(unique_local_ids, input_lat, input_lon, row_slice)
        for grid_index, local_id in enumerate(unique_local_ids):
            pos = pos_source[local_ids[pos_source] == local_id]
            for pos_idx in pos:
                sample_idx = sample_indices[pos_idx]
                mask = valid[sample_idx]
                out[sample_idx][mask] = grids[grid_index][mask]
    return out


def _references_for_slot_median_samples(
    ref_reader: ReferenceReader,
    input_lat: np.ndarray,
    input_lon: np.ndarray,
    row_slice: slice,
    day_values: np.ndarray,
    obs_min_list: list[np.ndarray],
) -> list[np.ndarray]:
    outputs = [np.full(obs_min.shape, np.float32(np.nan), dtype=np.float32) for obs_min in obs_min_list]
    requests: list[tuple[int, int, np.ndarray]] = []
    obs_seconds_all: list[int] = []

    for slot_idx, obs_min in enumerate(obs_min_list):
        valid = np.isfinite(obs_min)
        sample_has_obs = valid.reshape(valid.shape[0], -1).any(axis=1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            representative_min = np.nanmedian(obs_min, axis=(1, 2)).astype(np.float32)
        representative_min[~sample_has_obs] = np.float32(np.nan)
        sample_indices = np.flatnonzero(np.isfinite(representative_min))
        if sample_indices.size == 0:
            continue
        obs_seconds = (
            _day_to_epoch_seconds(day_values[sample_indices])
            + np.rint(representative_min[sample_indices].astype(np.float64) * 60.0).astype(np.int64)
        )
        for sample_idx, obs_second in zip(sample_indices, obs_seconds):
            requests.append((slot_idx, int(sample_idx), valid[int(sample_idx)]))
            obs_seconds_all.append(int(obs_second))

    if not requests:
        return outputs

    source_ids, local_ids = ref_reader.nearest_lookup(np.asarray(obs_seconds_all, dtype=np.int64))
    for source_id in np.unique(source_ids):
        src = ref_reader.sources[int(source_id)]
        pos_source = np.nonzero(source_ids == source_id)[0]
        unique_local_ids = np.unique(local_ids[pos_source])
        grids = src.read_block_times(unique_local_ids, input_lat, input_lon, row_slice)
        local_to_grid = {int(local_id): grids[i] for i, local_id in enumerate(unique_local_ids)}
        for pos in pos_source:
            slot_idx, sample_idx, mask = requests[int(pos)]
            grid = local_to_grid[int(local_ids[pos])]
            outputs[slot_idx][sample_idx][mask] = grid[mask]
    return outputs


def run_cdf_matching(config: CDFMatchConfig) -> None:
    input_files = list_input_files(config.input_path)
    if len(input_files) != 1:
        raise NotImplementedError("Grid-wise CDF matching currently expects one 0.5-degree AMSR NetCDF file.")

    reference_paths = load_reference_paths(config.reference, config.reference_dict)
    for path in reference_paths:
        if not path.exists():
            raise FileNotFoundError(f"Reference file not found: {path}")

    ref_reader = ReferenceReader(
        reference_paths,
        time_name=config.reference_time,
        lat_name=config.reference_lat,
        lon_name=config.reference_lon,
        ref_var_name=config.reference_var,
        reference_depth=config.reference_depth,
    )
    reference_start = _seconds_to_datetime(int(ref_reader.time_seconds.min()))
    reference_end = _seconds_to_datetime(int(ref_reader.time_seconds.max()))
    start = _to_datetime(config.start) or reference_start
    end = _to_datetime(config.end) or reference_end

    input_path = input_files[0]
    output_path = config.output_dir / config.output_file
    try:
        with Dataset(input_path, "r") as src:
            slot_pairs = collect_slot_pairs(src)
            if not slot_pairs:
                raise ValueError(
                    "Input file must contain soil_moisture/observation_time_min "
                    "or slot variables like soil_moisture_A_01."
                )

            times = np.asarray(src.variables["time"][:], dtype=np.int64)
            time_mask = _within_range(times, start, end)
            selected = np.flatnonzero(time_mask)
            if selected.size == 0:
                raise ValueError("No input days are inside the requested time range.")
            tqdm.write(
                "CDF matching period: "
                f"{_seconds_to_datetime(int(times[selected[0]]) * SEC_PER_DAY).date()} "
                f"to {_seconds_to_datetime(int(times[selected[-1]]) * SEC_PER_DAY).date()}"
            )

            lat = np.asarray(src.variables["lat"][:], dtype=np.float64)
            lon = np.asarray(src.variables["lon"][:], dtype=np.float64)
            nlat = lat.size

            with create_output_like(output_path, src, slot_pairs, config.overwrite) as dst:
                block_rows = max(1, int(config.block_rows))
                pbar = tqdm(range(0, nlat, block_rows), desc="grid-wise CDF", unit="block")
                for row0 in pbar:
                    row1 = min(row0 + block_rows, nlat)
                    row_slice = slice(row0, row1)

                    raw_samples: list[np.ndarray] = []
                    ref_samples: list[np.ndarray] = []
                    obs_min_samples: list[np.ndarray] = []
                    sample_targets: list[tuple[str, np.ndarray]] = []

                    for sm_name, tm_name in slot_pairs:
                        sm_var = src.variables[sm_name]
                        tm_var = src.variables[tm_name]
                        units = getattr(sm_var, "units", None)

                        sm = input_to_m3m3(sm_var[selected, row_slice, :], units=units)
                        tm = np.asarray(tm_var[selected, row_slice, :], dtype=np.float32)
                        raw_samples.append(sm)
                        obs_min_samples.append(tm)
                        sample_targets.append((sm_name, selected))

                    if config.reference_time_mode == "pixel":
                        ref_samples = [
                            _reference_for_samples(ref_reader, lat, lon, row_slice, times[selected], tm)
                            for tm in obs_min_samples
                        ]
                    elif config.reference_time_mode == "slot-median":
                        ref_samples = _references_for_slot_median_samples(
                            ref_reader, lat, lon, row_slice, times[selected], obs_min_samples
                        )
                    else:
                        raise ValueError(f"Unsupported reference_time_mode: {config.reference_time_mode}")

                    raw_stack = np.concatenate(raw_samples, axis=0)
                    ref_stack = np.concatenate(ref_samples, axis=0)
                    mapped_stack = _gridwise_map(raw_stack, ref_stack, config.workers, config.map_columns)

                    offset = 0
                    for (sm_name, selected_idx), sm in zip(sample_targets, raw_samples):
                        n = sm.shape[0]
                        dst.variables[sm_name][selected_idx, row_slice, :] = mapped_stack[offset:offset + n]
                        offset += n
    finally:
        ref_reader.close()
