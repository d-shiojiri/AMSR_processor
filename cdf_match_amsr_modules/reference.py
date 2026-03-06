from __future__ import annotations

import datetime as dt
import json
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from netCDF4 import Dataset, num2date

from .io_utils import nearest_index_1d, normalize_reference_lon

EPOCH = dt.datetime(1970, 1, 1)
SEC_PER_DAY = 24 * 60 * 60


def _to_datetime(value: dt.datetime | str) -> dt.datetime:
    if isinstance(value, dt.datetime):
        return value
    return dt.datetime.fromisoformat(value)


def _epoch_seconds(value: dt.datetime | np.datetime64 | object) -> int:
    if isinstance(value, np.datetime64):
        return int(value.astype("datetime64[s]").astype(np.int64))
    if isinstance(value, dt.datetime):
        return int((value - EPOCH).total_seconds())
    if hasattr(value, "year") and hasattr(value, "month") and hasattr(value, "day"):
        hour = int(getattr(value, "hour", 0))
        minute = int(getattr(value, "minute", 0))
        second = int(getattr(value, "second", 0))
        microsecond = int(getattr(value, "microsecond", 0))
        return int(
            (dt.datetime(int(value.year), int(value.month), int(value.day), hour, minute, second, microsecond) - EPOCH).total_seconds()
        )
    raise TypeError(f"Unsupported datetime type: {type(value)}")


def _infer_time_lat_lon(
    ds: Dataset,
    time_name: str | None,
    lat_name: str | None,
    lon_name: str | None,
) -> tuple[str, str, str]:
    tn = time_name or ("time" if "time" in ds.variables else None)
    ln = lat_name or ("lat" if "lat" in ds.variables else ("latitude" if "latitude" in ds.variables else None))
    lonn = lon_name or ("lon" if "lon" in ds.variables else ("longitude" if "longitude" in ds.variables else None))
    if tn is None or ln is None or lonn is None:
        raise ValueError(f"Cannot infer coordinate variable names in {ds.filepath()}")
    if tn not in ds.variables or ln not in ds.variables or lonn not in ds.variables:
        raise ValueError(f"Missing required variables {tn}, {ln}, {lonn} in {ds.filepath()}")
    return tn, ln, lonn


def _infer_reference_variable(
    ds: Dataset,
    ref_var: str | None,
    time_name: str,
    lat_name: str,
    lon_name: str,
) -> tuple[str, str | None, int]:
    def analyze(name: str) -> tuple[bool, int | None]:
        var = ds.variables[name]
        dims = tuple(var.dimensions)
        if len(dims) == 3:
            if dims == (time_name, lat_name, lon_name):
                return True, -1
            return False, -1

        if len(dims) != 4:
            return False, -1

        required = {time_name, lat_name, lon_name}
        if not set(required).issubset(dims):
            return False, -1

        depth_axes = [i for i, d in enumerate(dims) if d not in required]
        if len(depth_axes) != 1:
            return False, -1
        return True, depth_axes[0]

    if ref_var is not None:
        if ref_var not in ds.variables:
            raise ValueError(f"Reference variable '{ref_var}' not found in {ds.filepath()}")
        ok, depth_axis = analyze(ref_var)
        if not ok:
            raise ValueError(f"Reference variable '{ref_var}' has unsupported shape in {ds.filepath()}")
        if depth_axis >= 0:
            return ref_var, None, depth_axis
        return ref_var, None, -1

    for name in ds.variables:
        ok, depth_axis = analyze(name)
        if not ok:
            continue
        return name, None, depth_axis

    raise ValueError(f"No supported reference variable found in {ds.filepath()}")


def _extract_time_seconds(ds: Dataset, time_name: str) -> np.ndarray:
    var = ds.variables[time_name]
    raw = np.array(var[:], dtype=np.float64)
    if raw.ndim != 1:
        raise ValueError(f"Reference time variable must be 1D: {time_name}")

    units = getattr(var, "units", None)
    if units and "since" in str(units):
        try:
            dts = num2date(raw, units=units, calendar=getattr(var, "calendar", "standard"))
            out = np.empty(len(dts), dtype=np.int64)
            for i, t in enumerate(dts):
                out[i] = _epoch_seconds(t)
            return out
        except Exception:
            pass

    return np.rint(raw).astype(np.int64)


@dataclass
class ReferenceSource:
    path: Path
    ds: Dataset
    time_name: str
    lat_name: str
    lon_name: str
    var_name: str
    depth_axis: int
    depth_index: int | None
    time_seconds: np.ndarray
    time_sorted_indices: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    cache_size: int
    _cache: OrderedDict[int, np.ndarray]

    @classmethod
    def open(
        cls,
        path: Path,
        time_name: str | None,
        lat_name: str | None,
        lon_name: str | None,
        var_name: str | None,
        depth_index: int | None,
        cache_size: int,
    ) -> "ReferenceSource":
        ds = Dataset(path, mode="r")
        time_name, lat_name, lon_name = _infer_time_lat_lon(ds, time_name, lat_name, lon_name)
        var_name, depth_name, depth_axis = _infer_reference_variable(
            ds,
            var_name,
            time_name,
            lat_name,
            lon_name,
        )
        del depth_name

        var = ds.variables[var_name]
        if depth_axis >= 0:
            if depth_index is None:
                depth_index = 0
            if depth_index < 0 or depth_index >= var.shape[depth_axis]:
                ds.close()
                raise ValueError(
                    f"reference_depth out of range for {path}: {depth_index} (size={var.shape[depth_axis]})"
                )
        elif depth_index is not None and depth_index != 0:
            ds.close()
            raise ValueError(f"reference_depth must be omitted or 0 for 3D reference variable in {path}")

        time_seconds = _extract_time_seconds(ds, time_name)
        if time_seconds.ndim != 1:
            ds.close()
            raise ValueError(f"time variable must be 1D in {path}")

        order = np.argsort(time_seconds)
        lat = np.array(ds.variables[lat_name][:], dtype=np.float32)
        lon = np.array(ds.variables[lon_name][:], dtype=np.float32)
        if lat.size == 0 or lon.size == 0:
            ds.close()
            raise ValueError(f"Invalid lat/lon in reference file: {path}")

        return cls(
            path=path,
            ds=ds,
            time_name=time_name,
            lat_name=lat_name,
            lon_name=lon_name,
            var_name=var_name,
            depth_axis=depth_axis,
            depth_index=depth_index,
            time_seconds=time_seconds[order],
            time_sorted_indices=order.astype(np.int64),
            lat=lat,
            lon=lon,
            cache_size=max(1, int(cache_size)),
            _cache=OrderedDict(),
        )

    def close(self) -> None:
        self.ds.close()

    def get_slice(self, local_index: int) -> np.ndarray:
        cached = self._cache.get(local_index)
        if cached is not None:
            self._cache.move_to_end(local_index)
            return cached

        var = self.ds.variables[self.var_name]
        indexes = [slice(None)] * var.ndim
        time_axis = var.dimensions.index(self.time_name)
        indexes[time_axis] = int(local_index)

        if self.depth_axis >= 0:
            assert self.depth_index is not None
            indexes[self.depth_axis] = int(self.depth_index)

        raw = np.array(var[tuple(indexes)], dtype=np.float32)
        fill_value = getattr(var, "_FillValue", None)
        if fill_value is not None:
            raw = np.where(raw == np.float32(fill_value), np.float32(np.nan), raw)
        missing_value = getattr(var, "missing_value", None)
        if missing_value is not None:
            raw = np.where(raw == np.float32(missing_value), np.float32(np.nan), raw)
        remaining_axes = [i for i, indexer in enumerate(indexes) if isinstance(indexer, slice)]
        lat_axis = list(var.dimensions).index(self.lat_name)
        lon_axis = list(var.dimensions).index(self.lon_name)
        lat_pos = remaining_axes.index(lat_axis)
        lon_pos = remaining_axes.index(lon_axis)

        if (lat_pos, lon_pos) != (0, 1):
            raw = np.moveaxis(raw, (lat_pos, lon_pos), (0, 1))

        raw = np.squeeze(raw)
        if raw.ndim != 2:
            raise ValueError(f"Invalid reference slice shape at {self.path}")

        self._cache[local_index] = raw
        while len(self._cache) > self.cache_size:
            self._cache.popitem(last=False)
        return raw

    def nearest_grid_indices(self, latitudes: np.ndarray, longitudes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        lon_norm = np.array(normalize_reference_lon(longitudes, self.lon), dtype=np.float64)
        return nearest_index_1d(self.lat, latitudes), nearest_index_1d(self.lon, lon_norm)


class ReferenceReader:
    def __init__(
        self,
        paths: Iterable[str | Path],
        time_name: str | None = None,
        lat_name: str | None = None,
        lon_name: str | None = None,
        ref_var_name: str | None = None,
        reference_depth: int = 0,
        cache_size: int = 32,
    ) -> None:
        self.sources: list[ReferenceSource] = []
        self.time_seconds: np.ndarray = np.empty(0, dtype=np.int64)
        self.source_index: np.ndarray = np.empty(0, dtype=np.int64)
        self.time_local_index: np.ndarray = np.empty(0, dtype=np.int64)
        self.reference_depth = reference_depth
        self._collect(paths, time_name, lat_name, lon_name, ref_var_name, reference_depth, cache_size)

    def _collect(
        self,
        paths: Iterable[str | Path],
        time_name: str | None,
        lat_name: str | None,
        lon_name: str | None,
        ref_var_name: str | None,
        reference_depth: int,
        cache_size: int,
    ) -> None:
        entries: list[tuple[int, int, int]] = []
        for source_id, path in enumerate(map(Path, paths)):
            source = ReferenceSource.open(
                Path(path),
                time_name=time_name,
                lat_name=lat_name,
                lon_name=lon_name,
                var_name=ref_var_name,
                depth_index=reference_depth,
                cache_size=cache_size,
            )
            self.sources.append(source)

            for local_sorted_idx, t in enumerate(source.time_seconds):
                original_index = int(source.time_sorted_indices[local_sorted_idx])
                entries.append((int(t), source_id, original_index))

        if not entries:
            raise ValueError("No valid time entries found in reference files.")

        entries.sort(key=lambda x: x[0])
        self.time_seconds = np.array([entry[0] for entry in entries], dtype=np.int64)
        self.source_index = np.array([entry[1] for entry in entries], dtype=np.int64)
        self.time_local_index = np.array([entry[2] for entry in entries], dtype=np.int64)

    def close(self) -> None:
        for source in self.sources:
            source.close()

    def sample(
        self,
        obs_seconds: np.ndarray,
        latitudes: np.ndarray,
        longitudes: np.ndarray,
        require_within_reference: bool = False,
    ) -> np.ndarray:
        if obs_seconds.size == 0:
            return np.empty(0, dtype=np.float32)

        if self.time_seconds.size == 0:
            return np.full(obs_seconds.shape, np.float32(np.nan), dtype=np.float32)

        obs_seconds_arr = np.asarray(obs_seconds, dtype=np.int64)
        if require_within_reference:
            min_t = int(self.time_seconds[0]) // SEC_PER_DAY
            max_t = int(self.time_seconds[-1]) // SEC_PER_DAY
            obs_days = obs_seconds_arr // SEC_PER_DAY
            in_range = (obs_days >= min_t) & (obs_days <= max_t)
        else:
            in_range = np.ones(obs_seconds_arr.shape, dtype=bool)

        idx = np.searchsorted(self.time_seconds, obs_seconds_arr, side="left")
        left = np.clip(idx - 1, 0, self.time_seconds.size - 1)
        right = np.clip(idx, 0, self.time_seconds.size - 1)
        use_right = (obs_seconds_arr - self.time_seconds[left]) > (self.time_seconds[right] - obs_seconds_arr)
        nearest = np.where(use_right, right, left).astype(np.int64)

        out = np.full(obs_seconds.shape, np.float32(np.nan), dtype=np.float32)
        if not np.any(in_range):
            return out

        source_ids = self.source_index[nearest]
        local_ids = self.time_local_index[nearest]

        for source_id in np.unique(source_ids):
            pos = np.nonzero(source_ids == source_id)[0]
            pos = pos[in_range[pos]]
            if pos.size == 0:
                continue

            src = self.sources[int(source_id)]
            local_pos = local_ids[pos]
            for local_idx in np.unique(local_pos):
                target = pos[local_pos == local_idx]
                lat_idx, lon_idx = src.nearest_grid_indices(
                    latitudes[target].astype(np.float32),
                    longitudes[target].astype(np.float32),
                )
                grid = src.get_slice(int(local_idx))
                out[target] = grid[lat_idx, lon_idx]

        return out


def _normalize_year_value(value) -> int:
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.isdigit():
            return int(stripped)
    raise ValueError(f"Invalid year value in reference JSON: {value}")


def _append_reference_paths(paths: list[Path], value) -> None:
    if value is None:
        return
    if isinstance(value, (str, Path)):
        paths.append(Path(value))
        return
    if isinstance(value, list):
        for item in value:
            _append_reference_paths(paths, item)
        return
    if isinstance(value, dict):
        if "files" in value:
            _append_reference_paths(paths, value["files"])
            return
        if "paths" in value:
            _append_reference_paths(paths, value["paths"])
            return

        path_template = value.get("path_template")
        if path_template is not None:
            years = value.get("years")
            if years is None:
                start_year = value.get("start_year")
                end_year = value.get("end_year")
                if start_year is None or end_year is None:
                    raise ValueError("path_template requires years or start_year/end_year.")
                start_i = _normalize_year_value(start_year)
                end_i = _normalize_year_value(end_year)
                if end_i < start_i:
                    raise ValueError("start_year must be <= end_year.")
                years = list(range(start_i, end_i + 1))

            if not isinstance(years, list):
                raise ValueError("years must be a list.")
            for y in years:
                y_i = _normalize_year_value(y)
                template = str(path_template)
                paths.append(Path(template.format(year=y_i)))
            return

        raise ValueError("Unsupported reference-dict object format.")

        raise ValueError(f"Unsupported reference entry type in JSON: {type(value)}")


def _load_reference_json(raw: str):
    candidate = Path(raw)
    if candidate.exists():
        with candidate.open("r", encoding="utf-8") as f:
            return json.load(f)

    stripped = raw.strip()
    if stripped.startswith(("{", "[")) and stripped.endswith(("}", "]")):
        try:
            return json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON in --reference-dict: {raw}") from exc

    raise FileNotFoundError(f"Reference dict file not found: {candidate}")


def load_reference_paths(reference: list[str] | None, reference_dict: str | None) -> list[Path]:
    paths: list[Path] = []
    if reference_dict:
        data = _load_reference_json(reference_dict)
        if isinstance(data, list):
            for item in data:
                _append_reference_paths(paths, item)
        elif isinstance(data, dict):
            if "path_template" in data or "files" in data or "paths" in data:
                _append_reference_paths(paths, data)
            else:
                for _, value in sorted(data.items(), key=lambda kv: str(kv[0])):
                    _append_reference_paths(paths, value)
        else:
            raise ValueError("--reference-dict must be JSON array or object.")
    if reference:
        paths.extend(Path(p) for p in reference)
    if not paths:
        raise ValueError("No reference files specified.")
    return paths
