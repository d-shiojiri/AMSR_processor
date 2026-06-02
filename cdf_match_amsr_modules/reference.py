from __future__ import annotations

import datetime as dt
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import netCDF4
import numpy as np
from netCDF4 import Dataset

from .io_utils import SEC_PER_DAY, nearest_index_1d, normalize_lon


EPOCH = dt.datetime(1970, 1, 1)


def _epoch_seconds(value) -> int:
    if isinstance(value, np.datetime64):
        return int(value.astype("datetime64[s]").astype(np.int64))
    if isinstance(value, dt.datetime):
        return int((value - EPOCH).total_seconds())
    if hasattr(value, "year") and hasattr(value, "month") and hasattr(value, "day"):
        return int(
            (
                dt.datetime(
                    int(value.year),
                    int(value.month),
                    int(value.day),
                    int(getattr(value, "hour", 0)),
                    int(getattr(value, "minute", 0)),
                    int(getattr(value, "second", 0)),
                )
                - EPOCH
            ).total_seconds()
        )
    raise TypeError(f"Unsupported datetime type: {type(value)}")


def _infer_time_lat_lon(ds: Dataset, time_name: str | None, lat_name: str | None, lon_name: str | None):
    tn = time_name or ("time" if "time" in ds.variables else None)
    ln = lat_name or ("lat" if "lat" in ds.variables else ("latitude" if "latitude" in ds.variables else None))
    lonn = lon_name or ("lon" if "lon" in ds.variables else ("longitude" if "longitude" in ds.variables else None))
    if tn is None or ln is None or lonn is None:
        raise ValueError(f"Cannot infer coordinate names in {ds.filepath()}")
    return tn, ln, lonn


def _infer_reference_var(ds: Dataset, var_name: str | None, time_name: str, lat_name: str, lon_name: str):
    def ok(name: str):
        dims = tuple(ds.variables[name].dimensions)
        required = {time_name, lat_name, lon_name}
        if len(dims) == 3 and dims == (time_name, lat_name, lon_name):
            return True, -1
        if len(dims) == 4 and required.issubset(dims):
            depth_axes = [i for i, d in enumerate(dims) if d not in required]
            if len(depth_axes) == 1:
                return True, depth_axes[0]
        return False, -1

    if var_name is not None:
        if var_name not in ds.variables:
            raise ValueError(f"Reference variable not found: {var_name}")
        is_ok, depth_axis = ok(var_name)
        if not is_ok:
            raise ValueError(f"Unsupported reference variable shape: {var_name}")
        return var_name, depth_axis
    for name in ds.variables:
        is_ok, depth_axis = ok(name)
        if is_ok:
            return name, depth_axis
    raise ValueError(f"No supported reference variable found in {ds.filepath()}")


def _time_seconds(ds: Dataset, time_name: str) -> np.ndarray:
    var = ds.variables[time_name]
    raw = np.asarray(var[:], dtype=np.float64)
    units = getattr(var, "units", None)
    if units and "since" in str(units):
        dts = netCDF4.num2date(raw, units, calendar=getattr(var, "calendar", "standard"))
        return np.array([_epoch_seconds(t) for t in dts], dtype=np.int64)
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
    depth_index: int
    time_seconds: np.ndarray
    lat: np.ndarray
    lon: np.ndarray

    @classmethod
    def open(
        cls,
        path: Path,
        time_name: str | None,
        lat_name: str | None,
        lon_name: str | None,
        var_name: str | None,
        depth_index: int,
    ) -> "ReferenceSource":
        ds = Dataset(path, "r")
        time_name, lat_name, lon_name = _infer_time_lat_lon(ds, time_name, lat_name, lon_name)
        var_name, depth_axis = _infer_reference_var(ds, var_name, time_name, lat_name, lon_name)
        var = ds.variables[var_name]
        if depth_axis >= 0 and not (0 <= depth_index < var.shape[depth_axis]):
            ds.close()
            raise ValueError(f"reference_depth out of range for {path}: {depth_index}")
        return cls(
            path=path,
            ds=ds,
            time_name=time_name,
            lat_name=lat_name,
            lon_name=lon_name,
            var_name=var_name,
            depth_axis=depth_axis,
            depth_index=depth_index,
            time_seconds=_time_seconds(ds, time_name),
            lat=np.asarray(ds.variables[lat_name][:], dtype=np.float64),
            lon=np.asarray(ds.variables[lon_name][:], dtype=np.float64),
        )

    def close(self) -> None:
        self.ds.close()

    def nearest_time_indices(self, obs_seconds: np.ndarray) -> np.ndarray:
        obs_seconds = np.asarray(obs_seconds, dtype=np.int64)
        pos = np.searchsorted(self.time_seconds, obs_seconds, side="left")
        left = np.clip(pos - 1, 0, self.time_seconds.size - 1)
        right = np.clip(pos, 0, self.time_seconds.size - 1)
        choose_left = (obs_seconds - self.time_seconds[left]) <= (self.time_seconds[right] - obs_seconds)
        return np.where(choose_left, left, right).astype(np.int64)

    def read_block(self, time_index: int, input_lat: np.ndarray, input_lon: np.ndarray, row_slice: slice) -> np.ndarray:
        lat_idx = nearest_index_1d(self.lat, input_lat[row_slice])
        lon_idx = nearest_index_1d(self.lon, normalize_lon(input_lon, self.lon))

        var = self.ds.variables[self.var_name]
        indexes = [slice(None)] * var.ndim
        indexes[var.dimensions.index(self.time_name)] = int(time_index)
        if self.depth_axis >= 0:
            indexes[self.depth_axis] = int(self.depth_index)
        raw = np.asarray(var[tuple(indexes)], dtype=np.float32)

        lat_axis = list(var.dimensions).index(self.lat_name)
        lon_axis = list(var.dimensions).index(self.lon_name)
        remaining_axes = [i for i, indexer in enumerate(indexes) if isinstance(indexer, slice)]
        lat_pos = remaining_axes.index(lat_axis)
        lon_pos = remaining_axes.index(lon_axis)
        if (lat_pos, lon_pos) != (0, 1):
            raw = np.moveaxis(raw, (lat_pos, lon_pos), (0, 1))
        return raw[np.ix_(lat_idx, lon_idx)].astype(np.float32, copy=False)


class ReferenceReader:
    def __init__(
        self,
        paths: Iterable[str | Path],
        time_name: str | None = None,
        lat_name: str | None = None,
        lon_name: str | None = None,
        ref_var_name: str | None = None,
        reference_depth: int = 0,
    ) -> None:
        self.sources = [
            ReferenceSource.open(Path(p), time_name, lat_name, lon_name, ref_var_name, reference_depth)
            for p in paths
        ]
        if not self.sources:
            raise ValueError("No reference files specified")
        self.time_seconds = np.concatenate([s.time_seconds for s in self.sources])
        self.source_ids = np.concatenate([
            np.full(s.time_seconds.size, i, dtype=np.int64) for i, s in enumerate(self.sources)
        ])
        self.local_ids = np.concatenate([np.arange(s.time_seconds.size, dtype=np.int64) for s in self.sources])
        order = np.argsort(self.time_seconds)
        self.time_seconds = self.time_seconds[order]
        self.source_ids = self.source_ids[order]
        self.local_ids = self.local_ids[order]

    def close(self) -> None:
        for src in self.sources:
            src.close()

    def nearest_lookup(self, obs_seconds: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        obs_seconds = np.asarray(obs_seconds, dtype=np.int64)
        pos = np.searchsorted(self.time_seconds, obs_seconds, side="left")
        left = np.clip(pos - 1, 0, self.time_seconds.size - 1)
        right = np.clip(pos, 0, self.time_seconds.size - 1)
        choose_left = (obs_seconds - self.time_seconds[left]) <= (self.time_seconds[right] - obs_seconds)
        nearest = np.where(choose_left, left, right).astype(np.int64)
        return self.source_ids[nearest], self.local_ids[nearest]


def _append_reference_paths(paths: list[Path], value) -> None:
    if value is None:
        return
    if isinstance(value, (str, Path)):
        paths.append(Path(value))
    elif isinstance(value, list):
        for item in value:
            _append_reference_paths(paths, item)
    elif isinstance(value, dict):
        if "files" in value:
            _append_reference_paths(paths, value["files"])
        elif "paths" in value:
            _append_reference_paths(paths, value["paths"])
        elif "path_template" in value:
            template = str(value["path_template"])
            years = value.get("years")
            if years is None:
                years = range(int(value["start_year"]), int(value["end_year"]) + 1)
            for year in years:
                paths.append(Path(template.format(year=int(year))))
        else:
            for item in value.values():
                _append_reference_paths(paths, item)
    else:
        raise ValueError(f"Unsupported reference entry type: {type(value)}")


def load_reference_paths(reference: list[str] | None, reference_dict: str | None) -> list[Path]:
    paths: list[Path] = []
    if reference_dict:
        p = Path(reference_dict)
        if p.exists():
            data = json.loads(p.read_text())
        else:
            data = json.loads(reference_dict)
        _append_reference_paths(paths, data)
    if reference:
        paths.extend(Path(p) for p in reference)
    if not paths:
        raise ValueError("No reference files specified.")
    return paths
