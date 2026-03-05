#!/usr/bin/env python3
from __future__ import annotations

import datetime as dt
import re
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

import numpy as np
from netCDF4 import Dataset

SLOT_PATTERN = re.compile(r"^soil_moisture_([AD])_(\d{2})$")
SEC_PER_DAY = 86400
EPOCH = dt.datetime(1970, 1, 1)


def _to_datetime(value: dt.datetime | str) -> dt.datetime:
    if isinstance(value, dt.datetime):
        return value
    # Supports "YYYY-MM-DDTHH:MM:SS" and "YYYY-MM-DD HH:MM:SS"
    return dt.datetime.fromisoformat(value)


def _to_epoch_seconds(value: dt.datetime) -> float:
    return (value - EPOCH).total_seconds()


def _from_epoch_seconds_array(values: np.ndarray) -> np.ndarray:
    return values.astype("datetime64[s]")


@dataclass(frozen=True)
class DayChunkRef:
    file_path: Path
    time_index: int


@dataclass
class DayObservationArray:
    sm: np.ndarray
    obs_sec: np.ndarray
    lat: np.ndarray
    lon: np.ndarray


class AmsrL3ObservationReader:
    """
    Fast reader for repeated datetime-range queries.

    Design for iterative queries:
      1) Build day index once at initialization.
      2) Read one day to flattened arrays on demand.
      3) Keep day arrays in LRU cache and reuse across queries.
    """

    def __init__(
        self,
        path: str | Path,
        max_cache_days: int | None = None,
        interval_hours: float | None = None,
        window_hours: float | None = None,
    ):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"Path not found: {self.path}")

        if max_cache_days is None:
            if interval_hours is None:
                max_cache_days = 8
            else:
                max_cache_days = self.suggest_cache_days(
                    interval_hours=interval_hours,
                    window_hours=window_hours,
                )
        if max_cache_days < 1:
            raise ValueError("max_cache_days must be >= 1")

        self.max_cache_days = max_cache_days
        self.day_index: Dict[int, List[DayChunkRef]] = defaultdict(list)
        self.file_slot_pairs: Dict[Path, List[Tuple[str, str]]] = {}
        self.cache: OrderedDict[int, DayObservationArray] = OrderedDict()
        self._lat: np.ndarray | None = None
        self._lon: np.ndarray | None = None

        self._build_index()

    @staticmethod
    def suggest_cache_days(
        interval_hours: float,
        window_hours: float | None = None,
    ) -> int:
        """
        Recommend cache size (days) from fixed query interval/window.

        Heuristic:
          cache_days = ceil((step + window) / 24h) + 1
        where `step=interval_hours`, `window=window_hours or interval_hours`.
        """
        if interval_hours <= 0:
            raise ValueError("interval_hours must be > 0")
        win = interval_hours if window_hours is None else window_hours
        if win <= 0:
            raise ValueError("window_hours must be > 0")

        days = int(np.ceil((interval_hours + win) / 24.0) + 1)
        return max(1, days)

    def _list_nc_files(self) -> List[Path]:
        if self.path.is_file():
            return [self.path]
        files = sorted(self.path.glob("*.nc"))
        if not files:
            raise FileNotFoundError(f"No NetCDF files found in directory: {self.path}")
        return files

    @staticmethod
    def _collect_slot_pairs(ds: Dataset) -> List[Tuple[str, str]]:
        sm_names = [name for name in ds.variables if SLOT_PATTERN.match(name)]
        pairs: List[Tuple[str, str]] = []
        for sm_name in sorted(sm_names):
            suffix = sm_name.replace("soil_moisture_", "", 1)
            tm_name = f"observation_time_min_{suffix}"
            if tm_name in ds.variables:
                pairs.append((sm_name, tm_name))
        return pairs

    def _build_index(self) -> None:
        for fp in self._list_nc_files():
            with Dataset(fp) as ds:
                if "time" not in ds.variables or "lat" not in ds.variables or "lon" not in ds.variables:
                    continue

                slot_pairs = self._collect_slot_pairs(ds)
                if not slot_pairs:
                    continue
                self.file_slot_pairs[fp] = slot_pairs

                if self._lat is None:
                    self._lat = np.array(ds.variables["lat"][:], dtype=np.float32)
                    self._lon = np.array(ds.variables["lon"][:], dtype=np.float32)

                day_values = np.array(ds.variables["time"][:], dtype=np.int64)
                for t_idx, day in enumerate(day_values):
                    self.day_index[int(day)].append(DayChunkRef(file_path=fp, time_index=int(t_idx)))

        if self._lat is None or self._lon is None or not self.day_index:
            raise ValueError(f"No AMSR merged observation variables found under: {self.path}")

    def _put_cache(self, day_key: int, day_data: DayObservationArray) -> None:
        self.cache[day_key] = day_data
        self.cache.move_to_end(day_key)
        while len(self.cache) > self.max_cache_days:
            self.cache.popitem(last=False)

    def clear_cache(self) -> None:
        self.cache.clear()

    def _load_day(self, day_key: int) -> DayObservationArray:
        cached = self.cache.get(day_key)
        if cached is not None:
            self.cache.move_to_end(day_key)
            return cached

        if day_key not in self.day_index:
            empty = DayObservationArray(
                sm=np.empty(0, dtype=np.float32),
                obs_sec=np.empty(0, dtype=np.int64),
                lat=np.empty(0, dtype=np.float32),
                lon=np.empty(0, dtype=np.float32),
            )
            self._put_cache(day_key, empty)
            return empty

        sm_parts: List[np.ndarray] = []
        sec_parts: List[np.ndarray] = []
        lat_parts: List[np.ndarray] = []
        lon_parts: List[np.ndarray] = []

        lat_axis = self._lat
        lon_axis = self._lon
        assert lat_axis is not None and lon_axis is not None
        day_sec0 = int(day_key) * SEC_PER_DAY

        for chunk in self.day_index[day_key]:
            slot_pairs = self.file_slot_pairs.get(chunk.file_path, [])
            if not slot_pairs:
                continue

            with Dataset(chunk.file_path) as ds:
                for sm_name, tm_name in slot_pairs:
                    sm_2d = np.array(ds.variables[sm_name][chunk.time_index, :, :], dtype=np.float32)
                    tm_2d = np.array(ds.variables[tm_name][chunk.time_index, :, :], dtype=np.float32)

                    valid = np.isfinite(sm_2d) & np.isfinite(tm_2d)
                    if not np.any(valid):
                        continue

                    lat_idx, lon_idx = np.nonzero(valid)
                    sm_vals = sm_2d[valid].astype(np.float32, copy=False)
                    tm_vals = tm_2d[valid].astype(np.float64, copy=False)
                    obs_sec = (day_sec0 + np.rint(tm_vals * 60.0)).astype(np.int64, copy=False)

                    sm_parts.append(sm_vals)
                    sec_parts.append(obs_sec)
                    lat_parts.append(lat_axis[lat_idx])
                    lon_parts.append(lon_axis[lon_idx])

        if not sm_parts:
            day_data = DayObservationArray(
                sm=np.empty(0, dtype=np.float32),
                obs_sec=np.empty(0, dtype=np.int64),
                lat=np.empty(0, dtype=np.float32),
                lon=np.empty(0, dtype=np.float32),
            )
            self._put_cache(day_key, day_data)
            return day_data

        sm_all = np.concatenate(sm_parts)
        sec_all = np.concatenate(sec_parts)
        lat_all = np.concatenate(lat_parts)
        lon_all = np.concatenate(lon_parts)

        order = np.argsort(sec_all, kind="mergesort")
        day_data = DayObservationArray(
            sm=sm_all[order],
            obs_sec=sec_all[order],
            lat=lat_all[order],
            lon=lon_all[order],
        )
        self._put_cache(day_key, day_data)
        return day_data

    def _read_range_arrays(
        self,
        start_datetime: dt.datetime | str,
        end_datetime: dt.datetime | str,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        start = _to_datetime(start_datetime)
        end = _to_datetime(end_datetime)
        if end < start:
            raise ValueError("end_datetime must be >= start_datetime")

        start_sec = _to_epoch_seconds(start)
        end_sec = _to_epoch_seconds(end)
        start_day = int(np.floor(start_sec / SEC_PER_DAY))
        end_day = int(np.floor(end_sec / SEC_PER_DAY))

        sm_parts: List[np.ndarray] = []
        sec_parts: List[np.ndarray] = []
        lat_parts: List[np.ndarray] = []
        lon_parts: List[np.ndarray] = []

        for day in range(start_day, end_day + 1):
            day_data = self._load_day(day)
            if day_data.sm.size == 0:
                continue
            mask = (day_data.obs_sec >= start_sec) & (day_data.obs_sec <= end_sec)
            if not np.any(mask):
                continue
            sm_parts.append(day_data.sm[mask])
            sec_parts.append(day_data.obs_sec[mask])
            lat_parts.append(day_data.lat[mask])
            lon_parts.append(day_data.lon[mask])

        if not sm_parts:
            return (
                np.empty(0, dtype=np.float32),
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float32),
                np.empty(0, dtype=np.float32),
            )

        sm = np.concatenate(sm_parts)
        obs_sec = np.concatenate(sec_parts)
        lat = np.concatenate(lat_parts)
        lon = np.concatenate(lon_parts)

        order = np.argsort(obs_sec, kind="mergesort")
        return (
            sm[order].astype(np.float32, copy=False),
            obs_sec[order],
            lat[order].astype(np.float32, copy=False),
            lon[order].astype(np.float32, copy=False),
        )

    def preload_range(self, start_datetime: dt.datetime | str, end_datetime: dt.datetime | str) -> None:
        """
        Preload day cache for iterative queries over the range.
        """
        start = _to_datetime(start_datetime)
        end = _to_datetime(end_datetime)
        if end < start:
            raise ValueError("end_datetime must be >= start_datetime")

        start_day = int(np.floor(_to_epoch_seconds(start) / SEC_PER_DAY))
        end_day = int(np.floor(_to_epoch_seconds(end) / SEC_PER_DAY))
        for day in range(start_day, end_day + 1):
            self._load_day(day)

    def read_range(
        self,
        start_datetime: dt.datetime | str,
        end_datetime: dt.datetime | str,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Return all observations in [start_datetime, end_datetime] as arrays:
          (soil_moisture, observation_time, lat, lon)
        where observation_time is datetime64[s].
        """
        start = _to_datetime(start_datetime)
        end = _to_datetime(end_datetime)
        if end < start:
            raise ValueError("end_datetime must be >= start_datetime")
        sm, obs_sec, lat, lon = self._read_range_arrays(start, end)
        return sm, _from_epoch_seconds_array(obs_sec), lat, lon

    def iterate_windows(
        self,
        start_datetime: dt.datetime | str,
        end_datetime: dt.datetime | str,
        step: dt.timedelta,
        window: dt.timedelta | None = None,
    ) -> Iterator[Tuple[dt.datetime, dt.datetime, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        """
        Iterate range queries efficiently with shared cache.

        Examples:
          - 3-hour windows: step=timedelta(hours=3), window=None
          - 6-hour windows: step=timedelta(hours=6), window=None
          - daily windows:  step=timedelta(days=1), window=None
          - 3-hour step with 1-day window: step=timedelta(hours=3), window=timedelta(days=1)
        """
        start = _to_datetime(start_datetime)
        end = _to_datetime(end_datetime)
        if end < start:
            raise ValueError("end_datetime must be >= start_datetime")
        if step.total_seconds() <= 0:
            raise ValueError("step must be positive")

        win = step if window is None else window
        if win.total_seconds() <= 0:
            raise ValueError("window must be positive")

        cur = start
        eps = dt.timedelta(seconds=1)
        while cur <= end:
            cur_end = min(cur + win, end)
            # To avoid boundary duplication across consecutive windows,
            # intermediate windows are treated as [cur, cur_end).
            query_end = cur_end if cur_end >= end else (cur_end - eps)
            if query_end >= cur:
                sm, obs_time, lat, lon = self.read_range(cur, query_end)
            else:
                sm = np.empty(0, dtype=np.float32)
                obs_time = np.empty(0, dtype="datetime64[s]")
                lat = np.empty(0, dtype=np.float32)
                lon = np.empty(0, dtype=np.float32)
            yield cur, cur_end, sm, obs_time, lat, lon
            cur = cur + step


def read_amsr_observations_in_range(
    path: str | Path,
    start_datetime: dt.datetime | str,
    end_datetime: dt.datetime | str,
    max_cache_days: int | None = None,
    interval_hours: float | None = None,
    window_hours: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Convenience wrapper.
    """
    reader = AmsrL3ObservationReader(
        path,
        max_cache_days=max_cache_days,
        interval_hours=interval_hours,
        window_hours=window_hours,
    )
    return reader.read_range(start_datetime, end_datetime)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Read AMSR observations in datetime range.")
    parser.add_argument("path", type=Path, help="Merged NetCDF file or directory containing yearly NetCDFs.")
    parser.add_argument("start", type=str, help="Start datetime (ISO format).")
    parser.add_argument("end", type=str, help="End datetime (ISO format).")
    parser.add_argument("--limit", type=int, default=10, help="Print first N records.")
    parser.add_argument("--cache-days", type=int, default=None, help="LRU cache size in days.")
    parser.add_argument(
        "--interval-hours",
        type=float,
        default=None,
        help="Fixed query step (hours). Used to auto-decide cache-days when --cache-days is omitted.",
    )
    parser.add_argument(
        "--window-hours",
        type=float,
        default=None,
        help="Query window length (hours). Defaults to interval-hours when omitted.",
    )
    args = parser.parse_args()

    sm, obs_time, lat, lon = read_amsr_observations_in_range(
        args.path,
        args.start,
        args.end,
        max_cache_days=args.cache_days,
        interval_hours=args.interval_hours,
        window_hours=args.window_hours,
    )
    print(f"count={len(sm)}")
    for i in range(min(args.limit, len(sm))):
        print(sm[i], obs_time[i], lat[i], lon[i])
