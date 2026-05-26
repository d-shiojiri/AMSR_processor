#!/usr/bin/env python3
from __future__ import annotations

import datetime as dt
from collections import OrderedDict
from pathlib import Path
from typing import Iterator

import numpy as np
from netCDF4 import Dataset

from .amsr_l3_query import AmsrL3ObservationReader


EPOCH = dt.datetime(1970, 1, 1)
SEC_PER_DAY = 86400
ObservationArrayResult = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
EMPTY_RESULT: ObservationArrayResult = (
    np.empty(0, dtype=np.float32),
    np.empty(0, dtype="datetime64[s]"),
    np.empty(0, dtype=np.float32),
    np.empty(0, dtype=np.float32),
)


def _to_datetime(value: dt.datetime | str) -> dt.datetime:
    return value if isinstance(value, dt.datetime) else dt.datetime.fromisoformat(value)


def _to_epoch_seconds(value: dt.datetime) -> int:
    return int((value - EPOCH).total_seconds())


def _list_input_files(path: Path) -> list[Path]:
    return [path] if path.is_file() else sorted(path.glob("*.nc"))


def _is_averaged_file(path: Path) -> bool:
    with Dataset(path) as ds:
        return "soil_moisture" in ds.variables and "observation_count" in ds.variables


def _empty_day() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.empty(0, dtype=np.float32),
        np.empty(0, dtype=np.int64),
        np.empty(0, dtype=np.float32),
        np.empty(0, dtype=np.float32),
    )


class Amsr0p5AveragedReader:
    """Reader for 0.5-degree daily AMSR files."""

    def __init__(self, path: str | Path, max_cache_days: int | None = None) -> None:
        self.path = Path(path)
        self.input_files = _list_input_files(self.path)
        self.max_cache_days = max_cache_days or 8
        self.cache: OrderedDict[int, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = OrderedDict()
        self.day_index: dict[int, list[tuple[Path, int]]] = {}
        self._lat: np.ndarray | None = None
        self._lon: np.ndarray | None = None
        self._build_index()

    def _build_index(self) -> None:
        # Build day -> [(file, time_index)] and keep the shared 0.5-degree axes.
        for fp in self.input_files:
            with Dataset(fp) as ds:
                if self._lat is None:
                    self._lat = np.array(ds.variables["lat"][:], dtype=np.float32)
                    self._lon = np.array(ds.variables["lon"][:], dtype=np.float32)

                for t_idx, day in enumerate(np.array(ds.variables["time"][:], dtype=np.int64)):
                    self.day_index.setdefault(int(day), []).append((fp, int(t_idx)))

    def _put_cache(self, day: int, data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]) -> None:
        # Small LRU cache avoids reopening the same day in sliding-window queries.
        self.cache[day] = data
        self.cache.move_to_end(day)
        while len(self.cache) > self.max_cache_days:
            self.cache.popitem(last=False)

    def _load_day(self, day: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        cached = self.cache.get(day)
        if cached is not None:
            self.cache.move_to_end(day)
            return cached

        # Read all valid grid cells for one day and flatten them to observation arrays.
        sm_parts, sec_parts, lat_parts, lon_parts = [], [], [], []
        day_sec0 = day * SEC_PER_DAY
        for file_path, time_index in self.day_index.get(day, []):
            with Dataset(file_path) as ds:
                sm_2d = np.array(ds.variables["soil_moisture"][time_index], dtype=np.float32)
                count_2d = np.array(ds.variables["observation_count"][time_index], dtype=np.int32)
                valid = np.isfinite(sm_2d) & (count_2d > 0)
                if not np.any(valid):
                    continue

                lat_idx, lon_idx = np.nonzero(valid)
                obs_sec = np.full(lat_idx.size, day_sec0, dtype=np.int64)
                if "observation_time_min" in ds.variables:
                    minutes = np.array(ds.variables["observation_time_min"][time_index], dtype=np.float32)[valid]
                    finite = np.isfinite(minutes)
                    obs_sec[finite] += np.rint(minutes[finite] * 60.0).astype(np.int64, copy=False)

                sm_parts.append(sm_2d[valid])
                sec_parts.append(obs_sec)
                lat_parts.append(self._lat[lat_idx])
                lon_parts.append(self._lon[lon_idx])

        if sm_parts:
            data = (
                np.concatenate(sm_parts).astype(np.float32, copy=False),
                np.concatenate(sec_parts).astype(np.int64, copy=False),
                np.concatenate(lat_parts).astype(np.float32, copy=False),
                np.concatenate(lon_parts).astype(np.float32, copy=False),
            )
        else:
            data = _empty_day()

        self._put_cache(day, data)
        return data

    def read_range(self, start_datetime: dt.datetime | str, end_datetime: dt.datetime | str) -> ObservationArrayResult:
        """Return (soil_moisture[%], observation_time, lat, lon) in the requested time range."""
        start = _to_datetime(start_datetime)
        end = _to_datetime(end_datetime)
        start_sec = _to_epoch_seconds(start)
        end_sec = _to_epoch_seconds(end)
        start_day = start_sec // SEC_PER_DAY
        end_day = end_sec // SEC_PER_DAY

        # Load only the days touched by the query, then mask to the exact seconds.
        sm_parts, sec_parts, lat_parts, lon_parts = [], [], [], []
        for day in range(start_day, end_day + 1):
            sm, sec, lat, lon = self._load_day(day)
            mask = (sec >= start_sec) & (sec <= end_sec)
            if np.any(mask):
                sm_parts.append(sm[mask])
                sec_parts.append(sec[mask])
                lat_parts.append(lat[mask])
                lon_parts.append(lon[mask])

        if not sm_parts:
            return EMPTY_RESULT

        sm = np.concatenate(sm_parts)
        sec = np.concatenate(sec_parts)
        lat = np.concatenate(lat_parts)
        lon = np.concatenate(lon_parts)
        order = np.argsort(sec, kind="mergesort")
        return (
            sm[order].astype(np.float32, copy=False),
            sec[order].astype("datetime64[s]", copy=False),
            lat[order].astype(np.float32, copy=False),
            lon[order].astype(np.float32, copy=False),
        )

    def iterate_windows(
        self,
        start_datetime: dt.datetime | str,
        end_datetime: dt.datetime | str,
        step: dt.timedelta,
        window: dt.timedelta | None = None,
    ) -> Iterator[tuple[dt.datetime, dt.datetime, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        # Iterate [cur, cur+window) windows. The yielded cur_end keeps the original boundary.
        start = _to_datetime(start_datetime)
        end = _to_datetime(end_datetime)
        win = step if window is None else window
        cur = start
        eps = dt.timedelta(seconds=1)
        while cur <= end:
            cur_end = min(cur + win, end)
            query_end = cur_end if cur_end >= end else cur_end - eps
            yield (cur, cur_end, *self.read_range(cur, query_end))
            cur += step


class AmsrUpsampledL3ObservationReader:
    """Compatibility reader for averaged 0.5-degree files and legacy slot files."""

    def __init__(
        self,
        path: str | Path,
        max_cache_days: int | None = None,
        interval_hours: float | None = None,
        window_hours: float | None = None,
    ) -> None:
        path = Path(path)
        first = _list_input_files(path)[0]

        # New 0.5-degree averaged files are handled here; old slot files are delegated.
        if _is_averaged_file(first):
            cache_days = max_cache_days or (
                AmsrL3ObservationReader.suggest_cache_days(interval_hours, window_hours)
                if interval_hours is not None
                else 8
            )
            self._inner = Amsr0p5AveragedReader(path, max_cache_days=cache_days)
        else:
            self._inner = AmsrL3ObservationReader(
                path,
                max_cache_days=max_cache_days,
                interval_hours=interval_hours,
                window_hours=window_hours,
            )

    def read_range(self, start_datetime: dt.datetime | str, end_datetime: dt.datetime | str) -> ObservationArrayResult:
        return self._inner.read_range(start_datetime, end_datetime)

    def iterate_windows(
        self,
        start_datetime: dt.datetime | str,
        end_datetime: dt.datetime | str,
        step: dt.timedelta,
        window: dt.timedelta | None = None,
    ) -> Iterator[tuple[dt.datetime, dt.datetime, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        return self._inner.iterate_windows(start_datetime, end_datetime, step, window)


def read_upsampled_amsr_observations_in_range(
    path: str | Path = "processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc",
    start_datetime: dt.datetime | str | None = None,
    end_datetime: dt.datetime | str | None = None,
    max_cache_days: int | None = None,
    interval_hours: float | None = None,
    window_hours: float | None = None,
) -> ObservationArrayResult:
    reader = AmsrUpsampledL3ObservationReader(
        path,
        max_cache_days=max_cache_days,
        interval_hours=interval_hours,
        window_hours=window_hours,
    )
    return reader.read_range(start_datetime, end_datetime)


def _parse_args() -> "argparse.Namespace":
    import argparse

    parser = argparse.ArgumentParser(description="Read AMSR 0.5-degree observations in datetime range.")
    parser.add_argument("path", type=Path, nargs="?", default=Path("processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc"))
    parser.add_argument("start", type=str)
    parser.add_argument("end", type=str)
    parser.add_argument("--limit", type=int, default=10)
    parser.add_argument("--cache-days", type=int, default=None)
    parser.add_argument("--interval-hours", type=float, default=None)
    parser.add_argument("--window-hours", type=float, default=None)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    sm, obs_time, lat, lon = read_upsampled_amsr_observations_in_range(
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


if __name__ == "__main__":
    main()
