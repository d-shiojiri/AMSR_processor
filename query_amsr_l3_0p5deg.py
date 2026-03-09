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


def _to_datetime(value: dt.datetime | str) -> dt.datetime:
    if isinstance(value, dt.datetime):
        return value
    return dt.datetime.fromisoformat(value)


def _to_epoch_seconds(value: dt.datetime) -> float:
    return (value - EPOCH).total_seconds()


def _list_input_files(path: Path) -> list[Path]:
    if path.is_file():
        if path.suffix.lower() != ".nc":
            raise ValueError(f"Input file must be .nc: {path}")
        return [path]

    files = sorted(path.glob("*.nc"))
    if not files:
        raise FileNotFoundError(f"No .nc files found under: {path}")
    return files


def _is_averaged_file(path: Path) -> bool:
    with Dataset(path) as ds:
        return (
            "soil_moisture" in ds.variables
            and "observation_count" in ds.variables
            and "time" in ds.variables
        )


class Amsr0p5AveragedReader:
    """
    Reader for averaged 0.5-degree daily NetCDF output.
    """

    def __init__(
        self,
        path: str | Path,
        max_cache_days: int | None = None,
    ) -> None:
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"Path not found: {self.path}")

        self.input_files = _list_input_files(self.path)
        if max_cache_days is None:
            max_cache_days = 8
        if max_cache_days < 1:
            raise ValueError("max_cache_days must be >= 1")
        self.max_cache_days = max_cache_days

        self.cache: OrderedDict[int, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = OrderedDict()
        self.day_index: dict[int, list[tuple[Path, int]]] = {}
        self._lat: np.ndarray | None = None
        self._lon: np.ndarray | None = None
        self._build_index()

    def _build_index(self) -> None:
        for fp in self.input_files:
            with Dataset(fp) as ds:
                if "time" not in ds.variables or "lat" not in ds.variables or "lon" not in ds.variables:
                    continue
                if "soil_moisture" not in ds.variables:
                    raise ValueError(f"{fp} lacks soil_moisture variable")
                if "observation_count" not in ds.variables:
                    raise ValueError(f"{fp} lacks observation_count variable")

                if self._lat is None:
                    self._lat = np.array(ds.variables["lat"][:], dtype=np.float32)
                    self._lon = np.array(ds.variables["lon"][:], dtype=np.float32)

                times = np.array(ds.variables["time"][:], dtype=np.int64)
                if times.ndim != 1:
                    raise ValueError(f"{fp} time variable is not 1-D.")

                for t_idx, day in enumerate(times):
                    day_i = int(day)
                    self.day_index.setdefault(day_i, []).append((fp, int(t_idx)))

        if self._lat is None or self._lon is None or not self.day_index:
            raise ValueError(f"No averaged AMSR 0.5-degree files found under: {self.path}")

    def _put_cache(self, day_key: int, day_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]) -> None:
        self.cache[day_key] = day_data
        self.cache.move_to_end(day_key)
        while len(self.cache) > self.max_cache_days:
            self.cache.popitem(last=False)

    def _load_day(self, day_key: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        cached = self.cache.get(day_key)
        if cached is not None:
            self.cache.move_to_end(day_key)
            return cached

        refs = self.day_index.get(day_key, [])
        if not refs:
            empty = (
                np.empty(0, dtype=np.float32),
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float32),
                np.empty(0, dtype=np.float32),
            )
            self._put_cache(day_key, empty)
            return empty

        sm_parts: list[np.ndarray] = []
        sec_parts: list[np.ndarray] = []
        lat_parts: list[np.ndarray] = []
        lon_parts: list[np.ndarray] = []
        day_sec0 = day_key * SEC_PER_DAY

        for file_path, time_index in refs:
            with Dataset(file_path) as ds:
                sm_2d = np.array(ds.variables["soil_moisture"][time_index, :, :], dtype=np.float32)
                count_2d = np.array(ds.variables["observation_count"][time_index, :, :], dtype=np.int32)
                valid = np.isfinite(sm_2d) & (count_2d > 0)
                if not np.any(valid):
                    continue

                lat_idx, lon_idx = np.nonzero(valid)
                sm_parts.append(sm_2d[valid])
                obs_sec = np.full(len(lat_idx), int(day_sec0), dtype=np.int64)
                if "observation_time_min" in ds.variables:
                    tm_2d = np.array(ds.variables["observation_time_min"][time_index, :, :], dtype=np.float32)
                    tm_vals = tm_2d[valid]
                    finite_tm = np.isfinite(tm_vals)
                    if np.any(finite_tm):
                        obs_sec[finite_tm] = int(day_sec0) + np.rint(tm_vals[finite_tm] * 60.0).astype(
                            np.int64,
                            copy=False,
                        )
                sec_parts.append(obs_sec)
                lat_parts.append(self._lat[lat_idx] if self._lat is not None else np.empty(0, dtype=np.float32))
                lon_parts.append(self._lon[lon_idx] if self._lon is not None else np.empty(0, dtype=np.float32))

        if not sm_parts:
            day_data = (
                np.empty(0, dtype=np.float32),
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float32),
                np.empty(0, dtype=np.float32),
            )
            self._put_cache(day_key, day_data)
            return day_data

        sm_all = np.concatenate(sm_parts)
        sec_all = np.concatenate(sec_parts)
        lat_all = np.concatenate(lat_parts)
        lon_all = np.concatenate(lon_parts)
        day_data = (sm_all, sec_all, lat_all, lon_all)
        self._put_cache(day_key, day_data)
        return day_data

    def read_range(
        self,
        start_datetime: dt.datetime | str,
        end_datetime: dt.datetime | str,
    ) -> ObservationArrayResult:
        """
        Return all observations in [start_datetime, end_datetime] as arrays:
          (soil_moisture, observation_time, lat, lon)
        where observation_time is datetime64[s].
        """
        start = _to_datetime(start_datetime)
        end = _to_datetime(end_datetime)
        if end < start:
            raise ValueError("end_datetime must be >= start_datetime")

        start_sec = _to_epoch_seconds(start)
        end_sec = _to_epoch_seconds(end)
        start_day = int(np.floor(start_sec / SEC_PER_DAY))
        end_day = int(np.floor(end_sec / SEC_PER_DAY))
        sm_parts: list[np.ndarray] = []
        sec_parts: list[np.ndarray] = []
        lat_parts: list[np.ndarray] = []
        lon_parts: list[np.ndarray] = []

        for day in range(start_day, end_day + 1):
            day_data = self._load_day(day)
            if day_data[0].size == 0:
                continue

            mask = (day_data[1] >= start_sec) & (day_data[1] <= end_sec)
            if not np.any(mask):
                continue
            sm_parts.append(day_data[0][mask])
            sec_parts.append(day_data[1][mask])
            lat_parts.append(day_data[2][mask])
            lon_parts.append(day_data[3][mask])

        if not sm_parts:
            return (
                np.empty(0, dtype=np.float32),
                np.empty(0, dtype="datetime64[s]"),
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
            obs_sec[order].astype("datetime64[s]", copy=False),
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


class AmsrUpsampledL3ObservationReader:
    """
    Compatibility reader for:
      - averaged 0.5-degree output (soil_moisture + observation_count),
      - or legacy slot-based upsample/merged-style NetCDF files.
    """

    def __init__(
        self,
        path: str | Path,
        max_cache_days: int | None = None,
        interval_hours: float | None = None,
        window_hours: float | None = None,
    ) -> None:
        first = _list_input_files(Path(path))[0]
        if _is_averaged_file(first):
            if max_cache_days is None:
                if interval_hours is None:
                    max_cache_days = 8
                else:
                    max_cache_days = AmsrL3ObservationReader.suggest_cache_days(
                        interval_hours=interval_hours,
                        window_hours=window_hours,
                    )
            self._inner = Amsr0p5AveragedReader(Path(path), max_cache_days=max_cache_days)
        else:
            self._inner = AmsrL3ObservationReader(
                Path(path),
                max_cache_days=max_cache_days,
                interval_hours=interval_hours,
                window_hours=window_hours,
            )

    def read_range(
        self,
        start_datetime: dt.datetime | str,
        end_datetime: dt.datetime | str,
    ) -> ObservationArrayResult:
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
    """
    Return all upsampled/averaged AMSR observations in [start_datetime, end_datetime] as arrays:
      (soil_moisture, observation_time, lat, lon)
      where observation_time is datetime64[s].
    """
    if start_datetime is None or end_datetime is None:
        raise ValueError("start_datetime and end_datetime are required")
    reader = AmsrUpsampledL3ObservationReader(
        path,
        max_cache_days=max_cache_days,
        interval_hours=interval_hours,
        window_hours=window_hours,
    )
    return reader.read_range(start_datetime, end_datetime)


def _parse_args() -> "argparse.Namespace":
    import argparse

    parser = argparse.ArgumentParser(
        description="Read upsampled AMSR (0.5-degree) observations in datetime range."
    )
    parser.add_argument(
        "path",
        type=Path,
        nargs="?",
        default=Path("processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc"),
        help=(
            "Upsampled NetCDF file or directory (default: processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc). "
            "Legacy slot-based files are also supported."
        ),
    )
    parser.add_argument(
        "start",
        type=str,
        help="Start datetime (ISO format).",
    )
    parser.add_argument(
        "end",
        type=str,
        help="End datetime (ISO format).",
    )
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
