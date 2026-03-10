#!/usr/bin/env python3
"""
Merge AMSR-E / AMSR2 L3 SMC daily HDF5 files into NetCDF.

Output layout:
  - dimensions: (time, lat, lon)
  - variables:
      * time(time), lat(lat), lon(lon)
      * soil_moisture_A_XX(time, lat, lon), observation_time_min_A_XX(time, lat, lon)
      * soil_moisture_D_XX(time, lat, lon), observation_time_min_D_XX(time, lat, lon)
        where XX is slot number (01, 02, ...)

`observation_time_min` is minutes since 00:00 UTC of each `time` day.
If a source pixel has time >= 1440 min, it is moved to the next day.
If a source pixel has time < 0 min, it is moved to the previous day.
No averaging is applied. If multiple observations fall on the same
(time, lat, lon, direction), they are stored in additional slot variables.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import numpy as np
from netCDF4 import Dataset

DATE_PATTERN = re.compile(r"_(\d{8})_01D_")


def parse_date(text: str) -> dt.date:
    return dt.datetime.strptime(text, "%Y-%m-%d").date()


def discover_daily_files(
    input_dirs: List[Path],
    start_date: dt.date | None,
    end_date: dt.date | None,
) -> Dict[dt.date, List[Path]]:
    grouped: Dict[dt.date, List[Path]] = defaultdict(list)
    for root in input_dirs:
        for path in sorted(root.rglob("*.h5")):
            name = path.name
            if "_01D_" not in name:
                continue
            m = DATE_PATTERN.search(name)
            if not m:
                continue
            nominal_date = dt.datetime.strptime(m.group(1), "%Y%m%d").date()
            if start_date and nominal_date < start_date:
                continue
            if end_date and nominal_date > end_date:
                continue
            grouped[nominal_date].append(path)
    return grouped


def inspect_grid_shape(sample_file: Path) -> tuple[int, int]:
    with Dataset(sample_file) as ds:
        g = ds.variables["Geophysical Data"]
        if g.ndim == 3:
            return int(g.shape[0]), int(g.shape[1])
        if g.ndim == 2:
            return int(g.shape[0]), int(g.shape[1])
        raise ValueError(f"Unexpected Geophysical Data shape in {sample_file}: {g.shape}")


def build_lat_lon(
    nlat: int,
    nlon: int,
    resolution_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    # AMSR EQR grid is typically north -> south for latitude.
    lat = (90.0 - resolution_deg / 2.0) - np.arange(nlat, dtype=np.float32) * resolution_deg
    lon = (resolution_deg / 2.0) + np.arange(nlon, dtype=np.float32) * resolution_deg
    return lat.astype(np.float32), lon.astype(np.float32)


def parse_orbit_direction(name: str) -> str:
    if "_EQMA_" in name:
        return "A"
    if "_EQMD_" in name:
        return "D"
    raise ValueError(f"Cannot parse orbit direction from file name: {name}")


@dataclass
class SlotAccumulator:
    sm: np.ndarray
    obs_time: np.ndarray
    filled: np.ndarray

    @classmethod
    def create(cls, shape: tuple[int, int]) -> "SlotAccumulator":
        return cls(
            sm=np.full(shape, np.nan, dtype=np.float32),
            obs_time=np.full(shape, np.nan, dtype=np.float32),
            filled=np.zeros(shape, dtype=bool),
        )


@dataclass
class DayAccumulator:
    shape: tuple[int, int]
    slots_by_direction: dict[str, list[SlotAccumulator]] = field(
        default_factory=lambda: {"A": [], "D": []}
    )

    @classmethod
    def create(cls, shape: tuple[int, int]) -> "DayAccumulator":
        return cls(shape=shape)

    def add(self, direction: str, mask: np.ndarray, sm: np.ndarray, obs_time: np.ndarray) -> None:
        if direction not in self.slots_by_direction:
            self.slots_by_direction[direction] = []
        slots = self.slots_by_direction[direction]
        remaining = mask.copy()
        if not np.any(remaining):
            return

        for slot in slots:
            can_place = remaining & (~slot.filled)
            if np.any(can_place):
                slot.sm[can_place] = sm[can_place]
                slot.obs_time[can_place] = obs_time[can_place]
                slot.filled[can_place] = True
                remaining[can_place] = False
            if not np.any(remaining):
                return

        new_slot = SlotAccumulator.create(self.shape)
        new_slot.sm[remaining] = sm[remaining]
        new_slot.obs_time[remaining] = obs_time[remaining]
        new_slot.filled[remaining] = True
        slots.append(new_slot)

    def finalize(self) -> dict[str, list[tuple[np.ndarray, np.ndarray]]]:
        out: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {}
        for direction, slots in self.slots_by_direction.items():
            out[direction] = [(slot.sm, slot.obs_time) for slot in slots]
        return out


@dataclass
class OutputHandle:
    dataset: Dataset
    var_time: any
    nlat: int
    nlon: int
    next_index: int = 0
    sm_vars: dict[str, any] = field(default_factory=dict)
    time_vars: dict[str, any] = field(default_factory=dict)
    max_slots_by_direction: dict[str, int] = field(default_factory=lambda: {"A": 0, "D": 0})

    def _slot_key(self, direction: str, slot_id: int) -> str:
        return f"{direction}_{slot_id:02d}"

    def _ensure_slot_variables(self, direction: str, slot_id: int) -> tuple[any, any]:
        key = self._slot_key(direction, slot_id)
        if key in self.sm_vars:
            return self.sm_vars[key], self.time_vars[key]

        var_sm_name = f"soil_moisture_{direction}_{slot_id:02d}"
        var_tm_name = f"observation_time_min_{direction}_{slot_id:02d}"
        var_sm = self.dataset.createVariable(
            var_sm_name,
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, self.nlat), min(360, self.nlon)),
            fill_value=np.float32(np.nan),
        )
        var_tm = self.dataset.createVariable(
            var_tm_name,
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, self.nlat), min(360, self.nlon)),
            fill_value=np.float32(np.nan),
        )

        var_sm.units = "%"
        var_sm.long_name = f"Soil Moisture Content ({direction}, slot {slot_id})"
        var_sm.comment = "No averaging; duplicate observations are stored in additional slots."
        var_tm.units = "minutes since 00:00 UTC of corresponding time day"
        var_tm.long_name = f"Observation Time ({direction}, slot {slot_id})"

        self.sm_vars[key] = var_sm
        self.time_vars[key] = var_tm
        self.max_slots_by_direction[direction] = max(self.max_slots_by_direction.get(direction, 0), slot_id)
        return var_sm, var_tm

    def write_day(self, day: dt.date, day_payload: dict[str, list[tuple[np.ndarray, np.ndarray]]]) -> bool:
        idx = self.next_index
        epoch = dt.date(1970, 1, 1)
        self.var_time[idx] = (day - epoch).days

        for direction, slots in day_payload.items():
            for slot_id, (sm, obs_time) in enumerate(slots, start=1):
                var_sm, var_tm = self._ensure_slot_variables(direction, slot_id)
                var_sm[idx, :, :] = sm
                var_tm[idx, :, :] = obs_time

        self.next_index += 1
        return True

    def close(self) -> None:
        self.dataset.setncattr("max_slots_A", int(self.max_slots_by_direction.get("A", 0)))
        self.dataset.setncattr("max_slots_D", int(self.max_slots_by_direction.get("D", 0)))
        self.dataset.close()


class OutputManager:
    def __init__(
        self,
        mode: str,
        output_dir: Path,
        output_file: str | None,
        overwrite: bool,
        nlat: int,
        nlon: int,
        lat: np.ndarray,
        lon: np.ndarray,
    ) -> None:
        self.mode = mode
        self.output_dir = output_dir
        self.output_file = output_file
        self.overwrite = overwrite
        self.nlat = nlat
        self.nlon = nlon
        self.lat = lat
        self.lon = lon
        self.handles: dict[str, OutputHandle] = {}

        self.output_dir.mkdir(parents=True, exist_ok=True)
        if self.mode == "single":
            filename = self.output_file or "AMSR_SMC_daily_merged.nc"
            path = self.output_dir / filename
            self.handles["single"] = self._create_handle(path)

    def _create_handle(self, output_path: Path) -> OutputHandle:
        if output_path.exists():
            if self.overwrite:
                output_path.unlink()
            else:
                raise FileExistsError(f"{output_path} already exists. Use --overwrite.")

        ds = Dataset(output_path, "w", format="NETCDF4")
        ds.createDimension("time", None)
        ds.createDimension("lat", self.nlat)
        ds.createDimension("lon", self.nlon)

        var_time = ds.createVariable("time", "i4", ("time",))
        var_lat = ds.createVariable("lat", "f4", ("lat",))
        var_lon = ds.createVariable("lon", "f4", ("lon",))
        var_lat[:] = self.lat
        var_lon[:] = self.lon

        var_time.units = "days since 1970-01-01 00:00:00"
        var_time.calendar = "standard"
        var_time.long_name = "Time"
        var_lat.units = "degrees_north"
        var_lon.units = "degrees_east"

        ds.title = "Merged AMSR-E/AMSR2 daily L3 Soil Moisture Content"
        ds.history = (
            f"Created {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')} "
            f"by merge_amsr_l3_daily.py"
        )
        ds.source = "AQUA AMSR-E AMSR2Format L3 SMC + GCOM-W AMSR2 L3 SMC"
        ds.note = "A/D observations are separated; duplicate cells are stored in slot variables."

        return OutputHandle(dataset=ds, var_time=var_time, nlat=self.nlat, nlon=self.nlon)

    def _year_key(self, year: int) -> str:
        return f"year_{year:04d}"

    def _get_year_handle(self, year: int) -> OutputHandle:
        key = self._year_key(year)
        if key not in self.handles:
            output_path = self.output_dir / f"AMSR_SMC_daily_{year}.nc"
            self.handles[key] = self._create_handle(output_path)
        return self.handles[key]

    def write_day(self, day: dt.date, day_payload: dict[str, list[tuple[np.ndarray, np.ndarray]]]) -> bool:
        if self.mode == "single":
            return self.handles["single"].write_day(day, day_payload)
        handle = self._get_year_handle(day.year)
        return handle.write_day(day, day_payload)

    def close(self) -> None:
        for handle in self.handles.values():
            handle.close()


def extract_scale(dataset: any, default: float = 1.0) -> float:
    if "SCALE FACTOR" not in dataset.ncattrs():
        return float(default)
    raw = dataset.getncattr("SCALE FACTOR")
    val = np.array(raw).reshape(-1)[0]
    return float(val)


def process_file(
    file_path: Path,
    nominal_date: dt.date,
    accumulators: dict[dt.date, DayAccumulator],
    grid_shape: tuple[int, int],
) -> None:
    orbit_direction = parse_orbit_direction(file_path.name)

    with Dataset(file_path) as ds:
        ds.set_auto_maskandscale(False)
        ds_sm = ds.variables["Geophysical Data"]
        ds_tm = ds.variables["Time Information"]

        sm_raw = ds_sm[:, :, 0] if ds_sm.ndim == 3 else ds_sm[:, :]
        tm_raw = ds_tm[:, :]
        if sm_raw.shape != grid_shape or tm_raw.shape != grid_shape:
            raise ValueError(
                f"Grid shape mismatch in {file_path}: "
                f"sm={sm_raw.shape}, tm={tm_raw.shape}, expected={grid_shape}"
            )

        sm_scale = extract_scale(ds_sm, default=0.1)
        tm_scale = extract_scale(ds_tm, default=1.0)

        sm_raw = sm_raw.astype(np.int32, copy=False)
        tm_raw = tm_raw.astype(np.int32, copy=False)

        # Non-negative SMC indicates valid retrievals in this product.
        # Valid observation times can be negative (previous-day observations),
        # so we only remove obvious sentinel values.
        valid = (sm_raw >= 0) & (tm_raw > -30000)
        if not np.any(valid):
            return

        sm = sm_raw.astype(np.float32) * sm_scale
        time_min = np.rint(tm_raw.astype(np.float32) * tm_scale).astype(np.int32)

        day_offset = np.floor_divide(time_min, 1440)
        minute_in_day = np.mod(time_min, 1440).astype(np.float32)

        for off in np.unique(day_offset[valid]):
            mask = valid & (day_offset == off)
            if not np.any(mask):
                continue
            target_day = nominal_date + dt.timedelta(days=int(off))
            acc = accumulators.get(target_day)
            if acc is None:
                acc = DayAccumulator.create(grid_shape)
                accumulators[target_day] = acc

            acc.add(orbit_direction, mask, sm, minute_in_day)


def flush_ready_days(
    accumulators: dict[dt.date, DayAccumulator],
    cutoff_day: dt.date,
    writer: OutputManager,
) -> int:
    ready_days = sorted([d for d in accumulators.keys() if d <= cutoff_day])
    written = 0
    for day in ready_days:
        day_payload = accumulators[day].finalize()
        if writer.write_day(day, day_payload):
            written += 1
        del accumulators[day]
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge AMSR-E/AMSR2 L3 SMC daily files into daily NetCDF(s)."
    )
    parser.add_argument(
        "--input-dir",
        action="append",
        type=Path,
        default=None,
        help="Input root directory (can be specified multiple times).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/data02/shiojiri/DATA/AMSR/processed/l3_daily"),
        help="Output directory.",
    )
    parser.add_argument(
        "--output-mode",
        choices=["single", "yearly"],
        default="yearly",
        help="Write one file (`single`) or one file per year (`yearly`).",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        default="AMSR_SMC_daily_merged.nc",
        help="Output filename when --output-mode=single.",
    )
    parser.add_argument(
        "--start-date",
        type=parse_date,
        default=None,
        help="Filter by nominal date (inclusive), format: YYYY-MM-DD.",
    )
    parser.add_argument(
        "--end-date",
        type=parse_date,
        default=None,
        help="Filter by nominal date (inclusive), format: YYYY-MM-DD.",
    )
    parser.add_argument(
        "--resolution-deg",
        type=float,
        default=0.1,
        help="Grid resolution in degree.",
    )
    parser.add_argument(
        "--max-backward-days",
        type=int,
        default=1,
        help=(
            "Maximum day shift expected in backward direction. "
            "If non-zero, ready-day flush is delayed accordingly."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output file(s).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_dirs = (
        [Path(p) for p in args.input_dir]
        if args.input_dir
        else [
            Path("/data02/shiojiri/DATA/AMSR/download/AQUA.AMSR-E_AMSR2Format.L3.SMC_10.8"),
            Path("/data02/shiojiri/DATA/AMSR/download/GCOM-W.AMSR2_L3.SMC_10_3"),
        ]
    )
    for path in input_dirs:
        if not path.exists():
            raise FileNotFoundError(f"Input directory not found: {path}")

    grouped = discover_daily_files(input_dirs, args.start_date, args.end_date)
    if not grouped:
        print("No *_01D_*.h5 files found in the selected range.")
        return

    nominal_days = sorted(grouped.keys())
    sample_file = grouped[nominal_days[0]][0]
    nlat, nlon = inspect_grid_shape(sample_file)
    lat, lon = build_lat_lon(nlat, nlon, args.resolution_deg)
    grid_shape = (nlat, nlon)

    writer = OutputManager(
        mode=args.output_mode,
        output_dir=args.output_dir,
        output_file=args.output_file,
        overwrite=args.overwrite,
        nlat=nlat,
        nlon=nlon,
        lat=lat,
        lon=lon,
    )

    accumulators: dict[dt.date, DayAccumulator] = {}
    total_files = sum(len(v) for v in grouped.values())
    done_files = 0
    written_days = 0

    print(f"Found {total_files} daily files across {len(nominal_days)} nominal days.")
    print(f"Grid shape: lat={nlat}, lon={nlon}")

    for i, nominal_day in enumerate(nominal_days, start=1):
        for fp in grouped[nominal_day]:
            process_file(fp, nominal_day, accumulators, grid_shape)
            done_files += 1

        cutoff_day = nominal_day - dt.timedelta(days=args.max_backward_days)
        written_days += flush_ready_days(accumulators, cutoff_day, writer)

        if (i % 30 == 0) or (i == len(nominal_days)):
            print(
                f"[{i}/{len(nominal_days)} days] processed files={done_files}/{total_files}, "
                f"buffered_days={len(accumulators)}, written_days={written_days}"
            )

    # Flush remaining days
    for day in sorted(accumulators.keys()):
        day_payload = accumulators[day].finalize()
        if writer.write_day(day, day_payload):
            written_days += 1
    accumulators.clear()
    writer.close()

    print("Done.")
    print(f"Output mode: {args.output_mode}")
    print(f"Output directory: {args.output_dir}")
    print(f"Days written: {written_days}")


if __name__ == "__main__":
    main()
