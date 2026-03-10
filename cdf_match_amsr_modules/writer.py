from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from netCDF4 import Dataset
import numpy as np


@dataclass
class AveragedHandle:
    ds: Dataset
    var_time: object
    var_soil: object
    var_count: object
    var_obs_time_min: object
    next_index: int = 0


@dataclass
class Handle:
    ds: Dataset
    var_time: object
    lat_var: object
    lon_var: object
    sm_vars: dict[str, object] = field(default_factory=dict)
    tm_vars: dict[str, object] = field(default_factory=dict)
    max_slots: dict[str, int] = field(default_factory=lambda: {"A": 0, "D": 0})
    next_index: int = 0


class AveragedProcessedWriter:
    def __init__(
        self,
        output_dir: Path,
        output_file: Path,
        mode: str,
        overwrite: bool,
        nlat: int,
        nlon: int,
        lat: np.ndarray,
        lon: np.ndarray,
        time_units: str,
    ) -> None:
        self.output_dir = output_dir
        self.output_file = output_file
        self.mode = mode
        self.overwrite = overwrite
        self.nlat = nlat
        self.nlon = nlon
        self.lat = lat
        self.lon = lon
        self.time_units = time_units
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.handles: dict[Path, AveragedHandle] = {}

    def _create_handle(self, path: Path) -> AveragedHandle:
        if path.exists():
            if self.overwrite:
                path.unlink()
            else:
                raise FileExistsError(f"{path} already exists. Use --overwrite.")
        ds = Dataset(path, "w", format="NETCDF4")
        ds.createDimension("time", None)
        ds.createDimension("lat", self.nlat)
        ds.createDimension("lon", self.nlon)

        var_time = ds.createVariable("time", "i4", ("time",))
        var_lat = ds.createVariable("lat", "f4", ("lat",))
        var_lon = ds.createVariable("lon", "f4", ("lon",))
        var_soil = ds.createVariable(
            "soil_moisture",
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, self.nlat), min(360, self.nlon)),
            fill_value=np.float32(np.nan),
        )
        var_count = ds.createVariable(
            "observation_count",
            "i4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, self.nlat), min(360, self.nlon)),
            fill_value=np.int32(0),
        )
        var_obs_time_min = ds.createVariable(
            "observation_time_min",
            "f4",
            ("time", "lat", "lon"),
            zlib=True,
            complevel=4,
            chunksizes=(1, min(180, self.nlat), min(360, self.nlon)),
            fill_value=np.float32(np.nan),
        )

        var_lat[:] = self.lat
        var_lon[:] = self.lon

        var_time.units = self.time_units
        var_time.calendar = "standard"
        var_time.long_name = "Time"
        var_lat.units = "degrees_north"
        var_lon.units = "degrees_east"
        var_soil.units = "%"
        var_soil.long_name = "daily mean soil moisture (0.5 degree)"
        var_count.units = "1"
        var_count.long_name = "number of merged observations used in daily average"
        var_obs_time_min.units = "minutes since 00:00 UTC of corresponding time day"
        var_obs_time_min.long_name = "daily mean observation time (0.5 degree)"
        var_obs_time_min.comment = "Copied from input; CDF matching does not alter observation timing."

        ds.title = "CDF-matched 0.5-degree daily soil moisture"
        ds.source = "CDF matching against provided reference data"
        ds.note = "Output uses the same variable names as merged 0.5-degree averaged files."

        return AveragedHandle(
            ds=ds,
            var_time=var_time,
            var_soil=var_soil,
            var_count=var_count,
            var_obs_time_min=var_obs_time_min,
        )

    def _get_handle(self, year: int) -> AveragedHandle:
        if self.mode == "single":
            path = self.output_dir / self.output_file
        else:
            suffix = self.output_file.suffix if self.output_file.suffix else ".nc"
            stem = self.output_file.stem or "AMSR_SMC_daily_0p5deg_cdf"
            path = self.output_dir / f"{stem}_{year:04d}{suffix}"
        if path not in self.handles:
            self.handles[path] = self._create_handle(path)
        return self.handles[path]

    def write_day(
        self,
        day_value: int,
        year: int,
        soil_moisture_2d: np.ndarray,
        observation_count_2d: np.ndarray,
        observation_time_min_2d: np.ndarray,
    ) -> None:
        handle = self._get_handle(year)
        idx = handle.next_index
        handle.var_time[idx] = int(day_value)
        handle.next_index += 1
        handle.var_soil[idx, :, :] = soil_moisture_2d
        handle.var_count[idx, :, :] = observation_count_2d
        handle.var_obs_time_min[idx, :, :] = observation_time_min_2d

    def close(self) -> None:
        for handle in self.handles.values():
            handle.ds.close()


class YearlyProcessedWriter:
    def __init__(
        self,
        output_dir: Path,
        output_file: Path,
        mode: str,
        overwrite: bool,
        nlat: int,
        nlon: int,
        lat: np.ndarray,
        lon: np.ndarray,
        time_units: str,
    ) -> None:
        self.output_dir = output_dir
        self.output_file = output_file
        self.mode = mode
        self.overwrite = overwrite
        self.nlat = nlat
        self.nlon = nlon
        self.lat = lat
        self.lon = lon
        self.time_units = time_units
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.handles: dict[Path, Handle] = {}

    def _create_handle(self, path: Path) -> Handle:
        if path.exists():
            if self.overwrite:
                path.unlink()
            else:
                raise FileExistsError(f"{path} already exists. Use --overwrite.")
        ds = Dataset(path, "w", format="NETCDF4")
        ds.createDimension("time", None)
        ds.createDimension("lat", self.nlat)
        ds.createDimension("lon", self.nlon)

        var_time = ds.createVariable("time", "i4", ("time",))
        var_lat = ds.createVariable("lat", "f4", ("lat",))
        var_lon = ds.createVariable("lon", "f4", ("lon",))
        var_lat[:] = self.lat
        var_lon[:] = self.lon

        var_time.units = self.time_units
        var_time.calendar = "standard"
        var_time.long_name = "Time"
        var_lat.units = "degrees_north"
        var_lon.units = "degrees_east"

        ds.title = "CDF-matched merged AMSR-L3 daily Soil Moisture"
        ds.source = "CDF matching against provided reference data"
        ds.note = "A/D observations are separated; duplicate cells are stored in slot variables."
        ds.max_slots_A = 0
        ds.max_slots_D = 0

        return Handle(ds=ds, var_time=var_time, lat_var=var_lat, lon_var=var_lon)

    def _get_handle(self, year: int) -> Handle:
        if self.mode == "single":
            path = self.output_dir / self.output_file
        else:
            path = self.output_dir / f"AMSR_SMC_daily_{year:04d}.nc"
        if path not in self.handles:
            self.handles[path] = self._create_handle(path)
        return self.handles[path]

    @staticmethod
    def _copy_attrs(src_var, dst_var) -> None:
        for key, value in src_var.__dict__.items():
            if key.startswith("_"):
                continue
            try:
                setattr(dst_var, key, value)
            except Exception:
                continue

    def _ensure_slot_variables(self, handle: Handle, sm_name: str, tm_name: str, ds_src: Dataset) -> tuple[object, object]:
        if sm_name not in handle.sm_vars:
            src_sm = ds_src.variables[sm_name]
            src_tm = ds_src.variables[tm_name]

            var_sm = handle.ds.createVariable(
                sm_name,
                "f4",
                ("time", "lat", "lon"),
                zlib=True,
                complevel=4,
                chunksizes=(1, min(180, self.nlat), min(360, self.nlon)),
                fill_value=np.float32(np.nan),
            )
            var_tm = handle.ds.createVariable(
                tm_name,
                "f4",
                ("time", "lat", "lon"),
                zlib=True,
                complevel=4,
                chunksizes=(1, min(180, self.nlat), min(360, self.nlon)),
                fill_value=np.float32(np.nan),
            )
            self._copy_attrs(src_sm, var_sm)
            self._copy_attrs(src_tm, var_tm)
            handle.sm_vars[sm_name] = var_sm
            handle.tm_vars[tm_name] = var_tm

            if sm_name.startswith("soil_moisture_A_"):
                handle.max_slots["A"] = max(handle.max_slots["A"], int(sm_name[-2:]))
            elif sm_name.startswith("soil_moisture_D_"):
                handle.max_slots["D"] = max(handle.max_slots["D"], int(sm_name[-2:]))

        return handle.sm_vars[sm_name], handle.tm_vars[tm_name]

    def write_day(
        self,
        day_value: int,
        year: int,
        slot_payload: list[tuple[str, str, np.ndarray, np.ndarray]],
        src_ds: Dataset,
    ) -> None:
        handle = self._get_handle(year)
        idx = handle.next_index
        handle.var_time[idx] = int(day_value)
        handle.next_index += 1

        for sm_name, tm_name, sm_arr, tm_arr in slot_payload:
            var_sm, var_tm = self._ensure_slot_variables(handle, sm_name, tm_name, src_ds)
            var_sm[idx, :, :] = sm_arr
            var_tm[idx, :, :] = tm_arr

    def close(self) -> None:
        for handle in self.handles.values():
            handle.ds.max_slots_A = int(handle.max_slots["A"])
            handle.ds.max_slots_D = int(handle.max_slots["D"])
            handle.ds.close()
