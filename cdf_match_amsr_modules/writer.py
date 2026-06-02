from __future__ import annotations

from pathlib import Path

from netCDF4 import Dataset
import numpy as np


def create_output_like(
    path: Path,
    src: Dataset,
    slot_pairs: list[tuple[str, str]],
    overwrite: bool,
) -> Dataset:
    if path.exists():
        if overwrite:
            path.unlink()
        else:
            raise FileExistsError(f"{path} already exists. Use --overwrite.")
    path.parent.mkdir(parents=True, exist_ok=True)

    ds = Dataset(path, "w", format="NETCDF4")
    for name, dim in src.dimensions.items():
        ds.createDimension(name, None if dim.isunlimited() else len(dim))

    for coord in ("time", "lat", "lon"):
        src_var = src.variables[coord]
        dst = ds.createVariable(coord, src_var.datatype, src_var.dimensions)
        dst[:] = src_var[:]
        for key, value in src_var.__dict__.items():
            if key != "_FillValue":
                setattr(dst, key, value)

    for sm_name, tm_name in slot_pairs:
        src_sm = src.variables[sm_name]
        src_tm = src.variables[tm_name]
        chunks = src_sm.chunking()
        kwargs = {"zlib": True, "complevel": 4, "fill_value": np.float32(np.nan)}
        if isinstance(chunks, list):
            kwargs["chunksizes"] = tuple(chunks)
        sm = ds.createVariable(sm_name, "f4", src_sm.dimensions, **kwargs)
        sm.units = "m3 m-3"
        sm.long_name = "CDF-matched volumetric soil moisture"
        sm.comment = "AMSR input was converted from percent to m3/m3 before grid-wise CDF matching."

        tm_chunks = src_tm.chunking()
        tm_kwargs = {"zlib": True, "complevel": 4, "fill_value": np.float32(np.nan)}
        if isinstance(tm_chunks, list):
            tm_kwargs["chunksizes"] = tuple(tm_chunks)
        tm = ds.createVariable(tm_name, "f4", src_tm.dimensions, **tm_kwargs)
        tm[:] = src_tm[:]
        for key, value in src_tm.__dict__.items():
            if key != "_FillValue":
                setattr(tm, key, value)

    ds.title = "Grid-wise CDF-matched AMSR soil moisture"
    ds.source = "AMSR CDF matched to reference data independently at each grid cell"
    ds.soil_moisture_units = "m3 m-3"
    return ds
