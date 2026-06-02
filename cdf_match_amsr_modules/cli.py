from __future__ import annotations

import argparse
import os
from pathlib import Path

from .pipeline import CDFMatchConfig, run_cdf_matching


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply grid-wise CDF matching to AMSR 0.5-degree soil moisture."
    )
    parser.add_argument("input_path", type=Path, help="Input AMSR NetCDF file.")
    parser.add_argument("--reference", action="append", help="Reference NetCDF path. Repeatable.")
    parser.add_argument("--reference-dict", help="JSON file or inline JSON with reference paths.")
    parser.add_argument("--reference-time", default=None)
    parser.add_argument("--reference-lat", default=None)
    parser.add_argument("--reference-lon", default=None)
    parser.add_argument("--reference-var", default=None)
    parser.add_argument("--reference-depth", type=int, default=0)
    parser.add_argument("--start", default=None, help="ISO8601 start datetime.")
    parser.add_argument("--end", default=None, help="ISO8601 end datetime.")
    parser.add_argument("--output-dir", type=Path, default=Path("processed/l3_daily_0p5"))
    parser.add_argument("--output-file", type=Path, default=Path("AMSR_SMC_daily_0p5deg_cdf.nc"))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--block-rows",
        type=int,
        default=4,
        help="Latitude rows processed at once. Increase on high-memory nodes.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, os.cpu_count() or 1),
        help="Thread workers for per-grid CDF mapping inside each block.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_cdf_matching(
        CDFMatchConfig(
            input_path=args.input_path,
            reference=args.reference,
            reference_dict=args.reference_dict,
            reference_time=args.reference_time,
            reference_lat=args.reference_lat,
            reference_lon=args.reference_lon,
            reference_var=args.reference_var,
            start=args.start,
            end=args.end,
            output_dir=args.output_dir,
            output_file=args.output_file,
            overwrite=args.overwrite,
            reference_depth=args.reference_depth,
            block_rows=args.block_rows,
            workers=args.workers,
        )
    )
