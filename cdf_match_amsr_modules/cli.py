from __future__ import annotations

import argparse
from pathlib import Path

from .pipeline import CDFMatchConfig, run_cdf_matching


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply CDF matching to AMSR-L3 merged NetCDF files."
    )
    parser.add_argument("input_path", type=Path, help="Input merged AMSR directory or .nc file.")
    parser.add_argument(
        "--reference",
        action="append",
        help="Reference NetCDF file path. Repeat for multiple files.",
    )
    parser.add_argument(
        "--reference-dict",
        dest="reference_dict",
        help="JSON (file path or inline JSON string) containing reference paths as a list or dict.",
    )
    parser.add_argument("--reference-time", default=None, help="Reference time variable name (default: auto).")
    parser.add_argument("--reference-lat", default=None, help="Reference lat variable name (default: auto).")
    parser.add_argument("--reference-lon", default=None, help="Reference lon variable name (default: auto).")
    parser.add_argument("--reference-var", default=None, help="Reference data variable name (default: auto).")
    parser.add_argument(
        "--reference-depth",
        type=int,
        default=0,
        help="Depth index for 4D reference variables (default: 0, i.e., first depth).",
    )
    parser.add_argument("--start", default=None, help="ISO8601 start datetime.")
    parser.add_argument("--end", default=None, help="ISO8601 end datetime.")
    parser.add_argument("--step-hours", type=float, default=24.0, help="Iteration step in hours.")
    parser.add_argument("--window-hours", type=float, default=24.0, help="Window size in hours.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("processed/l3_daily_cdf"),
        help="Output directory.",
    )
    parser.add_argument(
        "--output-mode",
        choices=("single", "yearly"),
        default="yearly",
        help="single: one file, yearly: yearly files.",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default=Path("AMSR_SMC_daily_cdf.nc"),
        help="Output file name in single mode.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    parser.add_argument(
        "--reference-cache-size",
        type=int,
        default=32,
        help="Reference timestep cache size per file.",
    )
    return parser.parse_args()


def build_config(args: argparse.Namespace) -> CDFMatchConfig:
    reference_depth = 0 if args.reference_depth is None else args.reference_depth
    return CDFMatchConfig(
        input_path=args.input_path,
        reference=args.reference,
        reference_dict=args.reference_dict,
        reference_time=args.reference_time,
        reference_lat=args.reference_lat,
        reference_lon=args.reference_lon,
        reference_var=args.reference_var,
        start=args.start,
        end=args.end,
        step_hours=args.step_hours,
        window_hours=args.window_hours,
        output_dir=args.output_dir,
        output_mode=args.output_mode,
        output_file=args.output_file,
        overwrite=args.overwrite,
        reference_cache_size=args.reference_cache_size,
        reference_depth=reference_depth,
    )


def main() -> None:
    args = parse_args()
    config = build_config(args)
    run_cdf_matching(config)
