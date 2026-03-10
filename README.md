# AMSR Download / Merge / Query / Bias Correction

This repository provides a workflow to prepare AMSR-L3 products, generate 0.5-degree aggregated files, run CDF matching, and query by datetime range.

---

## 1. Requirements

- `bash`
- `lftp` (with SFTP support)
- Python 3.10+
- Python packages:
  - `numpy`
  - `netCDF4`
  - `tqdm` (used for progress output during upscaling)

Install Python packages:

```bash
python -m pip install numpy netCDF4 tqdm
```

---

## 2. Credentials file (required)

`AMSR_download.sh` reads credentials from `~/.gportal_sftp.env`.

Create the file as:

```bash
cat > ~/.gportal_sftp.env << 'EOF'
GPORTAL_USER="your_gportal_username"
GPORTAL_PASS="your_gportal_password"
EOF
chmod 600 ~/.gportal_sftp.env
```

Use the exact `KEY="value"` format.

---

## 3. Download AMSR data

```bash
cd /path/to/AMSR
bash AMSR_download.sh
```

Notes:

- Product settings (AMSR-E/AMSR2, L2/L3) are selected by editing the configuration block at the top of `AMSR_download.sh`.
- Current default product is AMSR2 L3 (`GCOM-W.AMSR2_L3.SMC_10_3`).
- Missing remote directories are written to `missing_dirs.log`.

Typical outputs:

- `./download/AQUA.AMSR-E_AMSR2Format.L3.SMC_10.8/...`
- `./download/GCOM-W.AMSR2_L3.SMC_10_3/...`

---

## 4. Merge to NetCDF

Merge daily L3 HDF5 files (`*_01D_*.h5`) into NetCDF:

```bash
cd /path/to/AMSR
python merge_amsr_l3_daily.py \
  --output-mode yearly \
  --output-dir processed/l3_daily \
  --overwrite
```

Output format:

- Dimensions: `(time, lat, lon)`
- Same-day duplicates on the same grid are not averaged.
- Orbit-direction slots are separated:
  - `soil_moisture_A_01`, `soil_moisture_A_02`, ...
  - `soil_moisture_D_01`, `soil_moisture_D_02`, ...
- Associated times:
  - `observation_time_min_A_XX`, `observation_time_min_D_XX`
- Day-shift correction:
  - `>= 1440` minutes is moved to next day
  - `< 0` minutes is moved to previous day

---

## 5. 0.5-degree upscaling

Upscale merged daily AMSR data to a 0.5-degree daily mean grid:

```bash
cd /path/to/AMSR
python upscale_amsr_l3_0p5deg.py \
  processed/l3_daily \
  processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc \
  --overwrite
```

Run with default paths:

```bash
cd /path/to/AMSR
python upscale_amsr_l3_0p5deg.py
```

Useful options:

- `input_path`: merged NetCDF file or directory
- `output_path`: output NetCDF path
- `--overwrite`: overwrite existing output
- `--workers N`: set parallel worker count for daily aggregation

Example with explicit options:

```bash
python upscale_amsr_l3_0p5deg.py \
  processed/l3_daily \
  processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc \
  --overwrite \
  --workers 8
```

Output variables:

- `soil_moisture(time, lat, lon)`
- `observation_count(time, lat, lon)`
- `observation_time_min(time, lat, lon)` (daily mean of in-day observation minutes)

`time` is daily (days since 1970-01-01).

---

## 6. CDF matching for AMSR-L3 (bias correction)

`cdf_match_amsr.py` performs bias correction by CDF matching.

### CLI

```bash
python cdf_match_amsr.py <amsr_input_path> \
  --reference <ref_nc_or_glob> \
  [--reference <another_ref_nc>] ... \
  [--reference-dict <json_file>] \
  [--reference-time <name>] \
  [--reference-lat <name>] \
  [--reference-lon <name>] \
  [--reference-var <name>] \
  [--reference-depth <depth_index>] \
  [--start <YYYY-MM-DDTHH:MM:SS>] \
  [--end <YYYY-MM-DDTHH:MM:SS>] \
  [--step-hours <hours>] \
  [--window-hours <hours>] \
  [--output-dir processed/l3_daily_cdf] \
  [--output-mode yearly|single] \
  [--output-file AMSR_SMC_daily_cdf.nc] \
  [--overwrite]
```

Example with yearly SoilMoistV references:

```bash
python cdf_match_amsr.py processed/l3_daily \
  --reference /path/to/reference/2016/SoilMoistV.nc \
  --reference /path/to/reference/2017/SoilMoistV.nc \
  --reference /path/to/reference/2018/SoilMoistV.nc \
  --start 2016-01-01T00:00:00 \
  --end 2018-12-31T23:59:59 \
  --output-dir processed/l3_daily_cdf \
  --output-mode yearly \
  --overwrite
```

0.5-degree daily file example (single and multi-year):

```bash
python cdf_match_amsr.py processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc \
  --reference-dict reference_paths_template.json \
  --start 2013-01-01T00:00:00 \
  --end 2022-12-31T23:59:59 \
  --output-dir processed/l3_daily_0p5 \
  --output-file AMSR_SMC_daily_0p5deg_cdf.nc \
  --overwrite
```

For a single output file over multiple years:

```bash
python cdf_match_amsr.py processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc \
  --reference-dict reference_paths_template.json \
  --start 2013-01-01T00:00:00 \
  --end 2022-12-31T23:59:59 \
  --output-dir processed/l3_daily_0p5 \
  --output-mode single \
  --output-file AMSR_SMC_daily_0p5deg_cdf_2013_2022.nc \
  --overwrite
```

Reference input supports:

- repeated `--reference`
- `--reference-dict` as a file path or JSON string
- yearly dictionary or `path_template` format

### Notes

- `--reference-depth` defaults to `0` if omitted.
- CDF mapping is built from AMSR-reference pairs inside the reference time span.
- For AMSR periods outside reference coverage, the same mapping table is applied by linear CDF interpolation.
- For 30-minute references starting at `00:30`, matching still works by same-day pairing at the 0.5-degree daily level.
- For 0.5-degree averaged input, `observation_time_min` is preserved in CDF output (soil moisture only is bias-corrected).

---

## 7. Query observations by datetime range

`query_amsr_l3_0p5deg.py` reads both averaged 0.5-degree files and legacy slot-style files. `amsr_l3_query.py` reads only legacy slot-style files.

### CLI examples

```bash
python query_amsr_l3_0p5deg.py \
  processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc \
  2012-07-03T00:00:00 \
  2012-07-03T23:59:59 \
  --interval-hours 3 \
  --limit 10
```

```bash
python amsr_l3_query.py \
  processed/l3_daily \
  2012-07-03T00:00:00 2012-07-03T23:59:59 \
  --interval-hours 3 \
  --limit 10
```

### Python usage

```python
from query_amsr_l3_0p5deg import AmsrUpsampledL3ObservationReader

reader = AmsrUpsampledL3ObservationReader(
    "processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc",
)
sm, obs_time, lat, lon = reader.read_range(
    "2012-07-03T00:00:00",
    "2012-07-03T23:59:59",
)
print(len(sm), sm[0], obs_time[0], lat[0], lon[0])

import datetime as dt
for ws, we, sm, obs_time, lat, lon in reader.iterate_windows(
    "2012-07-01T00:00:00",
    "2012-07-03T23:59:59",
    step=dt.timedelta(hours=24),
):
    print(ws, we, len(sm))
```

```python
from amsr_l3_query import AmsrL3ObservationReader

legacy = AmsrL3ObservationReader("processed/l3_daily")
sm, obs_time, lat, lon = legacy.read_range(
    "2012-07-03T00:00:00",
    "2012-07-03T05:59:59",
)
```
