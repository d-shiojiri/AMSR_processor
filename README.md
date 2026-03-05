# AMSR Download / Merge / Query

This repository provides a workflow to:

1. Download AMSR products from JAXA G-Portal with `AMSR_download.sh`
2. Merge daily L3 files into NetCDF with `merge_amsr_l3_daily.py`
3. Upscale merged daily data to 0.5-degree grid with `upscale_amsr_l3_0p5deg.py`
4. Query observations or averaged 0.5-degree grid data by datetime range with
   `amsr_l3_query.py` or `query_amsr_l3_0p5deg.py`

---

## 1. Requirements

- `bash`
- `lftp` (with SFTP support)
- Python 3.10+
- Python packages:
  - `numpy`
  - `netCDF4`
  - `tqdm` (for progress bar display in upscaling)

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
cd /data02/shiojiri/DATA/AMSR
bash AMSR_download.sh
```

Notes:

- Target product settings (AMSR-E/AMSR2, L2/L3) are selected by editing the config block at the top of `AMSR_download.sh`.
- Current default is AMSR2 L3 (`GCOM-W.AMSR2_L3.SMC_10_3`).
- Missing remote directories are logged to `missing_dirs.log`.

Typical output roots:

- `./download/AQUA.AMSR-E_AMSR2Format.L3.SMC_10.8/...`
- `./download/GCOM-W.AMSR2_L3.SMC_10_3/...`

---

## 4. Merge to NetCDF

Merge AMSR-E + AMSR2 daily L3 HDF5 files (`*_01D_*.h5`) into NetCDF:

```bash
cd /data02/shiojiri/DATA/AMSR
python merge_amsr_l3_daily.py \
  --output-mode yearly \
  --output-dir /data02/shiojiri/DATA/AMSR/processed/l3_daily \
  --overwrite
```

### Output behavior

- Dimensions: `(time, lat, lon)`
- Duplicate observations at the same day/grid are **not averaged**
- Observations are separated by orbit direction and slots:
  - `soil_moisture_A_01`, `soil_moisture_A_02`, ...
  - `soil_moisture_D_01`, `soil_moisture_D_02`, ...
- Observation times are also stored:
  - `observation_time_min_A_XX`, `observation_time_min_D_XX`
- Day-shift correction:
  - `>= 1440 min` is moved to next day
  - `< 0 min` is moved to previous day

---

## 5. Query observations by datetime range

`amsr_l3_query.py` returns merged L3 observations in a datetime range as:

- `(soil_moisture, observation_datetime, lat, lon)`

Grid shape is flattened, and unobserved (`NaN`) cells are excluded.

### Python usage

```python
import datetime as dt
from amsr_l3_query import AmsrL3ObservationReader

reader = AmsrL3ObservationReader(
    "/data02/shiojiri/DATA/AMSR/processed/l3_daily",
    interval_hours=3,  # auto-decide cache size from fixed interval
)

rows = reader.read_range(
    "2012-07-03T00:00:00",
    "2012-07-03T05:59:59",
)

print(len(rows))
print(rows[0])  # (sm, datetime, lat, lon)
```

Use `read_range()` when you only need one range extraction:

```python
rows = reader.read_range(
    "2012-07-03T00:00:00",
    "2012-07-03T05:59:59",
)
```

Read only one time window with `iterate_windows()`:

```python
import datetime as dt
from amsr_l3_query import AmsrL3ObservationReader

reader = AmsrL3ObservationReader("/data02/shiojiri/DATA/AMSR/processed/l3_daily")

for ws, we, rows in reader.iterate_windows(
    dt.datetime(2012, 7, 3, 0, 0, 0),
    dt.datetime(2012, 7, 3, 23, 59, 59),
    step=dt.timedelta(days=1),
    window=dt.timedelta(days=1),
):
    print(ws, we, len(rows))
    # executes exactly one loop
```

Useful reader utilities:

- `preload_range(start_datetime, end_datetime)`: preloads day data into cache before loops.
- `clear_cache()`: manually clear in-memory cache.
- `AmsrL3ObservationReader.suggest_cache_days(interval_hours, window_hours=None)`: convenience for `max_cache_days`.

### Iterative usage (3-hour / 6-hour / daily)

```python
import datetime as dt
from amsr_l3_query import AmsrL3ObservationReader

reader = AmsrL3ObservationReader(
    "/data02/shiojiri/DATA/AMSR/processed/l3_daily",
    interval_hours=6,
)

for ws, we, rows in reader.iterate_windows(
    "2012-07-01T00:00:00",
    "2012-07-03T23:59:59",
    step=dt.timedelta(hours=6),
):
    print(ws, we, len(rows))
```

`iterate_windows()` is preferred when you process many repeated intervals (e.g. 3-hour / 6-hour / daily windows) with one reader.

### CLI usage

```bash
python amsr_l3_query.py \
  /data02/shiojiri/DATA/AMSR/processed/l3_daily \
  2012-07-03T00:00:00 2012-07-03T23:59:59 \
  --interval-hours 3 \
  --limit 10
```

### 0.5-degree query usage

```bash
python query_amsr_l3_0p5deg.py \
  processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc \
  2012-07-03T00:00:00 \
  2012-07-03T23:59:59 \
  --interval-hours 3 \
  --limit 10
```

The output tuple is the same as `amsr_l3_query.py`:

- `(soil_moisture, observation_datetime, lat, lon)`

For 0.5-degree averaged files, `observation_datetime` is the day timestamp at `00:00:00`.

### 0.5-degree upscaling

Run after merge:

```bash
cd /data02/shiojiri/DATA/AMSR
python upscale_amsr_l3_0p5deg.py \
  processed/l3_daily \
  processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc \
  --overwrite
```

Relative defaults are also supported:

```bash
cd /data02/shiojiri/DATA/AMSR
python upscale_amsr_l3_0p5deg.py
```

### 0.5-degree query usage with Python

```python
from query_amsr_l3_0p5deg import read_upsampled_amsr_observations_in_range

rows = read_upsampled_amsr_observations_in_range(
    path="processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc",
    start_datetime="2012-07-03T00:00:00",
    end_datetime="2012-07-03T23:59:59",
)

print(len(rows))
print(rows[0])
```

```python
import datetime as dt
from query_amsr_l3_0p5deg import AmsrUpsampledL3ObservationReader

reader = AmsrUpsampledL3ObservationReader(
    "processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc",
    interval_hours=24,
)

for ws, we, rows in reader.iterate_windows(
    dt.datetime(2012, 7, 1),
    dt.datetime(2012, 7, 3, 23, 59, 59),
    step=dt.timedelta(days=1),
):
    print(ws, we, len(rows))
```

Single-step example:

```python
import datetime as dt
from query_amsr_l3_0p5deg import AmsrUpsampledL3ObservationReader

reader = AmsrUpsampledL3ObservationReader(
    "processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc",
)

for ws, we, rows in reader.iterate_windows(
    dt.datetime(2012, 7, 3),
    dt.datetime(2012, 7, 3, 23, 59, 59),
    step=dt.timedelta(days=1),
    window=dt.timedelta(days=1),
):
    print(ws, we, len(rows))
    # executes exactly one loop
```

`AmsrUpsampledL3ObservationReader` automatically detects the input type. It reads cell-level daily means for averaged 0.5-degree files (`soil_moisture` + `observation_count`), or falls back to `amsr_l3_query.py` behavior for legacy slot-format files.

For 0.5-degree averaged files, use `read_range()` through `AmsrUpsampledL3ObservationReader.read_range(...)` (via the class or `read_upsampled_amsr_observations_in_range`). `preload_range()` is available when the input is legacy slot format; for averaged 0.5-degree files, the class uses per-day cached arrays and repeated `read_range()`/`iterate_windows()` access is already efficient.

---

## 6. Frequently used options

### `merge_amsr_l3_daily.py`

- `--start-date YYYY-MM-DD`
- `--end-date YYYY-MM-DD`
- `--output-mode single|yearly`
- `--output-file` (for `single` mode)
- `--max-backward-days` (default: `1`)

### `upscale_amsr_l3_0p5deg.py`

- `input_path`: merged NetCDF file or directory (e.g., `processed/l3_daily`)
- `output_path`: output NetCDF path
- `--overwrite`: overwrite existing output
- `--workers N`: process per-day averaging in parallel (default: `min(4, cpu_count)`)
- `--no-progress`: disable progress bar (`tqdm` if installed)
- Output variables:
  - `soil_moisture`: daily mean soil moisture in each 0.5-degree cell
  - `observation_count`: number of observations used in each cell average
- `time` dimension is daily (`time` values are days since 1970-01-01), with variables shaped `(time, lat, lon)`

### `query_amsr_l3_0p5deg.py`

- `path`: upsampled 0.5-degree NetCDF file or directory (default:
  `processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg_avg.nc`)
- `start`: start datetime (`YYYY-MM-DDTHH:MM:SS` or ISO format supported)
- `end`: end datetime (`YYYY-MM-DDTHH:MM:SS` or ISO format supported)
- `--cache-days`: explicit cache size in days
- `--interval-hours`: fixed query step (hours), auto cache sizing if `--cache-days` is omitted
- `--window-hours`: query window length (hours, optional)
- `--cache-days`, `--interval-hours`, `--window-hours` are accepted for compatibility. For averaged files, returned `observation_datetime` is the day timestamp at `00:00:00`.
- Helper function: `read_upsampled_amsr_observations_in_range(path, start_datetime, end_datetime, ...)` is available for one-shot read.

### `amsr_l3_query.py`

- `read_amsr_observations_in_range(path, start_datetime, end_datetime, ...)` is available for one-shot read.
- `--cache-days` (explicit cache size)
- `--interval-hours` (auto cache sizing from fixed step)
- `--window-hours` (if query window differs from step)

---

## 7. Notes

- `AMSR_download.sh` is configured by editing its product block.
- `merge_amsr_l3_daily.py` processes only `*_01D_*.h5` files.
- Time handling is UTC-based.
