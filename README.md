# AMSR Download / Merge / Query

This repository provides a workflow to:

1. Download AMSR products from JAXA G-Portal with `AMSR_download.sh`
2. Merge daily L3 files into NetCDF with `merge_amsr_l3_daily.py`
3. Query observations by datetime range with `amsr_l3_query.py`

---

## 1. Requirements

- `bash`
- `lftp` (with SFTP support)
- Python 3.10+
- Python packages:
  - `numpy`
  - `netCDF4`

Install Python packages:

```bash
python -m pip install numpy netCDF4
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

`amsr_l3_query.py` returns all observations in a datetime range as:

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

### CLI usage

```bash
python amsr_l3_query.py \
  /data02/shiojiri/DATA/AMSR/processed/l3_daily \
  2012-07-03T00:00:00 2012-07-03T23:59:59 \
  --interval-hours 3 \
  --limit 10
```

### 0.5-degree upscaling

Run after merge:

```bash
cd /data02/shiojiri/DATA/AMSR
python upsample_amsr_l3_0p5deg.py \
  processed/l3_daily \
  processed/l3_daily_0p5/AMSR_SMC_daily_0p5deg.nc \
  --overwrite
```

Relative defaults are also supported:

```bash
cd /data02/shiojiri/DATA/AMSR
python upsample_amsr_l3_0p5deg.py
```

---

## 6. Frequently used options

### `merge_amsr_l3_daily.py`

- `--start-date YYYY-MM-DD`
- `--end-date YYYY-MM-DD`
- `--output-mode single|yearly`
- `--output-file` (for `single` mode)
- `--max-backward-days` (default: `1`)

### `upsample_amsr_l3_0p5deg.py`

- `input_path`: merged NetCDF file or directory (e.g., `processed/l3_daily`)
- `output_path`: output NetCDF path
- `--overwrite`: overwrite existing output

### `amsr_l3_query.py`

- `--cache-days` (explicit cache size)
- `--interval-hours` (auto cache sizing from fixed step)
- `--window-hours` (if query window differs from step)

---

## 7. Notes

- `AMSR_download.sh` is configured by editing its product block.
- `merge_amsr_l3_daily.py` processes only `*_01D_*.h5` files.
- Time handling is UTC-based.
