# AMSR Download / Merge / Query

このリポジトリは、以下の流れで AMSR データを扱うためのものです。

1. `AMSR_download.sh` で JAXA G-Portal からダウンロード  
2. `merge_amsr_l3_daily.py` で日次 NetCDF にマージ  
3. `amsr_l3_query.py` で指定 datetime 範囲の観測を抽出

---

## 1. 前提

- `bash`
- `lftp`（SFTP 対応）
- Python 3.10+
- Python パッケージ:
  - `numpy`
  - `netCDF4`

例:

```bash
python -m pip install numpy netCDF4
```

---

## 2. 認証ファイル（必須）

`AMSR_download.sh` は `~/.gportal_sftp.env` を読み込みます。  
まず次のファイルを作成してください。

```bash
cat > ~/.gportal_sftp.env << 'EOF'
GPORTAL_USER="your_gportal_username"
GPORTAL_PASS="your_gportal_password"
EOF
chmod 600 ~/.gportal_sftp.env
```

記法は **必ず** `KEY="value"` 形式にしてください。

---

## 3. ダウンロード

```bash
cd /data02/shiojiri/DATA/AMSR
bash AMSR_download.sh
```

### 補足

- ダウンロード対象（AMSR-E/AMSR2、L2/L3）は `AMSR_download.sh` 冒頭の設定ブロックで切り替えます。
- 現在のデフォルトは AMSR2 L3 (`GCOM-W.AMSR2_L3.SMC_10_3`) です。
- 欠損ディレクトリは `missing_dirs.log` に記録されます。

出力先の例:

- `./download/AQUA.AMSR-E_AMSR2Format.L3.SMC_10.8/...`
- `./download/GCOM-W.AMSR2_L3.SMC_10_3/...`

---

## 4. NetCDF へマージ

AMSR-E + AMSR2 の L3 日次 HDF5（`*_01D_*.h5`）を統合します。

```bash
cd /data02/shiojiri/DATA/AMSR
python merge_amsr_l3_daily.py \
  --output-mode yearly \
  --output-dir /data02/shiojiri/DATA/AMSR/processed/l3_daily \
  --overwrite
```

### 出力仕様（重要）

- 次元: `(time, lat, lon)`
- 同一日・同一格子で重複観測がある場合:
  - 平均化せず、`A/D` ごとに slot 変数へ保存
  - 例: `soil_moisture_A_01`, `soil_moisture_A_02`, `soil_moisture_D_01`, ...
- 観測時刻も保存:
  - `observation_time_min_A_XX`, `observation_time_min_D_XX`
- 時刻補正:
  - `>= 1440 min` は翌日へ
  - `< 0 min` は前日へ

---

## 5. 観測値の読み出し（datetime 範囲指定）

`amsr_l3_query.py` は、指定範囲の全観測を次のタプルで返します:

- `(soil_moisture, observation_datetime, lat, lon)`

格子形状は保持せず、未観測 (`NaN`) は返しません。

### Python から利用

```python
import datetime as dt
from amsr_l3_query import AmsrL3ObservationReader

reader = AmsrL3ObservationReader(
    "/data02/shiojiri/DATA/AMSR/processed/l3_daily",
    interval_hours=3,  # 固定インターバルから cache_days 自動決定
)

rows = reader.read_range(
    "2012-07-03T00:00:00",
    "2012-07-03T05:59:59",
)

print(len(rows))
print(rows[0])  # (sm, datetime, lat, lon)
```

### 反復利用（3時間 / 6時間 / 日単位）

```python
import datetime as dt
from amsr_l3_query import AmsrL3ObservationReader

reader = AmsrL3ObservationReader(
    "/data02/shiojiri/DATA/AMSR/processed/l3_daily",
    interval_hours=6,  # 6時間ステップ想定
)

for ws, we, rows in reader.iterate_windows(
    "2012-07-01T00:00:00",
    "2012-07-03T23:59:59",
    step=dt.timedelta(hours=6),
):
    print(ws, we, len(rows))
```

### CLI から確認

```bash
python amsr_l3_query.py \
  /data02/shiojiri/DATA/AMSR/processed/l3_daily \
  2012-07-03T00:00:00 2012-07-03T23:59:59 \
  --interval-hours 3 \
  --limit 10
```

---

## 6. よく使うオプション

### `merge_amsr_l3_daily.py`

- `--start-date YYYY-MM-DD`
- `--end-date YYYY-MM-DD`
- `--output-mode single|yearly`
- `--output-file`（`single` モード時）
- `--max-backward-days`（既定: `1`）

### `amsr_l3_query.py`

- `--cache-days`（明示指定）
- `--interval-hours`（自動 cache 設定用）
- `--window-hours`（窓長を step と分離する場合）

---

## 7. 注意点

- `AMSR_download.sh` は設定ブロックを書き換えて対象プロダクトを切替える設計です。
- `merge_amsr_l3_daily.py` は `*_01D_*.h5` を対象にしています。
- 時刻は UTC 前提で扱っています。
