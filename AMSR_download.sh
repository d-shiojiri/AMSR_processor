#!/usr/bin/env bash
# Requirements: lftp (with sftp support)
# Credentials: read from ~/.gportal_sftp.env (chmod 600)

HOST="ftp.gportal.jaxa.jp"
PORT=2051

##############################
# AMSR2 L2
##############################
# START_YEAR=2012
# END_YEAR=2024
# REMOTE_BASE="standard/GCOM-W/GCOM-W.AMSR2/L2.SMC/3"
# LOCAL_BASE="./download/GCOM-W.AMSR2.L2.SMC.3"

##############################
# AMSR-E L2
##############################
# START_YEAR=2002
# END_YEAR=2011
# REMOTE_BASE="standard/AQUA/AQUA.AMSR-E_AMSR2Format/L2.SMC/8/"
# LOCAL_BASE="./download/AQUA.AMSR-E_AMSR2Format.L2.SMC.8"

##############################
# AMSR-E L3 10 km
##############################
# START_YEAR=2002
# END_YEAR=2011
# REMOTE_BASE="standard/AQUA/AQUA.AMSR-E_AMSR2Format/L3.SMC_10/8/"
# LOCAL_BASE="./download/AQUA.AMSR-E_AMSR2Format.L3.SMC_10.8"

##############################
# AMSR2 L3
##############################
START_YEAR=2012
END_YEAR=2025
REMOTE_BASE="/standard/GCOM-W/GCOM-W.AMSR2/L3.SMC_10/3"
LOCAL_BASE="./download/GCOM-W.AMSR2_L3.SMC_10_3"

MISSING_LOG="./missing_dirs.log"
EXIT_MISSING=44

set -euo pipefail

# Load credentials without exposing them to history
# ~/.gportal_sftp.env should contain:
# GPORTAL_USER="your_username"
# GPORTAL_PASS="your_password"
source ~/.gportal_sftp.env

: > "$MISSING_LOG"

for ((y=START_YEAR; y<=END_YEAR; y++)); do
  for m in {01..12}; do
    remote="${REMOTE_BASE}/${y}/${m}"
    local="${LOCAL_BASE}/${y}/${m}"
    mkdir -p "$local"

    if lftp -e "
      set sftp:auto-confirm yes
      set net:timeout 30
      set net:max-retries 100
      set net:persist-retries 100
      set mirror:use-pget-n 4

      open -u ${GPORTAL_USER},${GPORTAL_PASS} -p ${PORT} sftp://${HOST}
      cd ${remote} || exit ${EXIT_MISSING}
      mirror -c --only-newer --parallel=4 --verbose . ${local}
      bye
    "; then
      :
    else
      rc=$?
      if [ "$rc" -eq "$EXIT_MISSING" ]; then
        echo "${y}/${m}" >> "$MISSING_LOG"
        continue
      else
        echo "Error on ${y}/${m} (exit ${rc})" >&2
        exit "$rc"
      fi
    fi
  done
done

echo "Done. Missing list -> ${MISSING_LOG}"
