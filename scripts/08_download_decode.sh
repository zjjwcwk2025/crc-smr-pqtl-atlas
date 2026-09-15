#!/bin/bash
# deCODE pQTL 批量下载 (17 files, 14 genes)
# 策略: 全并行 — 每文件独立 TCP 连接，避免 token 过期时某文件进度为 0
# 特性: curl单线程 + nohup防断连 + gzip -t验证 + 日志
set -euo pipefail

TOKEN="789e8bb8-6b9d-4493-a395-e8f13c5450ac"
DIR="/ifs1/User/zhouman/project9-v5-crc-atlas/data/decode"
LOGDIR="${DIR}/logs"
mkdir -p "$DIR" "$LOGDIR"

FILES=(
  "10047_12_NCF2_NCF_2.txt.gz"
  "10372_18_STAT6_STAT6.txt.gz"
  "11543_84_LIMA1_LIMA1.txt.gz"
  "12347_29_CCM2_CCM2.txt.gz"
  "13494_6_CERS5_CERS5.txt.gz"
  "14599_18_STAB1_STAB1.txt.gz"
  "14662_6_MZF1_MZF1.txt.gz"
  "15666_21_BMP2_BMP_2.txt.gz"
  "18291_8_CDKN1A_p21.txt.gz"
  "18419_20_ARPC5_p16_ARC.txt.gz"
  "19111_10_UBE2M_UBC12.txt.gz"
  "7105_7_CABLES2_CABL2.txt.gz"
  "5936_53_TNF_TNF_a.txt.gz"
  "4703_87_LTA_TNF_b.txt.gz"
  "3152_57_TNFRSF1B_TNF_sR_II.txt.gz"
  "2654_19_TNFRSF1A_TNF_sR_I.txt.gz"
  "8368_102_TNFRSF1B_TNF_sR_II.txt.gz"
)

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOGDIR}/master.log"; }

TOTAL=${#FILES[@]}
log "=== Starting deCODE download: ${TOTAL} files (ALL parallel) ==="

download_one() {
  local fname="$1"
  local OUTPUT="${DIR}/${fname}"
  local LOG="${LOGDIR}/${fname}.log"
  local URL="https://download.decode.is/s3/download?token=${TOKEN}&file=${fname}"

  echo "[$(date)] Starting: ${fname} (target: 955MB)" > "$LOG"

  if [ -f "$OUTPUT" ]; then
    if gzip -t "$OUTPUT" 2>/dev/null; then
      echo "[$(date)] ${fname}: SKIP (already valid, $(ls -lh "$OUTPUT" | awk '{print $5}'))" | tee -a "${LOGDIR}/master.log"
      return 0
    else
      local sz=$(stat -c%s "$OUTPUT" 2>/dev/null || echo 0)
      echo "[$(date)] ${fname}: Partial ${sz} bytes, will overwrite" | tee -a "${LOGDIR}/master.log"
      # Do NOT delete — let curl overwrite in place
    fi
  fi

  # No Range support on deCODE → can't resume. Just download full file.
  # --retry 5: retry on transient errors (connection drops)
  # --retry-delay 60: wait 60s between retries
  curl -sL --connect-timeout 300 \
       --retry 5 --retry-delay 60 \
       -o "$OUTPUT" \
       "$URL" >> "$LOG" 2>&1

  local rc=$?
  if [ $rc -eq 0 ] && gzip -t "$OUTPUT" 2>/dev/null; then
    log "✅ ${fname} DONE ($(ls -lh "$OUTPUT" | awk '{print $5}'))"
    return 0
  else
    local sz=$(stat -c%s "$OUTPUT" 2>/dev/null || echo 0)
    log "❌ ${fname} FAILED (curl rc=$rc, file size=${sz})"
    return 1
  fi
}

# Launch ALL 17 in parallel
PIDS=()
for f in "${FILES[@]}"; do
  download_one "$f" &
  PIDS+=($!)
done

log "All ${TOTAL} downloads launched (PIDs: ${PIDS[*]})"

# Wait for all
SUCCESS=0
FAILED=0
for pid in "${PIDS[@]}"; do
  wait $pid && SUCCESS=$((SUCCESS+1)) || FAILED=$((FAILED+1))
done

log "=== Complete: ${SUCCESS}/${TOTAL} OK, ${FAILED} FAILED ==="
