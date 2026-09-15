#!/bin/bash
set -euo pipefail

DECODE_DIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/decode
LOG=/ifs1/User/zhouman/project9-v5-crc-atlas/logs/decode_retry.log
TOKEN="30af6d7c-ebb3-4985-bf77-41ab862dd805"

mkdir -p "$(dirname "$LOG")"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "========================================="
log "Re-downloading 5 broken deCODE pQTL files"
log "========================================="

# Files to download: name, url
declare -A FILES
FILES["10047_12_NCF2_NCF_2.txt.gz"]="https://download.decode.is/s3/download?token=${TOKEN}&file=10047_12_NCF2_NCF_2.txt.gz"
FILES["13494_6_CERS5_CERS5.txt.gz"]="https://download.decode.is/s3/download?token=${TOKEN}&file=13494_6_CERS5_CERS5.txt.gz"
FILES["18291_8_CDKN1A_p21.txt.gz"]="https://download.decode.is/s3/download?token=${TOKEN}&file=18291_8_CDKN1A_p21.txt.gz"
FILES["3152_57_TNFRSF1B_TNF_sR_II.txt.gz"]="https://download.decode.is/s3/download?token=${TOKEN}&file=3152_57_TNFRSF1B_TNF_sR_II.txt.gz"
FILES["8368_102_TNFRSF1B_TNF_sR_II.txt.gz"]="https://download.decode.is/s3/download?token=${TOKEN}&file=8368_102_TNFRSF1B_TNF_sR_II.txt.gz"

for fname in "10047_12_NCF2_NCF_2.txt.gz" "13494_6_CERS5_CERS5.txt.gz" "18291_8_CDKN1A_p21.txt.gz" "3152_57_TNFRSF1B_TNF_sR_II.txt.gz" "8368_102_TNFRSF1B_TNF_sR_II.txt.gz"; do
    url="${FILES[$fname]}"
    out="$DECODE_DIR/$fname"
    
    log "--- $fname ---"
    
    # Remove existing broken file
    if [ -f "$out" ]; then
        log "  Removing existing broken file..."
        rm -f "$out"
    fi
    
    log "  Downloading..."
    log "  URL token: ${TOKEN:0:10}..."
    
    # curl with resume support, single-threaded (S3 signed URL limitation)
    curl -C - -sL --connect-timeout 60 --max-time 7200 \
        -o "$out" "$url" 2>&1 | tail -5 | tee -a "$LOG"
    
    if [ -f "$out" ]; then
        size=$(ls -lh "$out" | awk '{print $5}')
        log "  Downloaded: $size"
        
        # Verify gzip
        if gzip -t "$out" 2>/dev/null; then
            log "  $fname: VALID ✅"
        else
            log "  $fname: CORRUPT ❌ (may need new token)"
        fi
    else
        log "  $fname: DOWNLOAD FAILED ❌"
    fi
done

log "========================================="
log "deCODE re-download complete"
log "========================================="
