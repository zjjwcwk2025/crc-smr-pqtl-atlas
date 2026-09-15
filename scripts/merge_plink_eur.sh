#!/bin/bash
set -euo pipefail

PLINK2=/ifs1/User/zhouman/project9-v5-crc-atlas/tools/plink2
PDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/plink_eur
LOG=/ifs1/User/zhouman/project9-v5-crc-atlas/logs/merge_eur.log

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "========================================="
log "Merging 10 chr PLINK EUR files → merged_eur"
log "========================================="

# Create merge list
MERGELIST="${PDIR}/merge_list.txt"
rm -f "$MERGELIST"
for chr in 1 2 3 6 8 11 12 18 19 20; do
    echo "${PDIR}/chr${chr}_eur" >> "$MERGELIST"
done

log "Merge list:"
cat "$MERGELIST" | tee -a "$LOG"

OUT="${PDIR}/merged_eur"

if [ -f "${OUT}.bed" ]; then
    log "merged_eur already exists, removing old files..."
    rm -f "${OUT}.bed" "${OUT}.bim" "${OUT}.fam" "${OUT}.log"
fi

log "Starting PLINK2 pmerge..."
log "Total input variants: ~40M"
log "Total samples: 503 EUR"

$PLINK2 --pmerge-list "$MERGELIST" bfile \
    --set-all-var-ids @:#_\$r_\$a \
    --make-bed --out "$OUT" \
    --threads 16 --memory 64000 2>&1 | tee -a "$LOG"

if [ -f "${OUT}.bed" ]; then
    size=$(ls -lh "${OUT}.bed" | awk '{print $5}')
    variants=$(wc -l < "${OUT}.bim")
    log "merged_eur: DONE (${size}, ${variants} variants)"
else
    log "merged_eur: FAILED"
    exit 1
fi

log "========================================="
log "Merge complete"
log "========================================="
