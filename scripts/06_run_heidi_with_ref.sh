#!/bin/bash
set -euo pipefail

SMR=/ifs1/User/zhouman/project9-v5-crc-atlas/smr/build/smr
BFILE=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/ref_eur/1000G_EUR
GWAS=/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.ma
BESD=/ifs1/User/zhouman/project9-v5-crc-atlas/data/eqtlgen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense
PROBES=/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase1_75genes_probes.txt
OUT=/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase1_heidi_75genes
LOG=/ifs1/User/zhouman/project9-v5-crc-atlas/logs/heidi_ref_new.log

mkdir -p "$(dirname "$OUT")" "$(dirname "$LOG")"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "========================================="
log "SMR HEIDI: 75 significant genes"
log "Reference: 1000G EUR with rsIDs (Zenodo)"
log "========================================="
log "BFILE: $BFILE"
log "GWAS: $GWAS"
log "BESD: $BESD"
log "Probes: $(wc -l < "$PROBES") probes"
log "========================================="

time $SMR \
  --bfile "$BFILE" \
  --gwas-summary "$GWAS" \
  --beqtl-summary "$BESD" \
  --extract-probe "$PROBES" \
  --heidi-mtd 1 \
  --diff-freq 0.2 \
  --diff-freq-prop 0.05 \
  --out "$OUT" \
  --thread-num 16 2>&1 | tee -a "$LOG"

if [ -f "${OUT}.heidi" ]; then
    log "HEIDI completed! Rows: $(tail -n +2 ${OUT}.heidi | wc -l)"
else
    log "HEIDI output not found. Check log above for errors."
fi

log "DONE"
