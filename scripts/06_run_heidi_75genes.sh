#!/bin/bash
set -euo pipefail

SMR=/ifs1/User/zhouman/project9-v5-crc-atlas/smr/build/smr
PROJ=/ifs1/User/zhouman/project9-v5-crc-atlas
LOG="${PROJ}/logs/heidi_75genes.log"
OUT="${PROJ}/results/phase1_heidi_75genes"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "========================================="
log "SMR HEIDI: 75 Bonferroni-significant genes"
log "========================================="

log "Input files:"
log "  Bfile: ${PROJ}/data/1kg_phase3/plink_eur/merged_eur"
log "  GWAS:  ${PROJ}/data/gwas/GCST90255675.ma"
log "  BESD:  ${PROJ}/data/eqtlgen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd"
log "  Probes: ${PROJ}/results/phase1_75genes_probes.txt ($(wc -l < ${PROJ}/results/phase1_75genes_probes.txt) probes)"
log ""

$SMR \
  --bfile "${PROJ}/data/1kg_phase3/plink_eur/merged_eur" \
  --gwas-summary "${PROJ}/data/gwas/GCST90255675.ma" \
  --beqtl-summary "${PROJ}/data/eqtlgen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd" \
  --extract-probe "${PROJ}/results/phase1_75genes_probes.txt" \
  --heidi-mtd 1 \
  --smr --out "$OUT" \
  --thread-num 16 2>&1 | tee -a "$LOG"

log ""
log "========================================="
log "HEIDI analysis complete"
log "Output: ${OUT}.smr"
log "========================================="
