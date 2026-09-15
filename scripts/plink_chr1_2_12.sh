#!/bin/bash
set -euo pipefail

PLINK2=/ifs1/User/zhouman/project9-v5-crc-atlas/tools/plink2
VDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/vcf
PDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/plink_eur
KEEP=$PDIR/eur_samples.txt
LOG=/ifs1/User/zhouman/project9-v5-crc-atlas/logs/plink_chr1_2_12.log

mkdir -p "$PDIR" "$(dirname "$LOG")"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "========================================="
log "PLINK EUR conversion: chr1, chr2, chr12"
log "========================================="

for chr in 1 2 12; do
    VCF="$VDIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    OUT="$PDIR/chr${chr}_eur"
    
    if [ -f "${OUT}.bed" ] && [ -f "${OUT}.bim" ] && [ -f "${OUT}.fam" ]; then
        log "chr${chr}: PLINK files already exist, skip."
        continue
    fi
    
    log "--- chr${chr} ---"
    log "  VCF: $VCF"
    log "  OUT: $OUT"
    
    $PLINK2 --vcf "$VCF" --keep "$KEEP" \
        --max-alleles 2 --snps-only just-acgt --make-bed \
        --out "$OUT" --threads 16 --memory 32000 2>&1 | tee -a "$LOG"
    
    if [ -f "${OUT}.bed" ]; then
        log "  chr${chr}: DONE ($(ls -lh ${OUT}.bed | awk '{print $5}'))"
    else
        log "  chr${chr}: FAILED"
    fi
done

log "========================================="
log "ALL DONE"
log "========================================="
