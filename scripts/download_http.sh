#!/bin/bash
set -euo pipefail
PLINK2=/ifs1/User/zhouman/project9-v5-crc-atlas/tools/plink2
VDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/vcf
PDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/plink_eur
KEEP=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/plink_eur/eur_samples.txt
LOG=/ifs1/User/zhouman/project9-v5-crc-atlas/logs/plink_progress.log
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
mkdir -p "$VDIR" "$PDIR" "$(dirname "$LOG")"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

log "=== HTTP mode download + PLINK for chr1,2,12 ==="

for CHR in 1 2 12; do
    log "--- chr${CHR} ---"
    REMOTE="${BASE}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    VCF="$VDIR/chr${CHR}.vcf.gz"
    BED="$PDIR/chr${CHR}_eur.bed"
    
    # Download
    log "curl downloading chr${CHR} via HTTP..."
    curl -L --retry 3 --retry-delay 30 --connect-timeout 60 \
         -o "$VCF" "$REMOTE" 2>&1 | tail -3 | tee -a "$LOG"
    
    SZ=$(ls -lh "$VCF" | awk '{print $5}')
    log "  Downloaded chr${CHR}: $SZ"
    
    # Quick test: can zcat read it?
    if zcat "$VCF" 2>/dev/null | head -1 | grep -q "##fileformat=VCF"; then
        log "  chr${CHR}: zcat OK, running PLINK2..."
    else
        log "  chr${CHR}: zcat FAILED, skip"
        continue
    fi
    
    # PLINK
    $PLINK2 --vcf "$VCF" --keep "$KEEP" \
        --max-alleles 2 --snps-only just-acgt --make-bed \
        --out "$PDIR/chr${CHR}_eur" --threads 16 --memory 32000 2>&1 | tail -3 | tee -a "$LOG"
    
    if [ -f "$BED" ]; then
        log "  chr${CHR} PLINK: DONE $(ls -lh "$BED" | awk '{print $5}')"
    else
        log "  chr${CHR} PLINK: FAILED"
    fi
done
log "=== ALL DONE ==="
