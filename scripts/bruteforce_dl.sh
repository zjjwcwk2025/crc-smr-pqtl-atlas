#!/bin/bash
set -euo pipefail
PLINK2=/ifs1/User/zhouman/project9-v5-crc-atlas/tools/plink2
VDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/vcf
PDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/plink_eur
KEEP="$PDIR/eur_samples.txt"
LOG=/ifs1/User/zhouman/project9-v5-crc-atlas/logs/bruteforce.log
FTP="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

log() { echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

dl_one() {
    local chr=$1
    local vcf="$VDIR/chr${chr}.vcf.gz"
    local bed="$PDIR/chr${chr}_eur.bed"
    local url="${FTP}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

    if [ -f "$bed" ]; then
        log "chr${chr}: BED exists, SKIP"
        return 0
    fi

    for attempt in $(seq 1 20); do
        log "chr${chr}: attempt ${attempt}/20"
        rm -f "$vcf"

        if ! wget --no-proxy -t 3 --timeout=180 -O "$vcf" "$url" 2>&1 | tail -2 | tee -a "$LOG"; then
            log "chr${chr}: wget failed, retry in 30s"
            sleep 30
            continue
        fi
        SZ=$(ls -lh "$vcf" | awk '{print $5}')
        log "chr${chr}: downloaded $SZ"

        # 直接跑PLINK，不验证gzip
        log "chr${chr}: running PLINK2..."
        if $PLINK2 --vcf "$vcf" --keep "$KEEP" \
            --max-alleles 2 --snps-only just-acgt --make-bed \
            --out "$PDIR/chr${chr}_eur" --threads 16 --memory 32000 2>&1 | tail -3 | tee -a "$LOG"; then
            if [ -f "$bed" ]; then
                log "chr${chr}: PLINK DONE! $(ls -lh "$bed" | awk '{print $5}')"
                return 0
            fi
        fi
        log "chr${chr}: PLINK failed, retry in 30s"
        rm -f "$vcf"
        sleep 30
    done
    log "chr${chr}: FAILED after 20 attempts"
    return 1
}

log "===== BRUTEFORCE MODE ====="
for chr in 1 2 12; do
    bed="$PDIR/chr${chr}_eur.bed"
    if [ -f "$bed" ]; then
        log "chr${chr}: already done, skip"
    else
        log "chr${chr}: starting download loop..."
    fi
done

# 串行但最多20次重试
for chr in 1 2 12; do
    dl_one $chr || log "chr${chr}: GAVE UP"
done

log "===== FINAL ====="
for chr in 1 2 12; do
    if [ -f "$PDIR/chr${chr}_eur.bed" ]; then
        log "chr${chr}: ✅ $(ls -lh "$PDIR/chr${chr}_eur.bed" | awk '{print $5}')"
    else
        log "chr${chr}: ❌"
    fi
done
