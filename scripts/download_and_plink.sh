#!/bin/bash
# 1000G Phase 3 VCF 下载 + PLINK EUR 转换
# 用法: 在 screen/tmux 里跑，防止 SSH 断连杀进程
# screen -S plink_run && bash scripts/download_and_plink.sh

set -euo pipefail

PLINK2=/ifs1/User/zhouman/project9-v5-crc-atlas/tools/plink2
VDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/vcf
PDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/plink_eur
KEEP=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/plink_eur/eur_samples.txt
LOG=/ifs1/User/zhouman/project9-v5-crc-atlas/logs/plink_progress.log
FTP="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
mkdir -p "$PDIR" "$(dirname "$LOG")"

# 需要用 wget 重下的染色体（aria2c 多线程损坏的）
MISSING="1 2 12"
CHROMS="${MISSING}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "========================================="
log "START: 下载 chr${MISSING// /,} + PLINK EUR 转换"
log "========================================="

for CHR in $CHROMS; do
    log "--- chr${CHR} ---"

    # Step 1: 下载（wget 单线程 + 续传）
    VCF="$VDIR/chr${CHR}.vcf.gz"
    REMOTE="${FTP}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

    if [ -f "$VCF" ]; then
        log "  VCF exists ($(ls -lh "$VCF" | awk '{print $5}')), checking integrity..."
        if gzip -t "$VCF" 2>/dev/null; then
            log "  chr${CHR}.vcf.gz: ALREADY VALID, skip download"
        else
            log "  chr${CHR}.vcf.gz: CORRUPT, removing and re-downloading..."
            rm -f "$VCF"
        fi
    fi

    if [ ! -f "$VCF" ]; then
        TMP="${VCF}.tmp"
        rm -f "$TMP"
        log "  Downloading chr${CHR} from 1000G FTP..."
        log "  Output: $VCF"
        # wget 单线程，支持续传，重试3次
        wget -c -O "$TMP" "$REMOTE" 2>&1 | tail -3 | tee -a "$LOG"
        mv "$TMP" "$VCF"

        # 验证
        log "  Verifying chr${CHR}.vcf.gz..."
        if gzip -t "$VCF" 2>/dev/null; then
            log "  chr${CHR}.vcf.gz: VALID"
        else
            log "  chr${CHR}.vcf.gz: STILL CORRUPT after download!"
            continue
        fi
    fi

    # Step 2: PLINK EUR 转换
    BED="$PDIR/chr${CHR}_eur.bed"
    if [ -f "$BED" ]; then
        log "  chr${CHR}_eur.bed: ALREADY DONE, skip PLINK"
    else
        log "  Running PLINK2 for chr${CHR}..."
        $PLINK2 --vcf "$VCF" \
                --keep "$KEEP" \
                --max-alleles 2 \
                --snps-only just-acgt \
                --make-bed \
                --out "$PDIR/chr${CHR}_eur" \
                --threads 16 \
                --memory 32000 2>&1 | tail -5 | tee -a "$LOG"

        if [ -f "$BED" ]; then
            SIZE=$(ls -lh "$BED" | awk '{print $5}')
            log "  chr${CHR}: PLINK DONE ($SIZE)"
        else
            log "  chr${CHR}: PLINK FAILED"
        fi
    fi

    log "--- chr${CHR} complete ---"
done

log "========================================="
log "ALL DONE: ${MISSING// /, }"
log "Next: 合并 merged_eur → 跑 HEIDI"
log "========================================="
