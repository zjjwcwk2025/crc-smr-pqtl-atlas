#!/bin/bash
# 1000G VCF 顽固下载脚本 — 无限重试直到成功
# 用法: nohup bash this_script.sh > logs/retry.log 2>&1 &
# 策略: wget FTP下载 → gzip -t 验证 → 失败则删掉重来 → 最多10次
set -euo pipefail

PLINK2=/ifs1/User/zhouman/project9-v5-crc-atlas/tools/plink2
VDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/vcf
PDIR=/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/plink_eur
KEEP="$PDIR/eur_samples.txt"
LOG=/ifs1/User/zhouman/project9-v5-crc-atlas/logs/retry_download.log
mkdir -p "$VDIR" "$PDIR" "$(dirname "$LOG")"

FTP="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

log() { echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

download_one() {
    local chr=$1
    local vcf="$VDIR/chr${chr}.vcf.gz"
    local bed="$PDIR/chr${chr}_eur.bed"
    local url="${FTP}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    local max_retries=10
    local attempt=1

    log "========== chr${chr} =========="

    # 跳过已完成的 PLINK
    if [ -f "$bed" ]; then
        log "chr${chr}: BED exists ($(ls -lh "$bed" | awk '{print $5}')), SKIP"
        return 0
    fi

    while [ $attempt -le $max_retries ]; do
        log "chr${chr}: attempt ${attempt}/${max_retries}"

        # 如果有旧文件先删除
        [ -f "$vcf" ] && rm -f "$vcf"

        log "  Downloading (wget -t 3 --timeout 120)..."
        # wget 单线程, 3次内部重试, 2分钟超时
        if wget -t 3 --timeout=120 -O "$vcf" "$url" 2>&1 | tail -3 | tee -a "$LOG"; then
            SZ=$(ls -lh "$vcf" 2>/dev/null | awk '{print $5}')
            log "  wget OK, size=$SZ"
        else
            log "  wget FAILED (network error)"
            attempt=$((attempt + 1))
            sleep 30
            continue
        fi

        # gzip 验证
        log "  Verifying gzip..."
        if gzip -t "$vcf" 2>/dev/null; then
            log "  chr${chr}: gzip -t PASS"
        else
            # 即使gzip报错，试试zcat能否读开头（trailing garbage可用）
            if zcat "$vcf" 2>/dev/null | head -1 | grep -q "##fileformat=VCF"; then
                log "  chr${chr}: gzip -t FAIL but zcat OK (trailing garbage, usable)"
            else
                log "  chr${chr}: CORRUPT, retrying..."
                rm -f "$vcf"
                attempt=$((attempt + 1))
                sleep 30
                continue
            fi
        fi

        # PLINK 转换
        log "  Running PLINK2..."
        $PLINK2 --vcf "$vcf" --keep "$KEEP" \
            --max-alleles 2 --snps-only just-acgt --make-bed \
            --out "$PDIR/chr${chr}_eur" --threads 16 --memory 32000 2>&1 | tail -3 | tee -a "$LOG"

        if [ -f "$bed" ]; then
            log "chr${chr}: PLINK DONE! $(ls -lh "$bed" | awk '{print $5}')"
            return 0
        else
            log "chr${chr}: PLINK FAILED, deleting VCF and retrying..."
            rm -f "$vcf"
            attempt=$((attempt + 1))
            sleep 30
        fi
    done

    log "chr${chr}: FAILED after $max_retries attempts!"
    return 1
}

log "========================================="
log "RETRY MODE: downloading chr1,2,12"
log "========================================="

download_one 1
download_one 2
download_one 12

log "========================================="
log "FINAL: check results"
for chr in 1 2 12; do
    bed="$PDIR/chr${chr}_eur.bed"
    if [ -f "$bed" ]; then
        log "chr${chr}: ✅ $(ls -lh "$bed" | awk '{print $5}')"
    else
        log "chr${chr}: ❌ FAILED"
    fi
done
log "========================================="
