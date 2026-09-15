#!/bin/bash
# Phase 1b: Convert 1000G VCF → PLINK1 (EUR only) and run HEIDI test
# Project 9 v5: CRC Multi-Omics Drug-Target Atlas

set -e

PROJ=/ifs1/User/zhouman/project9-v5-crc-atlas
PLINK2=$PROJ/tools/plink2
SMR=$PROJ/smr/build/smr
VCF_DIR=$PROJ/data/1kg_phase3/vcf
PLINK_DIR=$PROJ/data/1kg_phase3/plink_eur
PANEL=$PROJ/data/1kg_phase3/integrated_call_samples_v3.20130502.ALL.panel

# Chromosomes with significant hits
CHROMS="1 2 3 6 8 11 12 18 19 20"

mkdir -p $PLINK_DIR

# --- Step 1: Create EUR sample list ---
echo "=== Step 1: Creating EUR sample list ==="

if [ ! -f "$PANEL" ]; then
    echo "ERROR: Panel file not found at $PANEL"
    exit 1
fi

# The panel file format: sample pop super_pop gender
# Extract EUR samples
awk '$3 == "EUR" {print $1}' $PANEL > $PLINK_DIR/eur_samples.txt
EUR_COUNT=$(wc -l < $PLINK_DIR/eur_samples.txt)
echo "EUR samples: $EUR_COUNT"

# --- Step 2: Convert each chromosome VCF → PLINK1 (EUR only) ---
echo ""
echo "=== Step 2: VCF → PLINK1 conversion ==="

for CHR in $CHROMS; do
    VCF=$VCF_DIR/chr${CHR}.vcf.gz
    OUT=$PLINK_DIR/chr${CHR}_eur

    if [ -f "${OUT}.bed" ] && [ -f "${OUT}.bim" ] && [ -f "${OUT}.fam" ]; then
        echo "chr${CHR}: Already converted, skipping"
        continue
    fi

    if [ ! -f "$VCF" ]; then
        echo "chr${CHR}: VCF not found, skipping"
        continue
    fi

    echo "chr${CHR}: Converting VCF → PLINK1 (EUR only)..."
    $PLINK2 --vcf $VCF \
            --keep $PLINK_DIR/eur_samples.txt \
            --max-alleles 2 \
            --snps-only just-acgt \
            --make-bed \
            --out $OUT \
            --threads 8 \
            --memory 16000 2>&1 | tail -5

    echo "chr${CHR}: $(wc -l < ${OUT}.bim) SNPs"
done

echo ""
echo "=== All conversions complete ==="
ls -lh $PLINK_DIR/*.bed $PLINK_DIR/*.bim $PLINK_DIR/*.fam
