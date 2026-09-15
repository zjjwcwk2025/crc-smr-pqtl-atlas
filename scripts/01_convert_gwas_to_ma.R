#!/usr/bin/env Rscript
# Phase 1 Step 1: Convert GCST90255675 GWAS to SMR COJO .ma format
# Format: SNP A1 A2 freq b se p N
# Project 9 v5: CRC Multi-Omics Drug-Target Atlas

library(data.table)

# --- Config ---
gwas_file <- "data/gwas/GCST90255675.h.tsv.gz"
out_file  <- "data/gwas/GCST90255675.ma"

# Sample size from Fernandez-Rozadilla et al. 2023 Nat Genet
# European subset: 78,473 cases + 107,143 controls
N <- 185616

cat("[1/3] Reading GWAS...\n")
gwas <- fread(gwas_file, sep = "\t", header = TRUE,
              select = c("rsid", "effect_allele", "other_allele",
                         "effect_allele_frequency", "beta", "standard_error", "p_value"))
cat(sprintf("  SNPs read: %s\n", format(nrow(gwas), big.mark=",")))

# --- Filter: remove NA essentials ---
cat("[2/3] Cleaning...\n")
before <- nrow(gwas)

# Remove rows with NA in critical columns
gwas <- gwas[!is.na(rsid) & rsid != ""]
gwas <- gwas[!is.na(beta) & !is.na(standard_error) & !is.na(p_value)]
gwas <- gwas[!is.na(effect_allele) & !is.na(other_allele)]

# Remove ambiguous/indel variants (keep only single-character alleles)
gwas <- gwas[nchar(effect_allele) == 1 & nchar(other_allele) == 1]

# Remove strand-ambiguous SNPs (A/T, G/C pairs)
is_palindromic <- (gwas$effect_allele == "A" & gwas$other_allele == "T") |
                  (gwas$effect_allele == "T" & gwas$other_allele == "A") |
                  (gwas$effect_allele == "G" & gwas$other_allele == "C") |
                  (gwas$effect_allele == "C" & gwas$other_allele == "G")
gwas <- gwas[!is_palindromic]

# Handle missing freq: set to NA (SMR accepts NA freq)
gwas[effect_allele_frequency == "NA", effect_allele_frequency := NA]

after <- nrow(gwas)
cat(sprintf("  After filtering: %s (%s removed)\n",
            format(after, big.mark=","), format(before - after, big.mark=",")))

# --- Build COJO .ma format ---
cat("[3/3] Writing .ma file...\n")
ma <- data.table(
    SNP  = gwas$rsid,
    A1   = gwas$effect_allele,
    A2   = gwas$other_allele,
    freq = gwas$effect_allele_frequency,
    b    = gwas$beta,
    se   = gwas$standard_error,
    p    = gwas$p_value,
    N    = N
)

fwrite(ma, out_file, sep = "\t", quote = FALSE, na = "NA",
       row.names = FALSE, col.names = TRUE)

cat(sprintf("\nDone! Output: %s\n", out_file))
cat(sprintf("  SNPs: %s\n", format(nrow(ma), big.mark=",")))
cat(sprintf("  File size: %s MB\n", round(file.size(out_file) / 1e6, 1)))
