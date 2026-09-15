#!/usr/bin/env Rscript
# Phase 3c: Bayesian colocalization (coloc) — pQTL level
# UKB-PPP 5 genes: CDKN1A, LTA, NCF2, TET2, TNF
# Input: UKB-PPP pQTL merged files, GWAS .ma, 1000G BIM
# Output: results/phase3c_pqtl_coloc.csv

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(dplyr)
})

setDTthreads(8)

# ===== Config =====
BIM_FILE   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/ref_eur/1000G_EUR.bim"
GWAS_FILE  <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.ma"
PQRTL_DIR  <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/ukbppp_merged"
OUT_CSV    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3c_pqtl_coloc.csv"

CIS_WINDOW <- 1e6  # ±1Mb

# ===== 1. Define 5 UKB-PPP genes (hg38 GENCODE v47 TSS, matches 07/14) =====
genes <- data.table(
  gene    = c("CDKN1A", "LTA",   "NCF2",  "TET2", "TNF"),
  ensg    = c("ENSG00000124762", "ENSG00000226979", "ENSG00000116701",
              "ENSG00000168769", "ENSG00000232810"),
  chr     = c(6, 6, 1, 4, 6),
  pos     = c(36652392, 31572054, 183555568, 105145428, 31575565),
  p_SMR   = c(1.51e-11, 9.12e-09, 4.09e-07, 1.46e-06, 1.31e-10),
  pqtl_file = c("CDKN1A_merged.txt.gz", "LTA_merged.txt.gz",
                "NCF2_merged.txt.gz", "TET2_merged.txt.gz", "TNF_merged.txt.gz")
)

cat(sprintf("pQTL coloc for %d UKB-PPP genes\n", nrow(genes)))
print(genes[, .(gene, chr, pos, p_SMR)])

# ===== 2. Load BIM for GWAS SNP → chr:pos mapping =====
cat("\nLoading BIM...\n")
bim <- fread(BIM_FILE, header = FALSE, select = 1:6,
             col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
bim[, chr := as.integer(chr)]
cat(sprintf("BIM: %d SNPs\n", nrow(bim)))

# ===== 3. Pre-load GWAS data (filter to relevant chromosomes) =====
cat("Loading GWAS data...\n")
target_chrs <- unique(genes$chr)
gwas_all <- fread(GWAS_FILE, header = TRUE,
                  select = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N"))
cat(sprintf("GWAS SNPs loaded: %d\n", nrow(gwas_all)))

# Map GWAS SNPs to chr:pos via BIM
gwas_all <- merge(gwas_all, bim[, .(rsid, chr, pos)],
                  by.x = "SNP", by.y = "rsid", all.x = FALSE)
gwas_all <- gwas_all[chr %in% target_chrs]
cat(sprintf("GWAS SNPs on target chromosomes: %d\n", nrow(gwas_all)))

# ===== 4. Run coloc per gene =====
results <- list()

for (i in 1:nrow(genes)) {
  gene     <- genes$gene[i]
  gchr     <- genes$chr[i]
  gpos     <- genes$pos[i]
  pqtl_f   <- file.path(PQRTL_DIR, genes$pqtl_file[i])

  cat(sprintf("\n=== %s (chr%d:%d) ===\n", gene, gchr, gpos))

  # ---- 4a. Load pQTL data for this chromosome only ----
  cmd <- sprintf("zcat %s | awk 'NR==1 || $1==%d'", pqtl_f, gchr)
  pqtl <- fread(cmd = cmd, header = TRUE)
  if (nrow(pqtl) < 10) {
    cat(sprintf("  <10 pQTL SNPs on chr%d\n", gchr)); next
  }

  # Parse chr:pos from ID column
  pqtl[, c("pqtl_chr", "pqtl_pos") := {
    parts <- tstrsplit(ID, ":")
    .(as.integer(parts[[1]]), as.integer(parts[[2]]))
  }]

  # Filter to cis window
  pqtl_cis <- pqtl[pqtl_pos >= gpos - CIS_WINDOW & pqtl_pos <= gpos + CIS_WINDOW]
  cat(sprintf("  pQTL SNPs in cis: %d\n", nrow(pqtl_cis)))
  if (nrow(pqtl_cis) < 10) {
    cat("  Too few cis pQTL SNPs\n"); next
  }

  # ---- 4b. GWAS SNPs in same region ----
  gwas_cis <- gwas_all[chr == gchr &
                       pos >= gpos - CIS_WINDOW &
                       pos <= gpos + CIS_WINDOW]
  cat(sprintf("  GWAS SNPs in cis: %d\n", nrow(gwas_cis)))
  if (nrow(gwas_cis) < 10) {
    cat("  Too few GWAS SNPs\n"); next
  }

  # ---- 4c. Merge by chr:pos ----
  # UKB-PPP: ALLELE0=ref, ALLELE1=effect, BETA is effect of ALLELE1
  # GWAS: A1=effect, A2=ref, b is effect of A1

  merged <- merge(
    pqtl_cis[, .(pqtl_pos, ALLELE0, ALLELE1, A1FREQ, BETA, SE, LOG10P, N)],
    gwas_cis[, .(pos, SNP, A1, A2, freq, b, se, p, N)],
    by.x = "pqtl_pos", by.y = "pos",
    suffixes = c("_pqtl", "_gwas")
  )

  cat(sprintf("  Overlapping SNPs (chr:pos): %d\n", nrow(merged)))
  if (nrow(merged) < 10) {
    cat("  Too few overlapping SNPs\n"); next
  }

  # ---- 4d. Deduplicate by position (keep best pQTL p-value) ----
  # Some positions have multiple GWAS SNPs (indels, multi-allelic) →
  # duplicate SNP names in coloc's snp vector
  merged[, dup_count := .N, by = pqtl_pos]
  if (any(merged$dup_count > 1)) {
    n_dup <- sum(merged$dup_count > 1)
    cat(sprintf("  Deduplicating: %d positions with >1 SNP\n", n_dup))
    # Keep the row with smallest pQTL p-value (most significant)
    merged <- merged[order(pqtl_pos, LOG10P)]
    merged <- merged[!duplicated(pqtl_pos)]
    merged$dup_count <- NULL
    cat(sprintf("  After dedup: %d unique positions\n", nrow(merged)))
  }
  merged$dup_count <- NULL

  # ---- 4e. Allele harmonization ----
  merged[, allele_match := (ALLELE1 == A1 & ALLELE0 == A2)]
  merged[, allele_flip  := (ALLELE1 == A2 & ALLELE0 == A1)]

  merged <- merged[allele_match | allele_flip]
  cat(sprintf("  After allele harmonization: %d SNPs\n", nrow(merged)))
  if (nrow(merged) < 10) next

  # ---- 4e. Compute p-values ----
  # pQTL p-value from LOG10P
  merged[, p_pqtl := 10^(-pmax(LOG10P, 0))]

  # Cap exact-zero p-values
  merged[, p_pqtl := pmax(p_pqtl, 1e-300)]
  merged[, p_gwas := pmax(p, 1e-300)]

  # MAF from pQTL A1FREQ
  merged[, maf := pmin(A1FREQ, 1 - A1FREQ, na.rm = TRUE)]
  merged[is.na(maf), maf := 0.1]

  # ---- 4f. Run coloc ----
  snp_ids <- merged$SNP
  N_gwas <- merged$N_gwas[1]

  # GWAS: case-control
  D1 <- list(
    snp = snp_ids,
    pvalues = merged$p_gwas,
    N = N_gwas,
    type = "cc",
    s = 0.42,  # CRC prevalence
    MAF = merged$maf
  )

  # pQTL: quantitative (protein levels)
  D2 <- list(
    snp = snp_ids,
    pvalues = merged$p_pqtl,
    N = median(merged$N_pqtl, na.rm = TRUE),
    type = "quant",
    MAF = merged$maf
  )

  tryCatch({
    col_res <- coloc.abf(dataset1 = D1, dataset2 = D2)
    summ <- col_res$summary
    pp_h4 <- as.numeric(summ[["PP.H4.abf"]])
    pp_h3 <- as.numeric(summ[["PP.H3.abf"]])

    # Get individual SNP posteriors
    snp_pp <- col_res$results
    top_snp_row <- snp_pp[which.max(snp_pp$SNP.PP.H4), ]
    top_snp <- top_snp_row$snp
    top_snp_pph4 <- top_snp_row$SNP.PP.H4

    flag <- ifelse(pp_h4 > 0.8, "✅ COLOC",
                   ifelse(pp_h4 > 0.5, "⚠ Marginal", "❌ No coloc"))

    cat(sprintf("  PPH4=%.4f, PPH3=%.4f, nsnp=%d, topSNP=%s(SNP.PP=%.3f) %s\n",
                pp_h4, pp_h3, length(snp_ids), top_snp, top_snp_pph4, flag))

    results[[length(results) + 1]] <- data.frame(
      gene    = gene,
      chr     = gchr,
      pos     = gpos,
      nsnps   = length(snp_ids),
      PPH4    = pp_h4,
      PPH3    = pp_h3,
      top_snp = top_snp,
      top_snp_pph4 = top_snp_pph4,
      p_SMR   = genes$p_SMR[i],
      verdict = flag,
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    cat(sprintf("  coloc error: %s\n", e$message))
  })
}

# ===== 5. Output =====
cat("\n\n=== pQTL Coloc Results ===\n")

if (length(results) > 0) {
  out <- rbindlist(results)
  out <- out[order(-PPH4)]
  fwrite(out, OUT_CSV)
  print(out)
  cat(sprintf("\nSaved: %s\n", OUT_CSV))
} else {
  cat("\nNo genes passed coloc.\n")

  # Create empty file so we know it ran
  fwrite(data.table(gene=character()), OUT_CSV)
}

cat("Done.\n")
