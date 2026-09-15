#!/usr/bin/env Rscript
#'
#' Phase 5c-5: MR Sensitivity Analysis for SMR-significant genes
#' Since full TwoSampleMR requires LD clumping references (unavailable),
#' we compute sensitivity directly from the SMR + HEIDI results:
#' - Wald ratio (single-instrument MR) for each gene
#' - Steiger directionality test proxy (check if eQTL -> GWAS signal consistent)
#' - Heterogeneity assessment via Cochran's Q from HEIDI test
#'
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

PROJ <- "/ifs1/User/zhouman/project9-v5-crc-atlas"
setwd(PROJ)
dir.create("results/phase5c_sensitivity", showWarnings = FALSE, recursive = TRUE)

# ── 1. Load SMR + HEIDI results ──
cat("[1/3] Loading SMR and HEIDI results...\n")

smr <- fread("results/phase1_bonferroni_significant_annotated.tsv")
heidi <- fread("results/phase1_heidi_75genes.smr")
setnames(heidi, c("probeID", "ProbeChr", "Gene", "Probe_bp", "topSNP",
                  "topSNP_chr", "topSNP_bp", "A1", "A2", "Freq",
                  "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL", "se_eQTL", "p_eQTL",
                  "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"))

# ── 2. Compute sensitivity metrics ──
cat(sprintf("[2/3] Computing sensitivity for %d Bonferroni-significant genes...\n", nrow(smr)))

results <- data.frame(
  gene = character(),
  smr_b = numeric(), smr_se = numeric(), smr_p = numeric(),
  heidi_p = numeric(), heidi_pass = character(),
  n_snps = integer(),
  # Steiger directionality: r²_exposure vs r²_outcome
  eQTL_b = numeric(), eQTL_r2 = numeric(),   # eQTL effect
  GWAS_b = numeric(), GWAS_r2 = numeric(),    # GWAS effect
  steiger_pass = character(),
  sensitivity_tier = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(smr))) {
  gene_i <- smr$SYMBOL[i]
  probe_i <- smr$Gene[i]

  # Match HEIDI result
  heidi_row <- heidi[Gene == probe_i]
  heidi_p_i <- if (nrow(heidi_row) > 0) heidi_row$p_HEIDI[1] else NA
  heidi_nsnp_i <- if (nrow(heidi_row) > 0) heidi_row$nsnp_HEIDI[1] else NA
  heidi_pass_i <- if (is.na(heidi_p_i)) "NA" else if (heidi_p_i >= 0.01) "PASS" else "FAIL"

  smr_b <- smr$b_SMR[i]
  smr_se <- smr$se_SMR[i]
  smr_p <- smr$p_SMR[i]

  # eQTL effect stats
  eQTL_b <- smr$b_eQTL[i]
  eQTL_se <- smr$se_eQTL[i]
  eQTL_p <- smr$p_eQTL[i]

  # For Steiger: compare eQTL effect size (on gene expression) vs GWAS effect (on CRC)
  # eQTL r² = b² / (b² + se² * N) — approximate with beta estimates
  # A proper Steiger test needs N (sample sizes), which we approximate:
  # eQTLGen has ~31,684 samples; GWAS has ~185,616 samples
  eQTL_N <- 31684
  GWAS_N <- 185616

  eQTL_r2 <- (eQTL_b^2) / (eQTL_b^2 + eQTL_se^2 * eQTL_N)
  GWAS_b <- smr$b_GWAS[i]
  GWAS_se <- smr$se_GWAS[i]
  GWAS_r2 <- (GWAS_b^2) / (GWAS_b^2 + GWAS_se^2 * GWAS_N)

  # Steiger: exposure (eQTL) should explain more variance than outcome (GWAS)
  # If eQTL_r2 > GWAS_r2, the direction is correct (gene -> disease)
  steiger_pass_i <- if (eQTL_r2 > GWAS_r2) "PASS" else if (is.na(GWAS_b) || is.na(eQTL_b)) "NA" else "WARN"

  # Sensitivity tier:
  # A: HEIDI PASS  + direction correct
  # B: HEIDI PASS  + direction unclear
  # C: HEIDI FAIL  + direction correct
  # D: HEIDI FAIL  + direction unclear
  tier_i <- if (heidi_pass_i == "PASS" && steiger_pass_i == "PASS") "A - Robust" else
    if (heidi_pass_i == "PASS") "B - HEIDI only" else
    if (steiger_pass_i == "PASS") "C - Direction only" else "D - Weak"

  results <- rbind(results, data.frame(
    gene = gene_i,
    smr_b = smr_b, smr_se = smr_se, smr_p = smr_p,
    heidi_p = heidi_p_i, heidi_pass = heidi_pass_i,
    n_snps = heidi_nsnp_i,
    eQTL_b = eQTL_b, eQTL_r2 = eQTL_r2,
    GWAS_b = GWAS_b, GWAS_r2 = GWAS_r2,
    steiger_pass = steiger_pass_i,
    sensitivity_tier = tier_i,
    stringsAsFactors = FALSE
  ))
}

# ── 3. Summarize ──
cat("\n[3/3] Saving results...\n")

# Add Tier 1 status
susie <- fread("results/phase3_susie_finemap.csv")
tier1_list <- susie[tier == 1, gene]
results$is_tier1 <- results$gene %in% tier1_list

fwrite(results, "results/phase5c_sensitivity/mr_sensitivity_summary.csv")

# Summary by tier
cat("\n── MR Sensitivity Summary ──\n")
cat(sprintf("Total genes analyzed: %d\n", nrow(results)))

tier_tbl <- table(results$sensitivity_tier)
for (nm in names(tier_tbl)) {
  cat(sprintf("  %s: %d genes\n", nm, tier_tbl[nm]))
}

# HEIDI
heidi_tbl <- table(results$heidi_pass)
cat(sprintf("\nHEIDI: PASS=%d  FAIL=%d  NA=%d\n",
            sum(results$heidi_pass == "PASS", na.rm = TRUE),
            sum(results$heidi_pass == "FAIL", na.rm = TRUE),
            sum(results$heidi_pass == "NA", na.rm = TRUE)))

# Steiger
steiger_tbl <- table(results$steiger_pass)
cat(sprintf("Steiger directionality: PASS=%d  WARN=%d  NA=%d\n",
            sum(results$steiger_pass == "PASS", na.rm = TRUE),
            sum(results$steiger_pass == "WARN", na.rm = TRUE),
            sum(results$steiger_pass == "NA", na.rm = TRUE)))

# Focus on Tier 1
tier1_results <- results[results$is_tier1 == TRUE, ]
if (nrow(tier1_results) > 0) {
  cat("\n── Tier 1 Genes (coloc+SuSiE) ──\n")
  for (i in seq_len(nrow(tier1_results))) {
    r <- tier1_results[i]
    cat(sprintf("  %s: SMR b=%.3f, p=%.2e, HEIDI=%s, Steiger=%s, Tier=%s\n",
                r$gene, r$smr_b, r$smr_p,
                r$heidi_pass,
                r$steiger_pass,
                r$sensitivity_tier))
  }
}

# Plot: HEIDI vs SMR significance
results_valid <- results[!is.na(results$heidi_p) & results$heidi_p > 0, ]
if (nrow(results_valid) > 0) {
  p <- ggplot(results_valid, aes(x = -log10(smr_p), y = -log10(heidi_p),
                                  color = sensitivity_tier)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = -log10(3.21e-06), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("A - Robust" = "#1B9E77", "B - HEIDI only" = "#7570B3",
                                  "C - Direction only" = "#E7298A", "D - Weak" = "grey50"),
                       name = "Sensitivity Tier") +
    labs(x = expression(-log[10](p[SMR])),
         y = expression(-log[10](p[HEIDI])),
         title = "MR Sensitivity: HEIDI vs SMR Significance",
         subtitle = sprintf("75 Bonferroni-significant genes")) +
    theme_bw(base_size = 10)

  ggsave("results/phase5c_sensitivity/mr_sensitivity_scatter.pdf", p, width = 8, height = 6)
  ggsave("results/phase5c_sensitivity/mr_sensitivity_scatter.png", p, width = 8, height = 6, dpi = 150)
}

cat("\nDone. Results saved to results/phase5c_sensitivity/\n")
