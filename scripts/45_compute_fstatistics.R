#!/usr/bin/env Rscript
# Phase 4b: Weak-instrument F-statistics for all MR analyses
# F computed per instrument as (beta_exposure / se_exposure)^2; standard threshold F > 10.
# Covers three MR groups:
#   (1) deCODE   pQTL MR  (5 genes)  -> per-instrument F from results/phase2_decode/*_instruments.csv
#   (2) UKB-PPP  pQTL MR  (4 genes)  -> per-instrument F from results/phase2_pqtl/*_instruments.csv
#   (3) FinnGen  eQTL MR  (8 genes)  -> F ~= Z^2 derived from eQTLGen Z-scores per scripts/12_finngen_replication.R
# Output: results/fstatistics_summary.csv

library(data.table)
library(dplyr)

setDTthreads(16)

proj <- "/ifs1/User/zhouman/project9-v5-crc-atlas"
out_file <- file.path(proj, "results/fstatistics_summary.csv")

results <- list()

compute_f_summary <- function(gene, source, n, F_values) {
  # F_values: numeric vector of per-instrument F
  Fv <- F_values[!is.na(F_values) & is.finite(F_values) & F_values > 0]
  n_eff <- length(Fv)
  if (n_eff == 0) {
    return(data.frame(
      gene = gene, source = source, n_instruments = n, n_valid_f = 0,
      f_mean = NA_real_, f_median = NA_real_, f_min = NA_real_, f_max = NA_real_,
      n_f_gt_10 = 0, pct_f_gt_10 = NA_real_, n_f_gt_100 = 0,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    gene = gene, source = source, n_instruments = n, n_valid_f = n_eff,
    f_mean = mean(Fv), f_median = median(Fv), f_min = min(Fv), f_max = max(Fv),
    n_f_gt_10 = sum(Fv > 10), pct_f_gt_10 = 100 * sum(Fv > 10) / n_eff,
    n_f_gt_100 = sum(Fv > 100),
    stringsAsFactors = FALSE
  )
}

# === 1. deCODE pQTL MR (5 genes) ===
cat("=== deCODE pQTL MR F-statistics ===\n")
decode_genes <- c("CCM2", "LIMA1", "STAT6", "TNFRSF1A", "TNFRSF1B")
for (g in decode_genes) {
  f <- file.path(proj, sprintf("results/phase2_decode/%s_instruments.csv", g))
  if (!file.exists(f)) {
    cat(sprintf("  [%s] MISSING %s\n", g, f))
    results[[length(results)+1]] <- data.frame(
      gene=g, source="deCODE", n_instruments=NA, n_valid_f=0,
      f_mean=NA, f_median=NA, f_min=NA, f_max=NA,
      n_f_gt_10=0, pct_f_gt_10=NA, n_f_gt_100=0, stringsAsFactors=FALSE
    )
    next
  }
  d <- fread(f, header=TRUE)
  # deCODE schema: pqtl.rsid, pqtl.beta, pqtl.se
  if (!"pqtl.beta" %in% names(d) || !"pqtl.se" %in% names(d)) {
    cat(sprintf("  [%s] schema missing pqtl.beta/pqtl.se\n", g))
    results[[length(results)+1]] <- data.frame(
      gene=g, source="deCODE", n_instruments=nrow(d), n_valid_f=0,
      f_mean=NA, f_median=NA, f_min=NA, f_max=NA,
      n_f_gt_10=0, pct_f_gt_10=NA, n_f_gt_100=0, stringsAsFactors=FALSE
    )
    next
  }
  d[, F := (pqtl.beta / pqtl.se)^2]
  cat(sprintf("  [%s] %d instruments\n", g, nrow(d)))
  results[[length(results)+1]] <- compute_f_summary(g, "deCODE", nrow(d), d$F)
}

# === 2. UKB-PPP pQTL MR (4 genes; TET2 has 0 instruments) ===
cat("\n=== UKB-PPP pQTL MR F-statistics ===\n")
ukb_genes <- c("CDKN1A", "LTA", "NCF2", "TNF")  # TET2 excluded (no cis-pQTL)
for (g in ukb_genes) {
  f <- file.path(proj, sprintf("results/phase2_pqtl/%s_instruments.csv", g))
  if (!file.exists(f)) {
    cat(sprintf("  [%s] MISSING %s\n", g, f))
    results[[length(results)+1]] <- data.frame(
      gene=g, source="UKB-PPP", n_instruments=NA, n_valid_f=0,
      f_mean=NA, f_median=NA, f_min=NA, f_max=NA,
      n_f_gt_10=0, pct_f_gt_10=NA, n_f_gt_100=0, stringsAsFactors=FALSE
    )
    next
  }
  d <- fread(f, header=TRUE)
  if (!"pqtl.beta" %in% names(d) || !"pqtl.se" %in% names(d)) {
    cat(sprintf("  [%s] schema missing pqtl.beta/pqtl.se\n", g))
    results[[length(results)+1]] <- data.frame(
      gene=g, source="UKB-PPP", n_instruments=nrow(d), n_valid_f=0,
      f_mean=NA, f_median=NA, f_min=NA, f_max=NA,
      n_f_gt_10=0, pct_f_gt_10=NA, n_f_gt_100=0, stringsAsFactors=FALSE
    )
    next
  }
  d[, F := (pqtl.beta / pqtl.se)^2]
  cat(sprintf("  [%s] %d instruments\n", g, nrow(d)))
  results[[length(results)+1]] <- compute_f_summary(g, "UKB-PPP", nrow(d), d$F)
}

# === 3. FinnGen eQTL MR (8 genes) — F ~= Z^2 ===
cat("\n=== FinnGen eQTL MR F-statistics (F ~= Z^2) ===\n")
fg <- fread(file.path(proj, "results/phase4_finngen_replication.csv"),
            select=c("gene","ensg","chr","pos"))
cat(sprintf("  FinnGen genes: %d\n", nrow(fg)))

eqtl_file <- file.path(proj,
  "data/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
fg_chrs <- sort(unique(as.integer(fg$chr)))
cat(sprintf("  Loading eQTLGen for chr: %s ...\n", paste(fg_chrs, collapse=",")))

eqtl_all <- fread(eqtl_file, header=TRUE,
                  select=c("SNP","SNPChr","SNPPos","Zscore","Gene","Pvalue"),
                  tmpdir="/tmp/zhouman_tmp")
eqtl_all <- eqtl_all[SNPChr %in% fg_chrs]
eqtl_all[, SNPChr := as.integer(SNPChr)]
cat(sprintf("  eQTL records loaded: %d\n", nrow(eqtl_all)))

CIS_WINDOW <- 1e6
for (i in 1:nrow(fg)) {
  gene <- fg$gene[i]; ensg <- fg$ensg[i]
  gchr <- as.integer(fg$chr[i]); gpos <- as.integer(fg$pos[i])

  # Per script 12: cis ±1Mb, then instruments with Pvalue < 5e-4
  eqtl_gene <- eqtl_all[Gene == ensg & SNPChr == gchr &
                        SNPPos >= gpos - CIS_WINDOW & SNPPos <= gpos + CIS_WINDOW]
  inst <- eqtl_gene[Pvalue < 5e-4]

  if (nrow(inst) == 0) {
    cat(sprintf("  [%s] no instruments (p<5e-4) in cis window\n", gene))
    results[[length(results)+1]] <- data.frame(
      gene=gene, source="FinnGen", n_instruments=0, n_valid_f=0,
      f_mean=NA, f_median=NA, f_min=NA, f_max=NA,
      n_f_gt_10=0, pct_f_gt_10=NA, n_f_gt_100=0, stringsAsFactors=FALSE
    )
    next
  }
  Fv <- inst$Zscore^2
  cat(sprintf("  [%s] %d instruments (p<5e-4)\n", gene, nrow(inst)))
  results[[length(results)+1]] <- compute_f_summary(gene, "FinnGen", nrow(inst), Fv)
}

# === 4. Output ===
cat("\n\n=== F-statistics Summary ===\n")
out <- rbindlist(results)
out <- out[order(source, -f_mean)]

# Overall weak-instrument assessment
n_total <- nrow(out)
n_valid <- sum(!is.na(out$f_mean))
n_weak <- sum(out$pct_f_gt_10 < 100, na.rm=TRUE)  # genes where not ALL instruments F>10

cat(sprintf("Genes tested: %d;  with valid F: %d\n", n_total, n_valid))
cat(sprintf("Genes with any instrument F<=10 (potential weak-IV): %d\n", n_weak))
cat(sprintf("Genes with ALL instruments F>100 (very strong): %d\n",
            sum(out$n_f_gt_100 > 0 & !is.na(out$n_f_gt_10), na.rm=TRUE)))

print(out)

fwrite(out, out_file)
cat(sprintf("\nSaved: %s\n", out_file))
