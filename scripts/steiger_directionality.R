#!/usr/bin/env Rscript
# Steiger directionality analysis of the protein-layer cis-pQTL instruments.
# Added 2026-09-21 (Supplementary Table S10). Requires TwoSampleMR (>= 0.5).
# Inputs : results/phase2_decode/*_instruments.csv, results/phase2_pqtl/*_instruments.csv
# Outputs: results/phase6_steiger/phase6_steiger_by_gene.csv
#          results/phase6_steiger/phase6_steiger_by_snp.csv
# Phase 6 (v2): Steiger directionality test for the protein-layer (cis-pQTL) MR
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })

PROJ   <- "/ifs1/User/zhouman/project9-v5-crc-atlas"
OUT    <- file.path(PROJ, "results/phase6_steiger")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
GWAS_N <- 185616
CRC    <- "Colorectal cancer (GCST90255675)"

panels <- list(
  list(name = "deCODE",  dir = file.path(PROJ, "results/phase2_decode"), n = 35278,
       genes = c("CCM2","LIMA1","STAT6","TNFRSF1A","TNFRSF1B")),
  list(name = "UKB-PPP", dir = file.path(PROJ, "results/phase2_pqtl"),   n = 33533,
       genes = c("CDKN1A","LTA","NCF2","TNF"))
)

r_from_bsen <- function(b, se, n) { Fv <- (b/se)^2; r <- sqrt(Fv/(Fv + n - 2)); r[!is.finite(r)] <- NA_real_; r }
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p == 0, "<1e-300", formatC(p, format = "g", digits = 3)))

per_gene <- list(); per_snp <- list()

for (pn in panels) for (g in pn$genes) {
  f <- file.path(pn$dir, paste0(g, "_instruments.csv"))
  if (!file.exists(f)) { cat("MISSING:", f, "\n"); next }
  d <- fread(f)
  snpcol <- if ("gwas.rsid" %in% names(d)) "gwas.rsid" else "gwas.snp"
  d <- d[!is.na(pqtl.beta) & !is.na(pqtl.se) & !is.na(gwas.beta) & !is.na(gwas.se)]

  r.exp <- r_from_bsen(d$pqtl.beta, d$pqtl.se, pn$n)
  r.out <- r_from_bsen(d$gwas.beta, d$gwas.se, GWAS_N)
  fwd   <- r.exp^2 > r.out^2

  per_snp[[paste(pn$name,g)]] <- data.table(panel = pn$name, gene = g, snp = d[[snpcol]],
    p_pqtl = d$pqtl.pval, p_crc = d$gwas.pval,
    r2_pqtl = r.exp^2, r2_crc = r.out^2, correct_direction = fwd)

  # aggregate test (as implemented in TwoSampleMR::directionality_test)
  rsum_exp <- sqrt(sum(r.exp^2, na.rm = TRUE)); rsum_out <- sqrt(sum(r.out^2, na.rm = TRUE))
  p_agg <- tryCatch({
    rt <- suppressWarnings(psych::r.test(n = pn$n, n2 = GWAS_N, r12 = rsum_exp, r34 = rsum_out))
    pnorm(-abs(rt[["z"]])) * 2
  }, error = function(e) NA_real_)
  if (length(p_agg) != 1 || !is.finite(p_agg)) {
    p_agg2 <- tryCatch({ rt <- suppressWarnings(psych::r.test(n = pn$n, n2 = GWAS_N, r12 = rsum_exp, r34 = rsum_out)); rt[["p"]] }, error = function(e) NA_real_)
    if (!is.null(p_agg2) && length(p_agg2) == 1 && is.finite(p_agg2) && p_agg2 < 1) p_agg <- p_agg2
  }

  # lead-variant test (LD-free single-instrument version)
  i <- which.min(d$pqtl.pval)
  lr_exp <- r.exp[i]; lr_out <- r.out[i]
  p_lead <- tryCatch({
    rt <- suppressWarnings(psych::r.test(n = pn$n, n2 = GWAS_N, r12 = lr_exp, r34 = lr_out))
    v <- pnorm(-abs(rt[["z"]])) * 2
    if (!is.finite(v)) { v2 <- rt[["p"]]; if (length(v2) == 1 && is.finite(v2)) v2 else NA_real_ } else v
  }, error = function(e) NA_real_)

  per_gene[[paste(pn$name,g)]] <- data.table(
    panel = pn$name, gene = g, n_instruments = nrow(d),
    n_correct = sum(fwd, na.rm = TRUE), pct_correct = round(100*mean(fwd, na.rm = TRUE), 1),
    median_r2_pqtl = median(r.exp^2, na.rm = TRUE), median_r2_crc = median(r.out^2, na.rm = TRUE),
    agg_r2_pqtl = rsum_exp^2, agg_r2_crc = rsum_out^2,
    steiger_p_aggregate = p_agg,
    lead_snp = d[[snpcol]][i], lead_p_pqtl = d$pqtl.pval[i], lead_p_crc = d$gwas.pval[i],
    lead_r2_pqtl = lr_exp^2, lead_r2_crc = lr_out^2, steiger_p_lead = p_lead,
    correct_direction = ifelse(is.finite(p_lead), lr_exp > lr_out, sum(fwd, na.rm=TRUE) == nrow(d)))

  cat(sprintf("%-8s %-9s n=%5d  ok=%5d (%.1f%%)  lead=%s  p_lead=%s  p_agg=%s\n",
              pn$name, g, nrow(d), sum(fwd, na.rm = TRUE), 100*mean(fwd, na.rm = TRUE),
              d[[snpcol]][i], fmt_p(p_lead), fmt_p(p_agg)))
}

pg <- rbindlist(per_gene, fill = TRUE); ps <- rbindlist(per_snp, fill = TRUE)
fwrite(pg, file.path(OUT, "phase6_steiger_by_gene.csv"))
fwrite(ps, file.path(OUT, "phase6_steiger_by_snp.csv"))

cat("\n=== BY GENE ===\n"); print(as.data.frame(pg))
cat(sprintf("\nSNP-level: %d pairs, %d correct (%.4f)\n", nrow(ps), sum(ps$correct_direction, na.rm=TRUE), mean(ps$correct_direction, na.rm=TRUE)))
cat("WROTE by_gene + by_snp in", OUT, "\n")
cat("STEIGERDONE\n")
