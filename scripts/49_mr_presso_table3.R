#!/usr/bin/env Rscript
# 49_mr_presso_table3.R
# MR-PRESSO diagnostics for every gene with >= 2 genome-wide significant cis-pQTL
# instruments, i.e. the nine genes listed in Table 3 of the manuscript.
#
# Why this script exists: 48_mr_methods_pleiotropy.R calls mr_presso() with the
# package default NbDistribution = 1000. At that draw count the smallest
# resolvable outlier-test p-value is n/1000 (0.061-0.228 for these genes), which is
# above 0.05, so the outlier test cannot flag anything. The manuscript reports
# NbDistribution = 5000, for which n/NbDistribution <= 0.046 for every gene with
# n <= 228 instruments. This script reproduces those numbers.
#
# Runtime: mr_presso() is roughly linear in n * NbDistribution. On a single core
# this is about 1.5 h per 60 instruments, so the seven feasible genes take several
# hours in total. Run them as separate jobs. TNF (n = 2288) and LTA (n = 10213)
# are not feasible: MR-PRESSO needs more than n/0.046 draws, i.e. > 5.0e4 for TNF
# and > 2.2e5 for LTA; both are recorded as skipped in Table 3.
#
# Output: results/phase5d_mr_presso_table3.csv
#   gene, source, n_inst, nb_distribution, presso_global_p, presso_n_outliers,
#   presso_n_outlier_p_na, presso_distortion_p, presso_raw_b, presso_raw_p,
#   presso_corrected_b, presso_corrected_p, seconds
suppressMessages({ library(data.table); library(MRPRESSO) })
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")

NB <- 5000          # draw count, matching the manuscript
MAX_N <- 228        # largest instrument count for which n / NB <= 0.046

genes <- list(
  list(g = "TNFRSF1A", src = "deCODE",  f = "results/phase2_decode/TNFRSF1A_instruments.csv"),
  list(g = "NCF2",     src = "UKB-PPP", f = "results/phase2_pqtl/NCF2_instruments.csv"),
  list(g = "STAT6",    src = "deCODE",  f = "results/phase2_decode/STAT6_instruments.csv"),
  list(g = "LIMA1",    src = "deCODE",  f = "results/phase2_decode/LIMA1_instruments.csv"),
  list(g = "TNFRSF1B", src = "deCODE",  f = "results/phase2_decode/TNFRSF1B_instruments.csv"),
  list(g = "CCM2",     src = "deCODE",  f = "results/phase2_decode/CCM2_instruments.csv"),
  list(g = "CDKN1A",   src = "UKB-PPP", f = "results/phase2_pqtl/CDKN1A_instruments.csv"),
  list(g = "TNF",      src = "UKB-PPP", f = "results/phase2_pqtl/TNF_instruments.csv"),
  list(g = "LTA",      src = "UKB-PPP", f = "results/phase2_pqtl/LTA_instruments.csv")
)

run_one <- function(z) {
  d <- fread(z$f, header = TRUE)
  d <- d[is.finite(pqtl.beta) & is.finite(pqtl.se) & pqtl.beta != 0 &
         is.finite(gwas.beta) & is.finite(gwas.se)]
  res <- data.table(gene = z$g, source = z$src, n_inst = nrow(d),
                    nb_distribution = NB, presso_global_p = NA_real_,
                    presso_n_outliers = NA_integer_, presso_n_outlier_p_na = NA_integer_,
                    presso_distortion_p = NA_real_, presso_raw_b = NA_real_,
                    presso_raw_p = NA_real_, presso_corrected_b = NA_real_,
                    presso_corrected_p = NA_real_, seconds = NA_real_)
  if (nrow(d) > MAX_N) { cat(sprintf("skip %s: n=%d\n", z$g, nrow(d))); return(res) }
  df <- data.frame(BetaOutcome = d$gwas.beta, BetaExposure = d$pqtl.beta,
                   SdOutcome = d$gwas.se, SdExposure = d$pqtl.se)
  cat(sprintf("=== %s n=%d start %s ===\n", z$g, nrow(d), Sys.time())); flush.console()
  t0 <- Sys.time()
  pr <- mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure",
                  SdOutcome = "SdOutcome", SdExposure = "SdExposure", data = df,
                  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                  NbDistribution = NB, SignifThreshold = 0.05)
  res[, seconds := round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)]
  res[, presso_global_p := as.numeric(pr$`MR-PRESSO results`$`Global Test`$Pvalue)]
  pv <- suppressWarnings(as.numeric(pr$`MR-PRESSO results`$`Outlier Test`$Pvalue))
  res[, presso_n_outliers := sum(pv < 0.05, na.rm = TRUE)]
  res[, presso_n_outlier_p_na := sum(is.na(pv))]
  ds <- pr$`MR-PRESSO results`$`Distortion Test`
  if (!is.null(ds)) res[, presso_distortion_p := suppressWarnings(as.numeric(ds$Pvalue)[1])]
  mn <- pr$`Main MR results`
  if (is.data.frame(mn) && "MR Analysis" %in% names(mn)) {
    k1 <- which(mn$`MR Analysis` == "Raw"); k2 <- which(mn$`MR Analysis` != "Raw")
    if (length(k1)) { res[, presso_raw_b := as.numeric(mn$`Causal Estimate`[k1[1]])]
                      res[, presso_raw_p := as.numeric(mn$`P-value`[k1[1]])] }
    if (length(k2)) { res[, presso_corrected_b := as.numeric(mn$`Causal Estimate`[k2[1]])]
                      res[, presso_corrected_p := as.numeric(mn$`P-value`[k2[1]])] }
  }
  print(as.data.frame(res)); flush.console()
  res
}

out <- rbindlist(lapply(genes, run_one), fill = TRUE)
fwrite(out, "results/phase5d_mr_presso_table3.csv")
cat("\nWROTE results/phase5d_mr_presso_table3.csv\n")