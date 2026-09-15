#!/usr/bin/env Rscript
# Append MR-PRESSO outlier-corrected estimates to results/phase5d_mr_presso.csv.
# Reads the per-gene .rds objects written by the NbDistribution = 5000 rerun
# (every <=300-instrument gene) and pulls the outlier-corrected MR estimate,
# its standard error and its p-value out of "Main MR results".
suppressMessages({ library(data.table) })
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")

f_csv <- "results/phase5d_mr_presso.csv"
pr <- fread(f_csv)

num <- function(v) {
  v <- as.character(v)[1]
  if (is.na(v) || v == "") return(NA_real_)
  suppressWarnings(as.numeric(sub("^<", "", v)))
}

rows <- rbindlist(lapply(pr$gene, function(g) {
  f <- file.path("results/presso3", paste0(g, ".rds"))
  if (!file.exists(f)) {
    return(data.table(gene = g, presso_corrected_b = NA_real_, presso_corrected_sd = NA_real_,
                      presso_corrected_p = NA_real_, presso_n_outliers_test = NA_integer_))
  }
  x    <- readRDS(f)
  main <- x[["Main MR results"]]
  r    <- x[["MR-PRESSO results"]]
  oc   <- main[main[["MR Analysis"]] == "Outlier-corrected", , drop = FALSE]
  idx  <- r[["Distortion Test"]][["Outliers Indices"]]
  data.table(gene = g,
             presso_corrected_b = if (nrow(oc)) num(oc[["Causal Estimate"]][1]) else NA_real_,
             presso_corrected_sd = if (nrow(oc)) num(oc[["Sd"]][1]) else NA_real_,
             presso_corrected_p = if (nrow(oc)) num(oc[["P-value"]][1]) else NA_real_,
             presso_n_outliers_test = if (is.null(idx)) NA_integer_ else length(idx))
}))

pr <- merge(pr, rows, by = "gene", all.x = TRUE)
ord <- c("gene", "source", "n_inst", "nb_distribution",
         "presso_global_rssobs", "presso_global_p", "presso_global_p_is_bound",
         "presso_outlier_floor", "presso_n_outliers", "presso_distortion_p",
         "presso_distortion_p_is_bound", "presso_corrected_b", "presso_corrected_sd",
         "presso_corrected_p", "presso_n_outliers_test", "seconds", "note")
setcolorder(pr, intersect(ord, names(pr)))
fwrite(pr, f_csv)

chk <- pr[!is.na(presso_n_outliers_test)]
if (nrow(chk) > 0 && any(chk$presso_n_outliers_test != chk$presso_n_outliers)) {
  cat(">>> WARNING: outlier counts in the CSV and in the .rds files disagree\n")
}
print(pr[, .(gene, n_inst, presso_n_outliers, presso_distortion_p,
             presso_corrected_b, presso_corrected_p)])
cat(">>> updated", f_csv, "\n")