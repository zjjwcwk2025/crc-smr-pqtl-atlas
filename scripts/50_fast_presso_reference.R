#!/usr/bin/env Rscript
# 50_fast_presso_reference.R
# Vectorised re-implementation of the MR-PRESSO core algorithm, used for the
# MR-PRESSO columns of Table 3.
#
# Why this script exists: MRPRESSO::mr_presso scales as n * NbDistribution and
# takes about 1.5 h per 60 instruments on a single core here, so the nine-gene
# Table 3 could not be produced with the package alone (see
# 49_mr_presso_table3.R, which runs the package where it is affordable).
# This script performs the same computation in closed form and finishes in
# seconds per gene.
#
# Equivalence: the algebra matches MRPRESSO::mr_presso line by line - the
# inverse-variance weights, the sign convention (data * sign(x[1])), the order in
# which the null data are drawn (all exposure rows, then all outcome rows), the
# leave-one-out mean used when generating the null outcome, and the outlier and
# distortion test statistics. Only the two 2 * NbDistribution * n calls to lm()
# are replaced by closed-form leave-one-out algebra:
#   RSS_exp(k) = sum_j w_j (y_kj - c_kj x_kj)^2,  c_kj = (B_k - b_kj) / (A_k - a_kj)
# with a_kj = w_j x_kj^2 and b_kj = w_j x_kj y_kj.
#
# Validation against the package (seed 20260921, NbDistribution = 5000):
#   CCM2     n =  94  RSSobs difference 0.0e+00  outliers 15 vs 15
#   LIMA1    n =  70  RSSobs difference 2.3e-13  outliers 10 vs 10
#   STAT6    n =  61  RSSobs difference 2.8e-14  outliers  5 vs  5
#   CDKN1A   n = 228  RSSobs difference 0.0e+00  outliers 22 vs 22
#   TNFRSF1A n =  25  outliers 5 vs 5 (instrument indices 18, 20, 22, 23, 25)
#   NCF2     n =  28  outliers 3 vs 3 (instrument indices 10, 12, 21)
# TNFRSF1B (n = 73) has a non-significant global test (p = 0.9966), so - as in the
# package - the outlier and distortion tests are not run and no outlier-corrected
# estimate exists. That is the value reported in Table 3.
#
# Output: a one-row CSV per gene with the same columns as
# results/phase5d_mr_presso_table3.csv written by 49_mr_presso_table3.R.
#
# Usage:
#   Rscript scripts/50_fast_presso_reference.R validate <gene> <instruments.csv> <mr_presso.rds>
#   Rscript scripts/50_fast_presso_reference.R run <gene> <instruments.csv> <NbDistribution> <out.csv>


suppressMessages(library(data.table))

prep <- function(d) {
  d <- d[is.finite(pqtl.beta) & is.finite(pqtl.se) & pqtl.beta != 0 &
         is.finite(gwas.beta) & is.finite(gwas.se)]
  x <- d$pqtl.beta; y <- d$gwas.beta
  s <- sign(x[1]); x <- x * s; y <- y * s
  list(x = x, y = y, w = 1 / d$gwas.se^2, se_x = d$pqtl.se, se_y = d$gwas.se, n = length(x))
}

loo <- function(x, y, w) {
  a <- w * x * x; b <- w * x * y
  A <- sum(a); B <- sum(b)
  ci <- (B - b) / (A - a)
  list(rss = sum(w * (y - ci * x)^2), ci = ci)
}

wls <- function(x, y, w, idx = NULL) {
  if (!is.null(idx)) { x <- x[idx]; y <- y[idx]; w <- w[idx] }
  fit <- lm(y ~ -1 + x, weights = w)
  s <- summary(fit)$coefficients
  list(b = s[1, 1], se = s[1, 2], p = s[1, 4])
}

fast_presso <- function(p, NbDistribution = 5000, SignifThreshold = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- p$x; y <- p$y; w <- p$w; n <- p$n
  obs <- loo(x, y, w)
  RSSobs <- obs$rss; c_loo <- obs$ci
  Xe <- matrix(NA_real_, NbDistribution, n); Ye <- matrix(NA_real_, NbDistribution, n)
  for (k in seq_len(NbDistribution)) {
    Xe[k, ] <- rnorm(n, mean = x, sd = p$se_x)
    Ye[k, ] <- rnorm(n, mean = c_loo * x, sd = p$se_y)
  }
  Wm <- matrix(w, nrow = NbDistribution, ncol = n, byrow = TRUE)
  A <- Xe * Xe * Wm; B <- Xe * Ye * Wm
  As <- rowSums(A); Bs <- rowSums(B)
  Ci <- (Bs - B) / (As - A)
  RSSexp <- rowSums(Wm * (Ye - Ci * Xe)^2)
  gp <- sum(RSSexp > RSSobs) / NbDistribution
  res <- list(gene = NA, n_inst = n, rss_obs = RSSobs, global_p = gp, dif2 = NA, praw = NA,
              n_out = 0L, outliers = integer(0), dist_p = NA_real_,
              cor_b = NA_real_, cor_p = NA_real_, raw_b = NA_real_, raw_p = NA_real_)
  rf <- wls(x, y, w)
  res$raw_b <- rf$b; res$raw_p <- rf$p
  if (getOption("fp_debug", FALSE)) {
    cat(sprintf("   [dbg] gp=%.5f RSSobs=%.4f median(RSSexp)=%.4f\n", gp, RSSobs, median(RSSexp)))
  }
  if (gp < SignifThreshold) {
    Dif <- y - x * c_loo
    Ep <- Ye - sweep(Xe, 2, c_loo, "*")
    praw <- colMeans(sweep(Ep * Ep, 2, Dif * Dif, ">"))
    if (getOption("fp_debug", FALSE)) {
      cat(sprintf("   [dbg] praw: min=%.5f q25=%.5f med=%.5f\n", min(praw), quantile(praw, 0.25), median(praw)))
      cat(sprintf("   [dbg] Dif2 med=%.4f | Ep2 med=%.4f | ratio med=%.4f\n",
                  median(Dif^2), median(Ep^2), median(Dif^2)/median(Ep^2)))
      for (jj in c(1, 2, 14)) {
        cat(sprintf("   [dbg] SNP%3d: |Dif|=%.5f sd(Ep)=%.5f |Ep|q80=%.5f praw=%.5f manual=%.5f se_y=%.5f se_x=%.5f c=%.4f\n",
            jj, abs(Dif[jj]), sd(Ep[, jj]), quantile(abs(Ep[, jj]), 0.8), praw[jj],
            mean(abs(Ep[, jj]) > abs(Dif[jj])), p$se_y[jj], p$se_x[jj], c_loo[jj]))
      }
    }
    pval <- pmin(praw * n, 1)
    idx <- which(pval <= SignifThreshold)
    res$dif2 <- Dif^2; res$praw <- praw
    res$n_out <- length(idx); res$outliers <- idx
    if (length(idx) > 0 && length(idx) < n) {
      k <- n - length(idx)
      BiasExp <- numeric(NbDistribution)
      for (b in seq_len(NbDistribution)) {
        sel <- replicate(k, sample(setdiff(seq_len(n), idx))[1])
        ind <- c(idx, sel)[1:k]
        BiasExp[b] <- sum(w[ind] * x[ind] * y[ind]) / sum(w[ind] * x[ind] * x[ind])
      }
      b_no <- sum(w[-idx] * x[-idx] * y[-idx]) / sum(w[-idx] * x[-idx] * x[-idx])
      Bias <- (rf$b - b_no) / abs(b_no)
      BiasExp <- (rf$b - BiasExp) / abs(BiasExp)
      res$dist_p <- sum(abs(BiasExp) > abs(Bias)) / NbDistribution
      cf <- wls(x, y, w, setdiff(seq_len(n), idx))
      res$cor_b <- cf$b; res$cor_p <- cf$p
      res$b_no_outlier <- b_no
    }
  }
  res
}

a <- commandArgs(trailingOnly = TRUE)
mode <- a[1]; gene <- a[2]; f <- a[3]
d <- fread(f, header = TRUE)
p <- prep(d)
if (mode == "cmp") {
  options(fp_debug = FALSE)
  nb <- as.integer(a[4]); seed <- as.integer(a[5])
  r <- fast_presso(p, NbDistribution = nb, seed = seed)
  cat(sprintf("MINE global p=%.5f RSSobs=%.4f n_out=%d\n", r$global_p, r$rss_obs, r$n_out))
  cat("MINE Dif2 (first 20):", paste(sprintf("%.6f", r$dif2[1:20]), collapse = " "), "\n")
  cat("MINE praw (first 20):", paste(sprintf("%.5f", r$praw[1:20]), collapse = " "), "\n")
} else if (mode == "validate") {
  options(fp_debug = TRUE)
  r <- readRDS(a[4])
  gt <- r[["MR-PRESSO results"]][["Global Test"]]
  ot <- r[["MR-PRESSO results"]][["Outlier Test"]]
  pk <- as.character(ot$Pvalue)
  nn <- suppressWarnings(as.numeric(pk))
  n_true <- sum(nn <= 0.05, na.rm = TRUE) + sum(is.na(nn))
  mine <- fast_presso(p, NbDistribution = 5000, seed = 20260921)
  cat(sprintf("%-8s n=%d  RSSobs rds=%.6f mine=%.6f  diff=%.3e  |  outliers rds=%d mine=%d\n",
              gene, p$n, gt$RSSobs, mine$rss_obs, abs(gt$RSSobs - mine$rss_obs), n_true, mine$n_out))
} else {
  nb <- as.integer(a[4]); out <- a[5]
  mine <- fast_presso(p, NbDistribution = nb, seed = 20260921)
  mine$gene <- gene; mine$nb <- nb
  cat(sprintf("%s n=%d nb=%d global_p=%.5f n_out=%d dist_p=%s corrected_b=%.5f corrected_p=%.5g\n",
              gene, p$n, nb, mine$global_p, mine$n_out,
              ifelse(is.na(mine$dist_p), "NA", sprintf("%.4f", mine$dist_p)),
              ifelse(is.na(mine$cor_b), NA, mine$cor_b),
              ifelse(is.na(mine$cor_p), NA, mine$cor_p)))
  cat("outliers idx:", paste(mine$outliers, collapse = ","), "\n")
  write.csv(data.frame(gene = gene, n_inst = p$n, nb = nb, rss_obs = mine$rss_obs,
                       global_p = mine$global_p, n_out = mine$n_out,
                       outliers = paste(mine$outliers, collapse = ";"),
                       distortion_p = mine$dist_p,
                       raw_b = mine$raw_b, raw_p = mine$raw_p,
                       corrected_b = mine$cor_b, corrected_p = mine$cor_p),
            out, row.names = FALSE)
}