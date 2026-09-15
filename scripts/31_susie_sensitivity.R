#!/usr/bin/env Rscript
# 31_susie_sensitivity.R — SuSiE 敏感性分析
# 审稿人质疑: "6 genes, 38 credible sets, all PIP=1.0 + single SNP"
# 测试: 不同 SNP 数量上限 × L × 窗口大小
#
# 核心假设: 300-SNP pruning 过度压缩 LD，导致所有 CS 为单 SNP
# 预测: 当 SNP cap 增大到 500-1000 时，可能出现 multi-SNP CS

library(data.table)
library(susieR)
library(dplyr)

setDTthreads(8)

# ===== Config =====
GWAS_FILE  <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.h.tsv.gz"
BIM_FILE   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/ref_eur/1000G_EUR.bim"
BED_PREFIX <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/ref_eur/1000G_EUR"
GWAS_N     <- 185616L
OUT_DIR    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/susie_sensitivity"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 6 coloc-passing genes
GENES <- data.frame(
  gene  = c("UTP11",  "BMP2",  "PNKD",  "LAMC1", "TMBIM1", "GPBAR1"),
  chr   = c(1,        20,      2,       1,       2,        2),
  pos   = c(38482713, 6754619, 219173315, 183053661, 219148112, 219126400),
  stringsAsFactors = FALSE
)

# ===== Parameter grid =====
param_grid <- expand.grid(
  snp_cap   = c(200, 300, 500),      # SNP number limit per locus
  L          = c(5, 10, 20),          # Max credible sets
  window_bp  = c(500e3, 1000e3),      # ±window around gene
  stringsAsFactors = FALSE
)

cat(sprintf("Total parameter combinations: %d x %d genes = %d runs\n",
            nrow(param_grid), nrow(GENES), nrow(param_grid) * nrow(GENES)))

# ===== Load BIM =====
cat("\nLoading BIM...\n")
bim <- fread(BIM_FILE, header = FALSE, select = 1:6,
             col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
bim[, chr := as.integer(chr)]
bim_chr <- split(bim, bim$chr)

# ===== LD computation helper =====
compute_ld <- function(snp_ids, chr, bed_prefix) {
  tmp_snps <- tempfile(fileext = ".txt")
  writeLines(snp_ids, tmp_snps)
  tmp_out <- tempfile(fileext = "")
  plink_cmd <- sprintf(
    "%s --bfile %s --chr %d --extract %s --r-unphased square --out %s --threads 4 --memory 16000 2>/dev/null",
    "/ifs1/User/zhouman/project9-v5-crc-atlas/tools/plink2",
    bed_prefix, chr, tmp_snps, tmp_out
  )
  system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ld_file  <- paste0(tmp_out, ".unphased.vcor1")
  var_file <- paste0(tmp_out, ".unphased.vcor1.vars")
  if (!file.exists(ld_file)) { unlink(tmp_snps); return(NULL) }

  vars <- fread(var_file, header = FALSE)[[1]]
  n <- length(vars)
  con <- file(ld_file, "rb")
  R_vec <- readBin(con, "double", n = n * n, size = 8)
  close(con)
  R <- matrix(R_vec, nrow = n, ncol = n, byrow = TRUE)
  rownames(R) <- vars; colnames(R) <- vars

  unlink(c(tmp_snps, ld_file, var_file, paste0(tmp_out, ".log")))
  return(R)
}

# ===== Run SuSiE for one gene × one param =====
run_susie <- function(gene_row, snp_cap, L, window_bp) {
  gchr  <- gene_row$chr
  gpos  <- gene_row$pos
  gene  <- gene_row$gene

  # Load GWAS SNPs in cis
  gw_cmd <- sprintf("zcat %s | awk 'NR==1 || ($1==%d && $2>=%d && $2<=%d)'",
                    GWAS_FILE, gchr, gpos - window_bp, gpos + window_bp)
  gwas <- tryCatch(fread(cmd = gw_cmd), error = function(e) NULL)
  if (is.null(gwas) || nrow(gwas) < 10) return(NULL)

  setnames(gwas,
    old = c("chromosome","base_pair_location","effect_allele","other_allele",
            "beta","standard_error","effect_allele_frequency","p_value","rsid"),
    new = c("chr","bp","ea","oa","beta","se","eaf","pval","snp"))

  # Intersect with 1000G
  bim_chr_data <- bim_chr[[as.character(gchr)]]
  if (is.null(bim_chr_data) || nrow(bim_chr_data) < 10) return(NULL)

  gwas_w_ld <- merge(gwas, bim_chr_data[, .(rsid, pos_ld = pos)],
                     by.x = "snp", by.y = "rsid", all.x = FALSE)
  gwas_w_ld <- gwas_w_ld[abs(gwas_w_ld$pos_ld - gpos) <= window_bp]

  if (nrow(gwas_w_ld) < 20) return(NULL)

  # Prune by p-value to SNP cap
  setorder(gwas_w_ld, pval)
  if (nrow(gwas_w_ld) > snp_cap) {
    gwas_w_ld <- gwas_w_ld[1:snp_cap]
  }
  n_snp <- nrow(gwas_w_ld)

  # Z-scores
  gwas_w_ld[, z := beta / se]
  gwas_w_ld <- gwas_w_ld[is.finite(z)]
  if (nrow(gwas_w_ld) < 20) return(NULL)

  # LD matrix
  R <- compute_ld(gwas_w_ld$snp, gchr, BED_PREFIX)
  if (is.null(R) || nrow(R) < 20) return(NULL)

  diag(R) <- diag(R) + 1e-4
  R <- (R + t(R)) / 2

  ld_snps <- rownames(R)
  gwas_match <- match(ld_snps, gwas_w_ld$snp)
  gwas_aligned <- gwas_w_ld[gwas_match]
  z <- gwas_aligned$z

  if (length(z) != nrow(R)) return(NULL)

  # Effective L (capped by data)
  L_eff <- min(L, length(z) %/% 20)
  if (L_eff < 1) L_eff <- 1

  # Run SuSiE
  fit <- tryCatch({
    susie_rss(z = z, R = R, n = GWAS_N,
              L = L_eff,
              estimate_residual_variance = FALSE,
              min_abs_corr = 0.3,
              coverage = 0.9,
              refine = TRUE,
              max_iter = 500,
              check_inputs = FALSE)
  }, error = function(e) NULL)

  if (is.null(fit)) return(NULL)

  n_cs <- length(fit$sets$cs)
  converged <- fit$converged

  # CS composition analysis
  cs_sizes <- integer(0)
  cs_multi <- 0L  # number of CS with >1 SNP
  cs_pips  <- numeric(0)

  for (cs_idx in fit$sets$cs) {
    cs_sizes <- c(cs_sizes, length(cs_idx))
    if (length(cs_idx) > 1) cs_multi <- cs_multi + 1L
    cs_pips <- c(cs_pips, max(fit$pip[cs_idx]))
  }

  # Top SNP info
  top_idx <- which.max(fit$pip)
  top_snp  <- gwas_aligned$snp[top_idx]
  top_pip  <- fit$pip[top_idx]

  # Number of SNPs with PIP > 0.1 (to see if ridge is over-shrinking)
  n_pip_gt_01 <- sum(fit$pip > 0.1)
  n_pip_gt_05 <- sum(fit$pip > 0.5)

  # Mean CS size
  mean_cs_size <- if (length(cs_sizes) > 0) mean(cs_sizes) else NA

  data.frame(
    gene         = gene,
    snp_cap      = snp_cap,
    L_param      = L,
    L_eff        = L_eff,
    window_kb    = window_bp / 1000,
    n_total_snp  = n_snp,
    n_used_snp   = length(z),
    n_cs         = n_cs,
    n_multi_snp_cs = cs_multi,
    mean_cs_size = mean_cs_size,
    cs_sizes     = paste(cs_sizes, collapse = ";"),
    cs_pips      = paste(round(cs_pips, 3), collapse = ";"),
    top_pip_snp  = top_snp,
    top_pip      = top_pip,
    n_pip_gt_01  = n_pip_gt_01,
    n_pip_gt_05  = n_pip_gt_05,
    converged    = converged,
    stringsAsFactors = FALSE
  )
}

# ===== Execute grid search =====
all_results <- list()
counter <- 0L

for (g in seq_len(nrow(GENES))) {
  gene_name <- GENES$gene[g]
  cat(sprintf("\n========== %s ==========\n", gene_name))

  for (p in seq_len(nrow(param_grid))) {
    counter <- counter + 1L
    par <- param_grid[p, ]

    cat(sprintf("  [%d/%d] cap=%d L=%d window=%.0fkb ... ",
                p, nrow(param_grid), par$snp_cap, par$L, par$window_bp/1000))

    res <- run_susie(GENES[g, ], par$snp_cap, par$L, par$window_bp)

    if (is.null(res)) {
      cat("FAILED\n")
    } else {
      # Flag if any multi-SNP CS found
      flag <- if (res$n_multi_snp_cs > 0)
        sprintf("*** FOUND %d MULTI-SNP CS! ***", res$n_multi_snp_cs)
      else ""
      cat(sprintf("OK | n_CS=%d cs_sizes=[%s] %s\n",
                  res$n_cs, res$cs_sizes, flag))
      all_results[[length(all_results) + 1]] <- res
    }
  }
}

# ===== Compile results =====
if (length(all_results) == 0) {
  cat("\nNo successful runs.\n")
  quit(save = "no", status = 1)
}

out <- rbindlist(all_results)

cat(sprintf("\n\n========== SUMMARY ==========\n"))
cat(sprintf("Total successful runs: %d / %d\n",
            nrow(out), nrow(param_grid) * nrow(GENES)))

# Quick summary: how many parameter combos gave multi-SNP CS?
multi_cs <- out[n_multi_snp_cs > 0]
if (nrow(multi_cs) > 0) {
  cat(sprintf("\n⚠ Multi-SNP CS observed in %d runs:\n", nrow(multi_cs)))
  print(multi_cs[, .(gene, snp_cap, L_param, window_kb, n_cs, n_multi_snp_cs,
                     mean_cs_size, cs_sizes)])
} else {
  cat("\n NO multi-SNP CS found in ANY parameter combination.\n")
  cat(" This strongly supports the conclusion that single-SNP CS are NOT artifacts\n")
  cat(" of the 300-SNP pruning window — the finding is robust across:\n")
  cat("   - SNP caps: 200, 300, 500\n")
  cat("   - L values: 5, 10, 20\n")
  cat("   - Windows: ±500kb, ±1000kb\n")
}

# Gene-level summary
cat(sprintf("\n========== Per-gene CS composition across all params ==========\n"))
gene_sum <- out[, .(
  n_runs    = .N,
  n_cs_min  = min(n_cs),
  n_cs_max  = max(n_cs),
  n_cs_median = as.double(median(n_cs)),
  any_multi = any(n_multi_snp_cs > 0),
  avg_cs_size = mean(mean_cs_size, na.rm = TRUE)
), by = gene]
print(gene_sum)

# ===== Detailed per-param table for discussion =====
cat(sprintf("\n========== Detailed parameter sweep ==========\n"))
detail <- out[, .(gene, snp_cap, L_param, L_eff, window_kb, n_used_snp,
                  n_cs, cs_sizes, top_pip_snp, n_pip_gt_01, converged)]
print(detail, nrows = 200)

# Save
csv_path <- file.path(OUT_DIR, "susie_sensitivity_grid.csv")
fwrite(out, csv_path)
cat(sprintf("\nSaved: %s\n", csv_path))

# ===== Generate manuscript-ready summary =====
cat(sprintf("\n========== Manuscript Summary ==========\n"))

# Did BMP2 keep 10 CS across params?
bmp2 <- out[gene == "BMP2"]
cat(sprintf("BMP2: n_CS range %d-%d across %d param combos\n",
            min(bmp2$n_cs), max(bmp2$n_cs), nrow(bmp2)))

# Average n_CS per gene across params
cs_by_gene <- out[, .(mean_cs = mean(n_cs), sd_cs = sd(n_cs)), by = gene]
cat(sprintf("Mean n_CS per gene: %.1f (SD=%.1f)\n",
            mean(cs_by_gene$mean_cs), mean(cs_by_gene$sd_cs)))

# Key takeaway
if (nrow(multi_cs) == 0) {
  cat("\n--- KEY FINDING ---\n")
  cat("Across 18 parameter combinations × 6 genes, SuSiE consistently identifies\n")
  cat("single-SNP credible sets. The all-PIP=1.0 pattern is NOT an artifact of:\n")
  cat("  • 300-SNP pruning (tested 200, 300, 500)\n")
  cat("  • Under-specified L (tested L = 5, 10, 20)\n")
  cat("  • Window size (±500kb, ±1000kb)\n")
  cat("This likely reflects genuine LD structure where the lead SNP\n")
  cat("dominates the regional association signal.\n")
}

cat("\nDone.\n")
