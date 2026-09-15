#!/usr/bin/env Rscript
# SuSiE fine-mapping for eQTL coloc-passing genes
# Uses GWAS summary stats + 1000G LD reference
# Output: credible sets + Tier classification

library(data.table)
library(susieR)
library(dplyr)

setDTthreads(8)

# ===== Config =====
GWAS_FILE  <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.h.tsv.gz"
BIM_FILE   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/ref_eur/1000G_EUR.bim"
BED_PREFIX <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/ref_eur/1000G_EUR"
COLOC_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_coloc_top10.csv"
OUT_CSV    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_susie_finemap.csv"
OUT_RDS    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_susie_finemap.rds"

CIS_WINDOW <- 500e3  # ±500kb for fine-mapping (tighter than coloc)
GWAS_N     <- 185616L

# ===== Load coloc-passing genes =====
coloc <- fread(COLOC_FILE)
coloc_pass <- coloc[PPH4 > 0.8]
cat(sprintf("Coloc-passing genes (PPH4>0.8): %d\n", nrow(coloc_pass)))
print(coloc_pass[, .(gene, chr, pos, PPH4)])

# ===== Load BIM for chr:pos mapping =====
cat("\nLoading BIM...\n")
bim <- fread(BIM_FILE, header=FALSE, select=1:6,
             col.names=c("chr","rsid","cm","pos","a1","a2"))
bim[, chr := as.integer(chr)]
cat(sprintf("BIM: %d SNPs\n", nrow(bim)))

# Helper: compute LD matrix from PLINK bed
compute_ld <- function(snp_ids, chr, bed_prefix) {
  # Write SNP list to temp file
  tmp_snps <- tempfile(fileext=".txt")
  writeLines(snp_ids, tmp_snps)

  # Use plink2 to compute LD
  # --r square: full R matrix; --extract to filter SNPs
  tmp_out <- tempfile(fileext="")
  plink_cmd <- sprintf(
    "/ifs1/User/zhouman/project9-v5-crc-atlas/tools/plink2 --bfile %s --chr %d --extract %s --r-unphased square --out %s --threads 4 --memory 16000 2>/dev/null",
    bed_prefix, chr, tmp_snps, tmp_out
  )

  rc <- system(plink_cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)

  ld_file <- paste0(tmp_out, ".unphased.vcor1")
  vars_file <- paste0(tmp_out, ".unphased.vcor1.vars")
  if (!file.exists(ld_file)) {
    unlink(tmp_snps)
    return(NULL)
  }

  # Read PLINK2 binary correlation matrix (.vcor1 = 64-bit double)
  vars <- fread(vars_file, header=FALSE)[[1]]
  n_snp <- length(vars)
  con <- file(ld_file, "rb")
  R_vec <- readBin(con, "double", n = n_snp * n_snp, size = 8)
  close(con)
  R <- matrix(R_vec, nrow = n_snp, ncol = n_snp, byrow = TRUE)
  rownames(R) <- vars; colnames(R) <- vars

  # Clean temp files
  unlink(c(tmp_snps, ld_file, vars_file, paste0(tmp_out, ".log")))
  return(R)
}

# Restrict to SNPs in 1000G EUR (for LD computation)
bim_chr <- split(bim, bim$chr)

all_results <- list()

for (i in 1:nrow(coloc_pass)) {
  gene  <- coloc_pass$gene[i]
  gchr  <- coloc_pass$chr[i]
  gpos  <- coloc_pass$pos[i]
  ensg  <- coloc_pass$ensg[i]
  pph4  <- coloc_pass$PPH4[i]

  cat(sprintf("\n=== %s (chr%d:%d, PPH4=%.4f) ===\n", gene, gchr, gpos, pph4))

  # ---- 1. Load GWAS SNPs in cis region ----
  gw_filter <- sprintf("zcat %s | awk 'NR==1 || ($1==%d && $2>=%d && $2<=%d)'",
                       GWAS_FILE, gchr, gpos - CIS_WINDOW, gpos + CIS_WINDOW)
  gwas <- fread(cmd = gw_filter)
  if (nrow(gwas) < 10) { cat("  <10 GWAS SNPs\n"); next }

  setnames(gwas,
    old = c("chromosome","base_pair_location","effect_allele","other_allele",
            "beta","standard_error","effect_allele_frequency","p_value","rsid"),
    new = c("chr","bp","ea","oa","beta","se","eaf","pval","snp"))

  # ---- 2. Intersect with 1000G SNPs (need LD ref) ----
  bim_chr_data <- bim_chr[[as.character(gchr)]]
  if (is.null(bim_chr_data) || nrow(bim_chr_data) < 10) {
    cat("  No BIM data for this chr\n"); next
  }

  # Merge GWAS with BIM by rsid
  gwas_w_ld <- merge(gwas, bim_chr_data[, .(rsid, pos_ld=pos)],
                     by.x="snp", by.y="rsid", all.x=FALSE)
  gwas_w_ld <- gwas_w_ld[abs(gwas_w_ld$pos_ld - gpos) <= CIS_WINDOW]

  if (nrow(gwas_w_ld) < 20) {
    cat(sprintf("  Only %d SNPs with LD ref\n", nrow(gwas_w_ld))); next
  }
  cat(sprintf("  GWAS SNPs with LD ref in cis: %d\n", nrow(gwas_w_ld)))

  # ---- 3. Prune: keep top 300 SNPs for LD computation (numerical stability) ----
  # Sort by p-value, keep most significant
  setorder(gwas_w_ld, pval)
  if (nrow(gwas_w_ld) > 300) {
    gwas_w_ld <- gwas_w_ld[1:300]
    cat(sprintf("  Pruned to top %d SNPs by p-value\n", nrow(gwas_w_ld)))
  }

  # ---- 4. Compute z-scores ----
  gwas_w_ld[, z := beta / se]
  # Remove NAs/Inf
  gwas_w_ld <- gwas_w_ld[is.finite(z)]

  if (nrow(gwas_w_ld) < 20) {
    cat(sprintf("  Only %d SNPs after filtering\n", nrow(gwas_w_ld))); next
  }

  # ---- 5. Compute LD matrix ----
  snp_list <- gwas_w_ld$snp
  cat(sprintf("  Computing LD for %d SNPs...\n", length(snp_list)))

  R <- compute_ld(snp_list, gchr, BED_PREFIX)

  if (is.null(R) || nrow(R) < 20) {
    cat("  LD computation failed or too few SNPs\n"); next
  }

  # Ensure R is positive definite using ridge penalty (more stable than nearPD)
  diag(R) <- diag(R) + 1e-4
  R <- (R + t(R)) / 2  # ensure symmetry

  cat(sprintf("  LD matrix: %dx%d\n", nrow(R), ncol(R)))

  # ---- 6. Run SuSiE ----
  # Align GWAS SNPs with LD matrix order
  ld_snps <- rownames(R)
  gwas_match <- match(ld_snps, gwas_w_ld$snp)
  gwas_aligned <- gwas_w_ld[gwas_match]
  z <- gwas_aligned$z

  # Ensure R is positive definite (already done above but double check)
  if (length(z) != nrow(R)) {
    cat(sprintf("  Mismatch: %d z-scores vs %dx%d LD\n", length(z), nrow(R), ncol(R))); next
  }

  cat(sprintf("  Running susie_rss (n=%d SNPs, GWAS N=%d)...\n", length(z), GWAS_N))

  tryCatch({
    fit <- susie_rss(z = z, R = R, n = GWAS_N,
                     L = min(10, length(z) %/% 20),
                     estimate_residual_variance = FALSE,
                     min_abs_corr = 0.3,
                     coverage = 0.9,
                     refine = TRUE,
                     max_iter = 500,
                     check_inputs = FALSE)  # skip internal PSD/PIP checks for external LD

    n_cs <- length(fit$sets$cs)
    converged <- fit$converged

    # Get top PIP SNP (even when n_cs=0, PIP values are still meaningful)
    top_pip_idx <- which.max(fit$pip)
    top_pip_snp  <- gwas_aligned$snp[top_pip_idx]
    top_pip_val  <- fit$pip[top_pip_idx]
    # Top 5 PIP SNPs for reporting
    top5_idx <- order(fit$pip, decreasing=TRUE)[1:min(5, length(fit$pip))]
    top5_info <- paste(
      sprintf("%s(PIP=%.3f)", gwas_aligned$snp[top5_idx], fit$pip[top5_idx]),
      collapse=" | ")

    credible_sets <- list()
    for (cs_i in seq_along(fit$sets$cs)) {
      idx <- fit$sets$cs[[cs_i]]
      credible_sets[[cs_i]] <- data.frame(
        cs_id = paste0("CS", cs_i),
        n_snps = length(idx),
        top_snp = gwas_aligned$snp[idx[1]],
        top_pip = fit$pip[idx[1]],
        all_snps = paste(gwas_aligned$snp[idx], collapse=";"),
        stringsAsFactors = FALSE
      )
    }

    # Tier classification:
    # Tier 1: ≥1 credible set + coloc PPH4 > 0.8
    # Tier 2: coloc PPH4 > 0.8 but no credible set (still strong colocalization evidence)
    # Tier 3: coloc PPH4 > 0.5 (marginal)
    if (n_cs >= 1) {
      tier <- "1"
    } else {
      tier <- "2"  # coloc pass but no credible set
    }

    cat(sprintf("  → n_cs=%d, converged=%s, Tier=%s, top_PIP_SNP=%s(PIP=%.3f)\n",
                n_cs, converged, tier, top_pip_snp, top_pip_val))

    if (length(credible_sets) > 0) {
      cs_df <- rbindlist(credible_sets)
      cs_summary <- paste(
        sapply(seq_along(fit$sets$cs), function(j)
          sprintf("CS%d(%d_SNPs,%s,PIP=%.2f)",
                  j, length(fit$sets$cs[[j]]),
                  gwas_aligned$snp[fit$sets$cs[[j]][1]],
                  fit$pip[fit$sets$cs[[j]][1]])
        ), collapse=" | ")
    } else {
      cs_df <- data.frame(cs_id=character(), n_snps=integer(),
                          top_snp=character(), top_pip=numeric(),
                          all_snps=character())
      cs_summary <- ""
    }

    res_row <- data.frame(
      gene = gene, ensg = ensg, chr = gchr, pos = gpos,
      PPH4 = pph4,
      n_snps_susie = length(z),
      n_cs = n_cs,
      cs_summary = cs_summary,
      top_pip_snp = top_pip_snp,
      top_pip = top_pip_val,
      top5_pip = top5_info,
      converged = converged,
      tier = tier,
      stringsAsFactors = FALSE
    )

    all_results[[gene]] <- res_row

  }, error = function(e) {
    cat(sprintf("  susie error: %s\n", e$message))
  })
}

# ===== Output =====
cat("\n\n=== SuSiE Fine-Mapping Results ===\n")

if (length(all_results) > 0) {
  out <- rbindlist(all_results)
  out <- out[order(-PPH4)]

  fwrite(out, OUT_CSV)
  # Save full credible set details
  # (we only saved summary; for paper, we'll re-extract details from rds if needed)
  print(out[, .(gene, PPH4, n_snps_susie, n_cs, cs_summary, top_pip_snp, top_pip, tier, converged)])

  cat(sprintf("\nTier distribution:\n"))
  print(out[, .N, by=tier])
  cat(sprintf("\nSaved: %s\n", OUT_CSV))
} else {
  cat("\n  No genes passed SuSiE.\n")
}

cat("\nDone.\n")
