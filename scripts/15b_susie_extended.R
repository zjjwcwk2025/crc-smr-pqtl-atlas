#!/usr/bin/env Rscript
# Phase 3c: SuSiE fine-mapping — EXTENDED (21 new genes: 6 MHC + 15 non-MHC)
# Goal: (1) Evaluate MHC coloc PPH4 inflation via SuSiE LD-aware fine-mapping
#       (2) Expand SuSiE coverage from 6 to 27 coloc-passing genes
# Based on scripts/15_susie_finemap.R — same methodology
# Output: results/phase3c_susie_extended.csv

library(data.table)
library(susieR)
library(dplyr)

setDTthreads(8)

# ===== Config =====
GWAS_FILE  <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.h.tsv.gz"
BIM_FILE   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/ref_eur/1000G_EUR.bim"
BED_PREFIX <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/ref_eur/1000G_EUR"
COLOC_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_coloc_all75.csv"
OUT_CSV    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3c_susie_extended.csv"
OUT_RDS    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3c_susie_extended.rds"

CIS_WINDOW <- 500e3
GWAS_N     <- 185616L

# MHC region (hg38/GRCh38)
MHC_CHR_MIN <- 29e6
MHC_CHR_MAX <- 33e6

# ===== Select target genes: MHC PASS + Top non-MHC PASS (exclude already-done 6) =====
already_done <- c("UTP11","BMP2","PNKD","LAMC1","TMBIM1","GPBAR1")

coloc_all <- fread(COLOC_FILE)
coloc_all <- coloc_all[PPH4 > 0.8]  # only PASS genes

# MHC: chr6 29-33Mb (hg38/GRCh38)
mhc_genes <- coloc_all[chr == 6 & pos >= MHC_CHR_MIN & pos <= MHC_CHR_MAX & !gene %in% already_done]

# Non-MHC: top PPH4, protein-coding only (exclude ENSG-only entries)
non_mhc <- coloc_all[!(chr == 6 & pos >= MHC_CHR_MIN & pos <= MHC_CHR_MAX) & !gene %in% already_done]
non_mhc <- non_mhc[!grepl("^ENSG", gene)]  # exclude non-coding
non_mhc <- non_mhc[order(-PPH4)][1:min(15, .N)]

target_genes <- rbind(mhc_genes, non_mhc)

cat(sprintf("Target genes for extended SuSiE: %d\n", nrow(target_genes)))
cat(sprintf("  MHC: %d genes\n", nrow(mhc_genes)))
cat(sprintf("  Non-MHC: %d genes\n", nrow(non_mhc)))
cat("\nTarget genes:\n")
for (i in 1:nrow(target_genes)) {
  cat(sprintf("  %d. %s (chr%d:%d, PPH4=%.4f)\n",
              i, target_genes$gene[i], target_genes$chr[i],
              target_genes$pos[i], target_genes$PPH4[i]))
}

# ===== Load BIM for chr:pos mapping =====
cat("\nLoading BIM...\n")
bim <- fread(BIM_FILE, header=FALSE, select=1:6,
             col.names=c("chr","rsid","cm","pos","a1","a2"))
bim[, chr := as.integer(chr)]
bim_chr <- split(bim, bim$chr)
cat(sprintf("BIM: %d SNPs across %d chromosomes\n", nrow(bim), length(bim_chr)))

# ===== Helper: compute LD matrix from PLINK bed =====
compute_ld <- function(snp_ids, chr, bed_prefix) {
  tmp_snps <- tempfile(fileext=".txt")
  writeLines(snp_ids, tmp_snps)
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
  vars <- fread(vars_file, header=FALSE)[[1]]
  n_snp <- length(vars)
  con <- file(ld_file, "rb")
  R_vec <- readBin(con, "double", n = n_snp * n_snp, size = 8)
  close(con)
  R <- matrix(R_vec, nrow = n_snp, ncol = n_snp, byrow = TRUE)
  rownames(R) <- vars; colnames(R) <- vars
  unlink(c(tmp_snps, ld_file, vars_file, paste0(tmp_out, ".log")))
  return(R)
}

# ===== Run SuSiE per gene =====
all_results <- list()

for (i in 1:nrow(target_genes)) {
  gene  <- target_genes$gene[i]
  gchr  <- target_genes$chr[i]
  gpos  <- target_genes$pos[i]
  ensg  <- target_genes$ensg[i]
  pph4  <- target_genes$PPH4[i]
  is_mhc <- (gchr == 6 && gpos >= MHC_CHR_MIN && gpos <= MHC_CHR_MAX)

  cat(sprintf("\n=== [%d/%d] %s (chr%d:%d, PPH4=%.4f, MHC=%s) ===\n",
              i, nrow(target_genes), gene, gchr, gpos, pph4, is_mhc))

  # ---- 1. Load GWAS SNPs in cis region ----
  gw_filter <- sprintf("zcat %s | awk 'NR==1 || ($1==%d && $2>=%d && $2<=%d)'",
                       GWAS_FILE, gchr, gpos - CIS_WINDOW, gpos + CIS_WINDOW)
  gwas <- fread(cmd = gw_filter)
  if (nrow(gwas) < 10) { cat("  <10 GWAS SNPs, skip\n"); next }

  setnames(gwas,
    old = c("chromosome","base_pair_location","effect_allele","other_allele",
            "beta","standard_error","effect_allele_frequency","p_value","rsid"),
    new = c("chr","bp","ea","oa","beta","se","eaf","pval","snp"))

  # ---- 2. Intersect with 1000G SNPs ----
  bim_chr_data <- bim_chr[[as.character(gchr)]]
  if (is.null(bim_chr_data) || nrow(bim_chr_data) < 10) {
    cat("  No BIM data for this chr\n"); next
  }

  gwas_w_ld <- merge(gwas, bim_chr_data[, .(rsid, pos_ld=pos)],
                     by.x="snp", by.y="rsid", all.x=FALSE)
  gwas_w_ld <- gwas_w_ld[abs(gwas_w_ld$pos_ld - gpos) <= CIS_WINDOW]

  if (nrow(gwas_w_ld) < 20) {
    cat(sprintf("  Only %d SNPs with LD ref\n", nrow(gwas_w_ld))); next
  }
  cat(sprintf("  GWAS SNPs in cis with LD ref: %d\n", nrow(gwas_w_ld)))

  # ---- 3. Prune to top 300 SNPs by p-value ----
  setorder(gwas_w_ld, pval)
  if (nrow(gwas_w_ld) > 300) {
    gwas_w_ld <- gwas_w_ld[1:300]
    cat(sprintf("  Pruned to top %d SNPs\n", nrow(gwas_w_ld)))
  }

  # ---- 4. Compute z-scores ----
  gwas_w_ld[, z := beta / se]
  gwas_w_ld <- gwas_w_ld[is.finite(z)]
  if (nrow(gwas_w_ld) < 20) { next }

  # ---- 5. Compute LD matrix ----
  snp_list <- gwas_w_ld$snp
  cat(sprintf("  Computing LD for %d SNPs...\n", length(snp_list)))
  R <- compute_ld(snp_list, gchr, BED_PREFIX)
  if (is.null(R) || nrow(R) < 20) { cat("  LD failed\n"); next }

  diag(R) <- diag(R) + 1e-4
  R <- (R + t(R)) / 2
  cat(sprintf("  LD matrix: %dx%d\n", nrow(R), ncol(R)))

  # ---- 6. Run SuSiE ----
  ld_snps <- rownames(R)
  gwas_match <- match(ld_snps, gwas_w_ld$snp)
  gwas_aligned <- gwas_w_ld[gwas_match]
  z <- gwas_aligned$z

  if (length(z) != nrow(R)) {
    cat(sprintf("  Mismatch: %d z vs %dx%d LD\n", length(z), nrow(R), ncol(R))); next
  }

  # For MHC: try L=10 and L=20 (MHC has more complex LD)
  L_values <- if (is_mhc) c(10, 20) else c(10)
  best_fit <- NULL

  for (L_val in L_values) {
    cat(sprintf("  Trying L=%d (n=%d SNPs)...\n", L_val, length(z)))
    tryCatch({
      fit <- susie_rss(z = z, R = R, n = GWAS_N,
                       L = L_val,
                       estimate_residual_variance = FALSE,
                       min_abs_corr = 0.3,
                       coverage = 0.9,
                       refine = TRUE,
                       max_iter = 500,
                       check_inputs = FALSE)

      if (is.null(best_fit) || fit$converged) {
        best_fit <- fit
        cat(sprintf("    converged=%s, n_cs=%d\n", fit$converged, length(fit$sets$cs)))
      }
    }, error = function(e) {
      cat(sprintf("    L=%d error: %s\n", L_val, e$message))
    })
  }

  if (is.null(best_fit)) { cat("  All L values failed\n"); next }

  fit <- best_fit
  n_cs <- length(fit$sets$cs)
  converged <- fit$converged

  # Top PIP
  top_pip_idx <- which.max(fit$pip)
  top_pip_snp  <- gwas_aligned$snp[top_pip_idx]
  top_pip_val  <- fit$pip[top_pip_idx]
  top5_idx <- order(fit$pip, decreasing=TRUE)[1:min(5, length(fit$pip))]
  top5_info <- paste(
    sprintf("%s(PIP=%.3f)", gwas_aligned$snp[top5_idx], fit$pip[top5_idx]),
    collapse=" | ")

  # CS summary
  if (n_cs > 0) {
    cs_summary <- paste(
      sapply(seq_along(fit$sets$cs), function(j)
        sprintf("CS%d(%d_SNPs,%s,PIP=%.2f)",
                j, length(fit$sets$cs[[j]]),
                gwas_aligned$snp[fit$sets$cs[[j]][1]],
                fit$pip[fit$sets$cs[[j]][1]])
      ), collapse=" | ")
  } else {
    cs_summary <- "no_credible_set"
  }

  # Multi-SNP CS check (for MHC evaluation)
  n_multi_snp_cs <- if (n_cs > 0) {
    sum(sapply(fit$sets$cs, length) > 1)
  } else 0

  # Tier (same as original)
  tier <- if (n_cs >= 1) "1" else "2"

  cat(sprintf("  → n_cs=%d, n_multi_SNP_CS=%d, converged=%s, Tier=%s, top_PIP=%s(%.3f)\n",
              n_cs, n_multi_snp_cs, converged, tier, top_pip_snp, top_pip_val))

  res_row <- data.frame(
    gene = gene, ensg = ensg, chr = gchr, pos = gpos,
    PPH4 = pph4, is_mhc = is_mhc,
    n_snps_susie = length(z),
    n_cs = n_cs, n_multi_snp_cs = n_multi_snp_cs,
    cs_summary = cs_summary,
    top_pip_snp = top_pip_snp,
    top_pip = top_pip_val,
    top5_pip = top5_info,
    converged = converged,
    tier = tier,
    stringsAsFactors = FALSE
  )

  all_results[[gene]] <- res_row
}

# ===== Output =====
cat("\n\n=== SuSiE Extended Results ===\n")

if (length(all_results) > 0) {
  out <- rbindlist(all_results)
  out <- out[order(-PPH4)]

  # Split into MHC vs non-MHC for reporting
  out_mhc <- out[is_mhc == TRUE]
  out_clean <- out[is_mhc == FALSE]

  cat(sprintf("\n--- MHC genes (%d) ---\n", nrow(out_mhc)))
  if (nrow(out_mhc) > 0) {
    print(out_mhc[, .(gene, PPH4, n_cs, n_multi_snp_cs, cs_summary, top_pip_snp, converged, tier)])
  }

  cat(sprintf("\n--- Non-MHC genes (%d) ---\n", nrow(out_clean)))
  if (nrow(out_clean) > 0) {
    print(out_clean[, .(gene, PPH4, n_cs, cs_summary, top_pip_snp, converged, tier)])
  }

  # MHC assessment summary
  if (nrow(out_mhc) > 0) {
    cat("\n--- MHC Assessment ---\n")
    # Good signs: converged + single-SNP CS or no CS
    # Bad signs: multi-SNP CS or non-converged → LD inflation likely
    n_mhc_converged <- sum(out_mhc$converged, na.rm=TRUE)
    n_mhc_multi_cs <- sum(out_mhc$n_multi_snp_cs > 0, na.rm=TRUE)
    cat(sprintf("  MHC genes with convergent SuSiE: %d/%d\n", n_mhc_converged, nrow(out_mhc)))
    cat(sprintf("  MHC genes with multi-SNP CS: %d (possible LD complexity)\n", n_mhc_multi_cs))
    if (n_mhc_multi_cs > 0) {
      cat("  ⚠️ Multi-SNP CS in MHC suggests complex LD → coloc PPH4 may be inflated\n")
    }
    if (n_mhc_converged == nrow(out_mhc) && n_mhc_multi_cs == 0) {
      cat("  ✅ All MHC genes converged with single-SNP CS → coloc PPH4 supported by SuSiE\n")
    }
  }

  # Overall tier distribution
  cat(sprintf("\nTier distribution (extended):\n"))
  print(out[, .N, by=tier])

  # Combine with original SuSiE results
  orig_file <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_susie_finemap.csv"
  if (file.exists(orig_file)) {
    orig <- fread(orig_file)
    orig[, is_mhc := (chr == 6 & pos >= MHC_CHR_MIN & pos <= MHC_CHR_MAX)]
    combined <- rbind(orig, out, fill=TRUE)
    combined <- combined[order(-PPH4)]
    fwrite(combined, "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_susie_combined.csv")
    cat(sprintf("Combined SuSiE results: %d genes → results/phase3_susie_combined.csv\n", nrow(combined)))
    cat(sprintf("  Tier 1: %d genes\n", sum(combined$tier == "1")))
    cat(sprintf("  Tier 2: %d genes\n", sum(combined$tier == "2")))
  }

  fwrite(out, OUT_CSV)
  cat(sprintf("\nSaved: %s\n", OUT_CSV))
} else {
  cat("\n⚠ No genes passed SuSiE.\n")
}

cat("\nDone.\n")
