#!/usr/bin/env Rscript
# Phase 2: UKB-PPP pQTL MR for 5 genes (CDKN1A, LTA, NCF2, TET2, TNF)
# UKB-PPP = hg38, GWAS = hg37 → match by chr + nearest position + allele harmonization

library(data.table)

# ===== Config =====
UKBPP_DIR <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/ukbppp_merged"
GWAS_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.h.tsv.gz"
OUT_DIR   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase2_pqtl"
CIS_WINDOW <- 1000000L   # ±1Mb
PQTL_THRESH <- 5e-8
MATCH_DIST <- 500L  # bp tolerance for hg37↔hg38 matching

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Gene info (hg38 TSS ± broad region from GENCODE v47)
genes <- list(
  CDKN1A = list(chr = 6,  pos = 36652392, file = "CDKN1A_merged.txt.gz", ensg = "ENSG00000124762"),
  LTA    = list(chr = 6,  pos = 31572054, file = "LTA_merged.txt.gz",    ensg = "ENSG00000226979"),
  NCF2   = list(chr = 1,  pos = 183555568,file = "NCF2_merged.txt.gz",   ensg = "ENSG00000116701"),
  TET2   = list(chr = 4,  pos = 105145428,file = "TET2_merged.txt.gz",   ensg = "ENSG00000168769"),
  TNF    = list(chr = 6,  pos = 31575565, file = "TNF_merged.txt.gz",    ensg = "ENSG00000232810")
)

# ===== Load GWAS =====
cat("Loading CRC GWAS...\n")
gwas <- fread(GWAS_FILE)
setnames(gwas,
  old = c("chromosome","base_pair_location","effect_allele","other_allele",
          "beta","standard_error","effect_allele_frequency","p_value","rsid"),
  new = c("chr","bp","ea","oa","beta","se","eaf","pval","snp"))
gwas <- gwas[, .(chr, bp, ea, oa, beta, se, eaf, pval, snp)]
setkey(gwas, chr, bp)
cat(sprintf("  GWAS SNPs: %d\n", nrow(gwas)))

# ===== Helper: complement =====
complement_allele <- function(x) chartr("ACGTacgt", "TGCAtgca", x)

# ===== Process each gene =====
all_results <- list()

for (gene_name in names(genes)) {
  cat(sprintf("\n========== %s ==========\n", gene_name))
  g <- genes[[gene_name]]

  # 1. Load UKB-PPP (pre-filter to cis chr:start-end to avoid loading genome-wide file)
  pqtl_file <- file.path(UKBPP_DIR, g$file)
  if (!file.exists(pqtl_file)) { cat(sprintf("SKIP: no file %s\n", pqtl_file)); next }
  cat(sprintf("Loading %s (cis chr%d:%d-%d)...\n", basename(pqtl_file),
              g$chr, g$pos - CIS_WINDOW, g$pos + CIS_WINDOW))

  # awk pre-filter: print header (NR==1) + rows where chr==g$chr and bp in cis window
  filter_cmd <- sprintf("zcat %s | awk 'NR==1 || ($1==%d && $2>=%d && $2<=%d)'",
                        pqtl_file, g$chr, g$pos - CIS_WINDOW, g$pos + CIS_WINDOW)
  pqtl <- fread(cmd = filter_cmd)

  if (nrow(pqtl) == 0) { cat("No SNPs in cis window\n"); next }

  # Columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
  setnames(pqtl,
    old = c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","BETA","SE","LOG10P","N"),
    new = c("chr","bp","snp_id","oa","ea","eaf","beta","se","log10p","n"))

  # 2. P-value threshold
  pqtl[, pval := 10^(-log10p)]
  total_cis <- nrow(pqtl)
  pqtl <- pqtl[pval < PQTL_THRESH]
  cat(sprintf("Cis SNPs: %d total, %d significant (p<5e-8)\n", total_cis, nrow(pqtl)))
  if (nrow(pqtl) == 0) next

  # 3. Match to GWAS by nearest position
  gw_chr <- gwas[chr == g$chr]
  setkey(gw_chr, bp)
  setkey(pqtl, bp)

  matched <- gw_chr[pqtl, on = .(chr, bp), roll = "nearest",
    .(chr, bp.pqtl = i.bp, bp.gwas = bp,
      pqtl.snp = i.snp_id, gwas.snp = snp,
      pqtl.ea = i.ea, pqtl.oa = i.oa, pqtl.beta = i.beta, pqtl.se = i.se,
      pqtl.pval = i.pval, pqtl.eaf = i.eaf, pqtl.n = i.n,
      gwas.ea = ea, gwas.oa = oa, gwas.beta = beta, gwas.se = se,
      gwas.pval = pval, gwas.eaf = eaf)]

  matched[, dist := abs(bp.pqtl - bp.gwas)]
  matched <- matched[dist <= MATCH_DIST & !is.na(gwas.snp)]
  cat(sprintf("Matched to GWAS (<=%dbp): %d\n", MATCH_DIST, nrow(matched)))
  if (nrow(matched) == 0) next

  # 4. Allele harmonization (vectorized)
  p_ea <- matched$pqtl.ea; p_oa <- matched$pqtl.oa
  g_ea <- matched$gwas.ea; g_oa <- matched$gwas.oa
  g_beta <- matched$gwas.beta

  direct  <- (p_ea == g_ea & p_oa == g_oa)
  flip    <- (p_ea == g_oa & p_oa == g_ea)
  swgdam  <- (p_ea == complement_allele(g_ea) & p_oa == complement_allele(g_oa))
  swgdam_flip <- (p_ea == complement_allele(g_oa) & p_oa == complement_allele(g_ea))

  # Flip GWAS beta for strand-flip cases
  flip_idx <- which(flip | swgdam_flip)
  g_beta[flip_idx] <- -g_beta[flip_idx]

  keep <- direct | flip  # SWGDAM cases: alleles are same on opposite strand, beta doesn't flip
  matched <- matched[keep]
  matched$gwas.beta <- g_beta[keep]
  cat(sprintf("After harmonization: %d (direct=%d, flipped=%d)\n",
    nrow(matched), sum(direct[keep]), sum(flip[keep])))

  if (nrow(matched) == 0) { cat("No instruments after harmonization.\n"); next }

  # 5. Select top cis-pQTL per gene
  # cis-pQTLs for a single protein are typically in high LD → top SNP is standard (Wang/Tian 2025)
  setorder(matched, pqtl.pval)
  matched <- matched[!duplicated(gwas.snp)]

  # Also report total matched for reference
  cat(sprintf("Matched instruments (before LD pruning): %d\n", nrow(matched)))

  # For multi-instrument: also compute IVW using all matched (for supplementary)
  # But primary analysis uses top SNP Wald ratio
  top_snp <- matched[1]
  cat(sprintf("Top instrument: %s p=%.2e\n", top_snp$gwas.snp, top_snp$pqtl.pval))

  # 6. MR analysis
  # Primary: Wald ratio with top cis-pQTL
  wald_b <- top_snp$gwas.beta / top_snp$pqtl.beta
  wald_se <- abs(wald_b) * sqrt((top_snp$pqtl.se / top_snp$pqtl.beta)^2 +
                                 (top_snp$gwas.se / top_snp$gwas.beta)^2)
  wald_p <- 2 * pnorm(-abs(wald_b / wald_se))

  mr_result <- data.table(
    gene = gene_name, ensg = g$ensg, chr = g$chr,
    n_instruments = nrow(matched),
    top_snp = top_snp$gwas.snp,
    top_pqtl_p = top_snp$pqtl.pval,
    method = "Wald ratio",
    b = wald_b, se = wald_se, pval = wald_p
  )

  # Supplementary: IVW if multiple instruments
  if (nrow(matched) > 1) {
    ratio <- matched$gwas.beta / matched$pqtl.beta
    ratio_se <- abs(ratio) * sqrt((matched$pqtl.se / matched$pqtl.beta)^2 +
                                   (matched$gwas.se / matched$gwas.beta)^2)
    w <- 1 / ratio_se^2
    ivw_b <- sum(w * ratio) / sum(w)
    ivw_se <- 1 / sqrt(sum(w))
    ivw_p <- 2 * pnorm(-abs(ivw_b / ivw_se))

    q_stat <- sum(w * (ratio - ivw_b)^2)
    q_pval <- pchisq(q_stat, df = nrow(matched) - 1, lower.tail = FALSE)

    mr_result[, `:=`(ivw_b = ivw_b, ivw_se = ivw_se, ivw_pval = ivw_p,
                     cochran_q = q_stat, cochran_q_pval = q_pval)]
  }

  cat(sprintf("=> %s: b=%.4f se=%.4f p=%.3e (%s, %d instruments)\n",
    gene_name, mr_result$b, mr_result$se, mr_result$pval, mr_result$method, nrow(matched)))
  all_results[[gene_name]] <- mr_result

  # Save per-gene
  fwrite(matched[, .(pqtl.snp, gwas.snp, chr, bp.pqtl, bp.gwas,
    pqtl.ea, pqtl.oa, pqtl.beta, pqtl.se, pqtl.pval,
    gwas.ea, gwas.oa, gwas.beta, gwas.se, gwas.pval)],
    file.path(OUT_DIR, sprintf("%s_instruments.csv", gene_name)))
  fwrite(mr_result, file.path(OUT_DIR, sprintf("%s_mr_results.csv", gene_name)))
}

# ===== Summary =====
cat("\n\n========== SUMMARY ==========\n")
combined <- rbindlist(all_results)
print(combined, digits = 4)
fwrite(combined, file.path(OUT_DIR, "phase2_pqtl_mr_combined.csv"))
cat(sprintf("\nResults saved to %s\n", file.path(OUT_DIR, "phase2_pqtl_mr_combined.csv")))
