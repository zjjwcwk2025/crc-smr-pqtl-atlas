#!/usr/bin/env Rscript
# Phase 2: deCODE pQTL MR for 12 unique genes (not covered by UKB-PPP)
# deCODE & GWAS both hg38 → direct position matching
library(data.table)

# ===== Config =====
DECODE_DIR <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/decode"
GWAS_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.h.tsv.gz"
OUT_DIR   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase2_decode"
CIS_WINDOW <- 1000000L
PQTL_THRESH <- 5e-8
MATCH_DIST <- 500L

dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# ===== Gene map (12 deCODE-only genes, hg38 TSS) =====
genes <- list(
  STAT6      = list(chr=12, pos=57491338,  pqtl_file="10372_18_STAT6_STAT6.txt.gz",      ensg="ENSG00000166888"),
  LIMA1      = list(chr=12, pos=50623450,  pqtl_file="11543_84_LIMA1_LIMA1.txt.gz",      ensg="ENSG00000050405"),
  CCM2       = list(chr=7,  pos=44976160,  pqtl_file="12347_29_CCM2_CCM2.txt.gz",        ensg="ENSG00000136280"),
  CERS5      = list(chr=12, pos=50522623,  pqtl_file="13494_6_CERS5_CERS5.txt.gz",       ensg="ENSG00000139624"),
  STAB1      = list(chr=3,  pos=52545033,  pqtl_file="14599_18_STAB1_STAB1.txt.gz",      ensg="ENSG00000010327"),
  MZF1       = list(chr=19, pos=59073298,  pqtl_file="14662_6_MZF1_MZF1.txt.gz",         ensg="ENSG00000099326"),
  BMP2       = list(chr=20, pos=6754619,   pqtl_file="15666_21_BMP2_BMP_2.txt.gz",       ensg="ENSG00000125845"),
  ARPC5      = list(chr=1,  pos=183620013, pqtl_file="18419_20_ARPC5_p16_ARC.txt.gz",    ensg="ENSG00000162704"),
  UBE2M      = list(chr=19, pos=59059813,  pqtl_file="19111_10_UBE2M_UBC12.txt.gz",      ensg="ENSG00000130725"),
  CABLES2    = list(chr=20, pos=60963688,  pqtl_file="7105_7_CABLES2_CABL2.txt.gz",       ensg="ENSG00000149679"),
  TNFRSF1B   = list(chr=1,  pos=12155565,  pqtl_file="3152_57_TNFRSF1B_TNF_sR_II.txt.gz", ensg="ENSG00000028137"),
  TNFRSF1A   = list(chr=12, pos=6291050,   pqtl_file="2654_19_TNFRSF1A_TNF_sR_I.txt.gz",  ensg="ENSG00000067182")
)
cat(sprintf("DeCODE-only genes: %d\n", length(genes)))

# ===== Load GWAS =====
cat("\nLoading GWAS (hg38, GCST90255675)...\n")
gwas <- fread(GWAS_FILE)
setnames(gwas,
  old=c("chromosome","base_pair_location","effect_allele","other_allele",
        "beta","standard_error","effect_allele_frequency","p_value","rsid"),
  new=c("chr","bp","ea","oa","beta","se","eaf","pval","snp"))
gwas <- gwas[, .(chr, bp, ea, oa, beta, se, eaf, pval, snp)]
setkey(gwas, chr, bp)
cat(sprintf("  GWAS SNPs: %d\n", nrow(gwas)))

# ===== Helpers =====
complement <- function(x) chartr("ACGTacgt", "TGCAtgca", x)

# ===== Process each gene =====
all_results <- list()

for (gname in names(genes)) {
  cat(sprintf("\n========== %s ==========\n", gname))
  g <- genes[[gname]]
  pqtl_file <- file.path(DECODE_DIR, g$pqtl_file)
  if (!file.exists(pqtl_file)) { cat(sprintf("SKIP: %s not found\n", pqtl_file)); next }

  # 1. Load deCODE pQTL (cis-region only via awk pre-filter)
  # deCODE columns: Chrom Pos Name rsids effectAllele otherAllele Beta Pval minus_log10_pval SE N ImpMAF
  # Chrom has "chr" prefix; strip it for matching
  cat(sprintf("Loading %s (cis chr%d:%d-%d)...\n", basename(pqtl_file),
              g$chr, g$pos-CIS_WINDOW, g$pos+CIS_WINDOW))
  filter_cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==1 || ($1==\"chr%d\" && $2>=%d && $2<=%d)'",
                        pqtl_file, g$chr, g$pos-CIS_WINDOW, g$pos+CIS_WINDOW)
  pqtl <- fread(cmd=filter_cmd)
  if (nrow(pqtl) <= 1) { cat("No SNPs in cis window\n"); next }

  # Rename columns
  setnames(pqtl, old=c("Chrom","Pos","Name","rsids","effectAllele","otherAllele",
                        "Beta","Pval","SE","N","ImpMAF"),
           new=c("chr_raw","bp","name","rsid","ea","oa","beta","pval","se","n","maf"))

  total_cis <- nrow(pqtl)
  cat(sprintf("Cis SNPs total: %d\n", total_cis))

  # 2. Filter significant cis-pQTLs
  pqtl[, chr := as.integer(gsub("chr","",chr_raw))]
  pqtl_sig <- pqtl[pval < PQTL_THRESH]
  cat(sprintf("Significant cis-pQTLs (p<5e-8): %d\n", nrow(pqtl_sig)))
  if (nrow(pqtl_sig) == 0) next

  # 3. Match to GWAS by chr + position (both hg38)
  setkey(pqtl_sig, chr, bp)
  gw_chr <- gwas[chr == g$chr]
  setkey(gw_chr, bp)

  matched <- gw_chr[pqtl_sig, on=.(chr, bp), roll="nearest",
    .(chr, bp.pqtl=i.bp, bp.gwas=bp,
      pqtl.rsid=i.rsid, gwas.rsid=snp,
      pqtl.ea=i.ea, pqtl.oa=i.oa, pqtl.beta=i.beta, pqtl.se=i.se, pqtl.pval=i.pval,
      gwas.ea=ea, gwas.oa=oa, gwas.beta=beta, gwas.se=se, gwas.pval=pval)]

  matched[, dist := abs(bp.pqtl - bp.gwas)]
  matched <- matched[dist <= MATCH_DIST & !is.na(gwas.rsid)]
  cat(sprintf("Matched to GWAS (<=%dbp): %d\n", MATCH_DIST, nrow(matched)))
  if (nrow(matched) == 0) next

  # 4. Allele harmonization
  p_ea <- matched$pqtl.ea; p_oa <- matched$pqtl.oa
  g_ea <- matched$gwas.ea; g_oa <- matched$gwas.oa
  g_beta <- matched$gwas.beta

  direct <- (p_ea == g_ea & p_oa == g_oa)
  flip   <- (p_ea == g_oa & p_oa == g_ea)
  # Strand swap cases
  sw_direct <- (p_ea == complement(g_ea) & p_oa == complement(g_oa))
  sw_flip   <- (p_ea == complement(g_oa) & p_oa == complement(g_ea))

  flip_idx <- which(flip | sw_flip)
  g_beta[flip_idx] <- -g_beta[flip_idx]

  keep <- direct | flip | sw_direct | sw_flip
  matched <- matched[keep]
  matched$gwas.beta <- g_beta[keep]
  cat(sprintf("After harmonization: %d (direct=%d, flipped=%d)\n",
              nrow(matched), sum(direct[keep]), sum(flip[keep])))

  if (nrow(matched) == 0) { cat("No instruments after harmonization.\n"); next }

  # 5. Select top cis-pQTL (and all matched for IVW)
  setorder(matched, pqtl.pval)
  matched <- matched[!duplicated(gwas.rsid)]

  top_snp <- matched[1]
  cat(sprintf("Instruments: %d (top: %s p=%.2e b_pqtl=%.4f)\n",
              nrow(matched), top_snp$gwas.rsid, top_snp$pqtl.pval, top_snp$pqtl.beta))

  # 6. Wald ratio
  wald_b <- top_snp$gwas.beta / top_snp$pqtl.beta
  wald_se <- abs(wald_b) * sqrt((top_snp$pqtl.se/top_snp$pqtl.beta)^2 +
                                (top_snp$gwas.se/top_snp$gwas.beta)^2)
  wald_p <- 2 * pnorm(-abs(wald_b/wald_se))

  mr_result <- data.table(
    gene = gname, ensg = g$ensg, chr = g$chr,
    top_snp = top_snp$gwas.rsid,
    top_pqtl_p = top_snp$pqtl.pval,
    n_instruments = nrow(matched),
    wald_b = wald_b, wald_se = wald_se, wald_pval = wald_p
  )

  # 7. IVW (if >1 instrument)
  if (nrow(matched) > 1) {
    ratio <- matched$gwas.beta / matched$pqtl.beta
    ratio_se <- abs(ratio) * sqrt((matched$pqtl.se/matched$pqtl.beta)^2 +
                                  (matched$gwas.se/matched$gwas.beta)^2)
    w <- 1/ratio_se^2
    ivw_b <- sum(w*ratio)/sum(w)
    ivw_se <- 1/sqrt(sum(w))
    ivw_p <- 2*pnorm(-abs(ivw_b/ivw_se))
    q_stat <- sum(w*(ratio-ivw_b)^2)
    q_pval <- pchisq(q_stat, df=nrow(matched)-1, lower.tail=FALSE)

    mr_result[, `:=`(ivw_b=ivw_b, ivw_se=ivw_se, ivw_pval=ivw_p,
                     cochran_q=q_stat, cochran_q_pval=q_pval)]
  }

  cat(sprintf("=> %s: Wald b=%.4f p=%.2e (IVW b=%.4f p=%.2e, %d instr)\n",
              gname, wald_b, wald_p,
              ifelse(exists("ivw_b"), ivw_b, NA_real_),
              ifelse(exists("ivw_p"), ivw_p, NA_real_),
              nrow(matched)))

  all_results[[gname]] <- mr_result

  # Save per-gene details
  if (nrow(matched) > 0) {
    fwrite(matched[, .(pqtl.rsid, gwas.rsid, chr, bp.pqtl, bp.gwas,
                       pqtl.ea, pqtl.oa, pqtl.beta, pqtl.se, pqtl.pval,
                       gwas.ea, gwas.oa, gwas.beta, gwas.se, gwas.pval)],
           file.path(OUT_DIR, sprintf("%s_instruments.csv", gname)))
  }
}

# ===== Summary =====
cat("\n\n========== DECODE pQTL MR SUMMARY ==========\n")
if (length(all_results) > 0) {
  combined <- rbindlist(all_results, fill=TRUE)
  combined[, `:=`(
    wald_significant = wald_pval < 0.05,
    ivw_significant  = ifelse(!is.na(ivw_pval), ivw_pval < 0.05, NA)
  )]

  # Sort by Wald p-value
  setorder(combined, wald_pval)

  cols_show <- c("gene","n_instruments","wald_b","wald_se","wald_pval","wald_significant",
                 "ivw_b","ivw_se","ivw_pval","ivw_significant","cochran_q_pval")
  print(combined[, ..cols_show], digits=4)

  fwrite(combined, file.path(OUT_DIR, "phase2_decode_pqtl_mr_combined.csv"))

  cat(sprintf("\nWald significant (p<0.05): %d/%d\n",
              sum(combined$wald_significant, na.rm=TRUE), nrow(combined)))
  if ("ivw_significant" %in% names(combined)) {
    cat(sprintf("IVW significant (p<0.05): %d/%d\n",
                sum(combined$ivw_significant, na.rm=TRUE), nrow(combined)))
  }
} else {
  cat("No results. All genes may lack cis-pQTL instruments.\n")
}

cat(sprintf("\nResults saved to %s\n", file.path(OUT_DIR, "phase2_decode_pqtl_mr_combined.csv")))
cat("Done.\n")
