#!/usr/bin/env Rscript
# Phase 2 supplement: deCODE pQTL MR for 4 UKB-PPP overlapping genes + CERS5/ARPC5 relaxed-threshold
library(data.table)

DECODE_DIR <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/decode"
GWAS_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.h.tsv.gz"
OUT_DIR   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase2_decode"
CIS_WINDOW <- 1000000L
PQTL_THRESH <- 5e-8
RELAXED_THRESH <- 5e-6
MATCH_DIST <- 500L

dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# Part A: 4 UKB-PPP overlapping genes with standard threshold
overlap_genes <- list(
  CDKN1A = list(chr=6,  pos=36675805, pqtl_file="18291_8_CDKN1A_p21.txt.gz",   ensg="ENSG00000124762"),
  NCF2   = list(chr=1,  pos=183555562,pqtl_file="10047_12_NCF2_NCF_2.txt.gz",    ensg="ENSG00000116701"),
  TNF    = list(chr=6,  pos=31575563, pqtl_file="5936_53_TNF_TNF_a.txt.gz",     ensg="ENSG00000232810"),
  LTA    = list(chr=6,  pos=31572095, pqtl_file="4703_87_LTA_TNF_b.txt.gz",     ensg="ENSG00000226979")
)

# Part B: CERS5/ARPC5 with relaxed threshold
relaxed_genes <- list(
  CERS5 = list(chr=12, pos=50522623, pqtl_file="13494_6_CERS5_CERS5.txt.gz",     ensg="ENSG00000139624"),
  ARPC5 = list(chr=1,  pos=183620013,pqtl_file="18419_20_ARPC5_p16_ARC.txt.gz",  ensg="ENSG00000162704")
)

complement <- function(x) chartr("ACGTacgt", "TGCAtgca", x)

# Load GWAS
cat("Loading GWAS...\n")
gwas <- fread(GWAS_FILE)
setnames(gwas,
  old=c("chromosome","base_pair_location","effect_allele","other_allele",
        "beta","standard_error","effect_allele_frequency","p_value","rsid"),
  new=c("chr","bp","ea","oa","beta","se","eaf","pval","snp"))
gwas <- gwas[, .(chr, bp, ea, oa, beta, se, eaf, pval, snp)]
setkey(gwas, chr, bp)

# ===== MR runner function =====
run_mr <- function(gname, g, threshold) {
  pqtl_file <- file.path(DECODE_DIR, g$pqtl_file)
  if (!file.exists(pqtl_file)) { cat(sprintf("SKIP: %s not found\n", pqtl_file)); return(NULL) }

  filter_cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==1 || ($1==\"chr%d\" && $2>=%d && $2<=%d)'",
                        pqtl_file, g$chr, g$pos-CIS_WINDOW, g$pos+CIS_WINDOW)
  pqtl <- fread(cmd=filter_cmd)
  if (nrow(pqtl) <= 1) { cat(sprintf("No cis SNPs\n")); return(NULL) }

  setnames(pqtl, old=c("Chrom","Pos","Name","rsids","effectAllele","otherAllele",
                        "Beta","Pval","SE","N","ImpMAF"),
           new=c("chr_raw","bp","name","rsid","ea","oa","beta","pval","se","n","maf"))

  total_cis <- nrow(pqtl)
  pqtl_sig <- pqtl[pval < threshold]
  min_p <- min(pqtl$pval)
  cat(sprintf("  Cis SNPs: %d, sig(p<%.0e): %d, min p=%.2e\n",
              total_cis, threshold, nrow(pqtl_sig), min_p))
  if (nrow(pqtl_sig) == 0) return(NULL)

  pqtl_sig[, chr := as.integer(gsub("chr","",chr_raw))]
  setkey(pqtl_sig, chr, bp)
  gw_chr <- gwas[chr == g$chr]
  setkey(gw_chr, bp)

  matched <- gw_chr[pqtl_sig, on=.(chr, bp), nomatch=0L,
    .(chr, bp.pqtl=i.bp, bp.gwas=bp,
      pqtl.rsid=i.rsid, gwas.rsid=snp,
      pqtl.ea=i.ea, pqtl.oa=i.oa, pqtl.beta=i.beta, pqtl.se=i.se, pqtl.pval=i.pval,
      gwas.ea=ea, gwas.oa=oa, gwas.beta=beta, gwas.se=se, gwas.pval=pval)]

  cat(sprintf("  Exact matched: %d\n", nrow(matched)))
  if (nrow(matched) == 0) return(NULL)

  # Harmonize
  p_ea <- matched$pqtl.ea; p_oa <- matched$pqtl.oa
  g_ea <- matched$gwas.ea; g_oa <- matched$gwas.oa
  g_beta <- matched$gwas.beta
  direct <- (p_ea == g_ea & p_oa == g_oa)
  flip   <- (p_ea == g_oa & p_oa == g_ea)
  sw_direct <- (p_ea == complement(g_ea) & p_oa == complement(g_oa))
  sw_flip   <- (p_ea == complement(g_oa) & p_oa == complement(g_ea))
  g_beta[flip|sw_flip] <- -g_beta[flip|sw_flip]
  keep <- direct|flip|sw_direct|sw_flip
  matched <- matched[keep]
  matched$gwas.beta <- g_beta[keep]
  cat(sprintf("  After harmonization: %d\n", nrow(matched)))
  if (nrow(matched) == 0) return(NULL)

  setorder(matched, pqtl.pval)
  matched <- matched[!duplicated(gwas.rsid)]

  top_snp <- matched[1]

  # Wald ratio
  wald_b <- top_snp$gwas.beta / top_snp$pqtl.beta
  wald_se <- abs(wald_b) * sqrt((top_snp$pqtl.se/top_snp$pqtl.beta)^2 +
                                (top_snp$gwas.se/top_snp$gwas.beta)^2)
  wald_p <- 2 * pnorm(-abs(wald_b/wald_se))

  result <- data.table(
    gene = gname, ensg = g$ensg, chr = g$chr,
    top_snp = top_snp$gwas.rsid,
    top_pqtl_p = top_snp$pqtl.pval,
    n_instruments = nrow(matched),
    wald_b = wald_b, wald_se = wald_se, wald_pval = wald_p,
    wald_significant = wald_p < 0.05,
    min_cis_p = min_p, threshold_used = threshold
  )

  # IVW
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
    result[, `:=`(ivw_b=ivw_b, ivw_se=ivw_se, ivw_pval=ivw_p,
                  cochran_q=q_stat, cochran_q_pval=q_pval,
                  ivw_significant = ivw_p < 0.05)]
  }

  cat(sprintf("  => Wald b=%.4f p=%.2e | %d instruments\n", wald_b, wald_p, nrow(matched)))
  return(result)
}

# ===== Part A: Overlap genes (standard threshold) =====
cat("\n========== Part A: Overlap genes (p<5e-8) ==========\n")
overlap_results <- list()
for (gname in names(overlap_genes)) {
  cat(sprintf("\n--- %s ---\n", gname))
  res <- run_mr(gname, overlap_genes[[gname]], PQTL_THRESH)
  if (!is.null(res)) overlap_results[[gname]] <- res
}

# ===== Part B: Relaxed threshold (CERS5, ARPC5) =====
cat("\n\n========== Part B: Relaxed threshold (p<5e-6) ==========\n")
relaxed_results <- list()
for (gname in names(relaxed_genes)) {
  cat(sprintf("\n--- %s ---\n", gname))
  res <- run_mr(gname, relaxed_genes[[gname]], RELAXED_THRESH)
  if (!is.null(res)) relaxed_results[[gname]] <- res
}

# ===== Summary =====
cat("\n\n========== SUMMARY ==========\n")

print_part <- function(results, title) {
  cat(sprintf("\n--- %s ---\n", title))
  if (length(results) > 0) {
    dt <- rbindlist(results, fill=TRUE)
    print(dt[, .(gene, n_instruments, wald_b, wald_pval, wald_significant, min_cis_p, threshold_used)], digits=4)
    return(dt)
  } else {
    cat("No results.\n")
    return(NULL)
  }
}

all_mr <- list()
a <- print_part(overlap_results, "Overlap genes (standard)")
if (!is.null(a)) all_mr[["overlap"]] <- a
b <- print_part(relaxed_results, "Relaxed threshold")
if (!is.null(b)) all_mr[["relaxed"]] <- b

if (length(all_mr) > 0) {
  combined <- rbindlist(all_mr, fill=TRUE)
  fwrite(combined, file.path(OUT_DIR, "phase2_decode_supplement_mr.csv"))
  cat(sprintf("\nSaved: %s\n", file.path(OUT_DIR, "phase2_decode_supplement_mr.csv")))
}

cat("\nDone.\n")
