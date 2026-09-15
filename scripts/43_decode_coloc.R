#!/usr/bin/env Rscript
# deCODE pQTL coloc — 5 genes with cis-pQTL instruments
library(data.table)
library(coloc)

GWAS_FILE <- "data/gwas/GCST90255675.h.tsv.gz"
DECODE_DIR <- "data/decode"
CIS_WINDOW <- 1000000L

gene_list <- list(
  STAT6    = list(chr=12, pos=57491338, file="10372_18_STAT6_STAT6.txt.gz"),
  LIMA1    = list(chr=12, pos=50623450, file="11543_84_LIMA1_LIMA1.txt.gz"),
  CCM2     = list(chr=7,  pos=44976160, file="12347_29_CCM2_CCM2.txt.gz"),
  TNFRSF1A = list(chr=12, pos=6291050,  file="2654_19_TNFRSF1A_TNF_sR_I.txt.gz"),
  TNFRSF1B = list(chr=1,  pos=12155565, file="3152_57_TNFRSF1B_TNF_sR_II.txt.gz")
)

cat("Loading GWAS...\n")
gwas <- fread(GWAS_FILE)
setnames(gwas, old=c("chromosome","base_pair_location","effect_allele","other_allele",
                      "beta","standard_error","effect_allele_frequency","p_value"),
          new=c("chr","bp","ea","oa","beta","se","eaf","pval"))
gwas[, varbeta := se^2]
comp <- function(x) chartr("ACGTacgt","TGCAtgca",x)

coloc_results <- list()

for (gname in names(gene_list)) {
  g <- gene_list[[gname]]
  cat(sprintf("\n=== %s (chr%d:%d) ===\n", gname, g$chr, g$pos))

  filter_cmd <- sprintf(
    "zcat %s | awk -F'\\t' 'NR==1 || ($1==\"chr%d\" && $2>=%d && $2<=%d)'",
    file.path(DECODE_DIR, g$file), g$chr, g$pos-CIS_WINDOW, g$pos+CIS_WINDOW)
  pqtl <- fread(cmd=filter_cmd)
  cat(sprintf("cis SNPs in pQTL file: %d\n", nrow(pqtl) - 1))
  if (nrow(pqtl) <= 1) { cat("  -> No cis SNPs\n"); next }

  setnames(pqtl, old=c("Chrom","Pos","effectAllele","otherAllele","Beta","SE","Pval","N","ImpMAF"),
           new=c("chr_raw","bp","ea","oa","beta","se","pval","n","maf"))
  pqtl[, chr := as.integer(gsub("chr","",chr_raw))]
  pqtl[, varbeta := se^2]
  pqtl[, MAF := pmax(maf, 0.001, na.rm=TRUE)]

  setkey(pqtl, chr, bp)
  gw_chr <- gwas[chr == g$chr]
  setkey(gw_chr, bp)

  # Merge by position with tolerance
  matched <- gw_chr[pqtl, on=.(chr, bp), nomatch=0L,
    .(chr, bp.pqtl=i.bp, bp.gwas=bp,
      pqtl.beta=i.beta, pqtl.se=i.se, pqtl.pval=i.pval,
      gwas.beta=beta, gwas.se=se, gwas.pval=pval,
      gwas.ea=ea, gwas.oa=oa, pqtl.ea=i.ea, pqtl.oa=i.oa,
      gwas.maf=eaf, pqtl.maf=i.MAF)]
  cat(sprintf("Exact-position matched: %d\n", nrow(matched)))
  if (nrow(matched) < 10) {
    cat("  -> < 10 matched SNPs, skipping\n")
    next
  }

  # Harmonize
  p_ea <- matched$pqtl.ea; p_oa <- matched$pqtl.oa
  g_ea <- matched$gwas.ea; g_oa <- matched$gwas.oa
  g_beta <- matched$gwas.beta
  direct <- (p_ea==g_ea & p_oa==g_oa)
  flip   <- (p_ea==g_oa & p_oa==g_ea)
  sw_direct <- (p_ea==comp(g_ea) & p_oa==comp(g_oa))
  sw_flip   <- (p_ea==comp(g_oa) & p_oa==comp(g_ea))
  g_beta[flip|sw_flip] <- -g_beta[flip|sw_flip]
  keep <- direct|flip|sw_direct|sw_flip
  matched <- matched[keep]
  matched$gwas.beta <- g_beta[keep]
  cat(sprintf("After harmonization: %d\n", nrow(matched)))
  if (nrow(matched) < 10) {
    cat("  -> < 10 SNPs after harmonization, skipping\n")
    next
  }

  # Filter SNPs with valid MAF (needed for sdY estimation in quant trait)
  maf_ok <- !is.na(matched$pqtl.maf)
  cat(sprintf("With MAF: %d / %d\n", sum(maf_ok), nrow(matched)))
  if (sum(maf_ok) < 10) {
    cat("  -> < 10 SNPs with valid MAF, skipping\n")
    next
  }
  matched_maf <- matched[maf_ok]

  # Build coloc datasets
  d1 <- list(beta=matched_maf$pqtl.beta, varbeta=matched_maf$pqtl.se^2,
             snp=paste0("chr",g$chr,":",matched_maf$bp.pqtl),
             position=matched_maf$bp.pqtl, type="quant",
             MAF=matched_maf$pqtl.maf,
             N=max(pqtl$n,na.rm=TRUE))
  d2 <- list(beta=matched_maf$gwas.beta, varbeta=matched_maf$gwas.se^2,
             snp=paste0("chr",g$chr,":",matched_maf$bp.gwas),
             position=matched_maf$bp.gwas, type="cc",
             N=185616, s=78473/185616)  # GWAS_N = 78473 cases + 107143 controls

  res <- tryCatch(coloc.abf(d1, d2), error=function(e){ cat("ERROR:", e$message, "\n"); NULL })
  if (is.null(res)) {
    cat("  -> coloc.abf failed\n")
    next
  }

  s <- res$summary
  cat(sprintf("PPH0=%.4f PPH1=%.4f PPH2=%.4f PPH3=%.4f PPH4=%.4f\n",
              s["PP.H0.abf"], s["PP.H1.abf"], s["PP.H2.abf"],
              s["PP.H3.abf"], s["PP.H4.abf"]))

  coloc_results[[gname]] <- data.table(
    gene=gname, nsnps=nrow(matched_maf),
    PPH0=s["PP.H0.abf"], PPH1=s["PP.H1.abf"], PPH2=s["PP.H2.abf"],
    PPH3=s["PP.H3.abf"], PPH4=s["PP.H4.abf"])
}

cat("\n\n=== deCODE pQTL coloc Summary ===\n")
if (length(coloc_results) > 0) {
  coloc_combined <- rbindlist(coloc_results)
  coloc_combined[, coloc_pass := PPH4 > 0.8]
  print(coloc_combined, digits=4)
  fwrite(coloc_combined, "results/phase2_decode/phase2_decode_pqtl_coloc.csv")
  cat(sprintf("\nSaved: results/phase2_decode/phase2_decode_pqtl_coloc.csv\n"))
} else {
  cat("No coloc results.\n")
}

cat("\nDone.\n")
