#!/usr/bin/env Rscript
# Phase 3b: Bayesian colocalization (coloc) — pQTL level (UKB-PPP)
# 5 genes: CDKN1A, LTA, NCF2, TET2, TNF
# UKB-PPP = hg38, GWAS = hg38 → direct chr:pos match

library(data.table)
library(coloc)

setDTthreads(8)

# ===== Config =====
UKBPP_DIR <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/ukbppp_merged"
GWAS_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.h.tsv.gz"
OUT_CSV   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3b_coloc_pqtl.csv"
OUT_RDS   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3b_coloc_pqtl.rds"

CIS_WINDOW <- 1e6   # ±1Mb

# UKB-PPP gene info (hg38 positions from script 07)
genes <- list(
  CDKN1A = list(chr = 6,  pos = 36652392, file = "CDKN1A_merged.txt.gz", ensg = "ENSG00000124762"),
  LTA    = list(chr = 6,  pos = 31572054, file = "LTA_merged.txt.gz",    ensg = "ENSG00000226979"),
  NCF2   = list(chr = 1,  pos = 183555568, file = "NCF2_merged.txt.gz",   ensg = "ENSG00000116701"),
  TET2   = list(chr = 4,  pos = 105145428, file = "TET2_merged.txt.gz",   ensg = "ENSG00000168769"),
  TNF    = list(chr = 6,  pos = 31575565, file = "TNF_merged.txt.gz",    ensg = "ENSG00000232810")
)

# ===== Helper: complement allele =====
complement_allele <- function(x) chartr("ACGTacgt", "TGCAtgca", x)

cat("=== coloc @ pQTL level: UKB-PPP → CRC GWAS ===\n\n")

all_coloc <- list()

for (gene_name in names(genes)) {
  g <- genes[[gene_name]]
  cat(sprintf("\n========== %s (chr%d:%d) ==========\n", gene_name, g$chr, g$pos))

  # ---- 1. Load UKB-PPP pQTL (cis ±1Mb) ----
  pqtl_file <- file.path(UKBPP_DIR, g$file)
  if (!file.exists(pqtl_file)) {
    cat(sprintf("  SKIP: missing %s\n", pqtl_file)); next
  }

  filter_cmd <- sprintf("zcat %s | awk 'NR==1 || ($1==%d && $2>=%d && $2<=%d)'",
                        pqtl_file, g$chr, g$pos - CIS_WINDOW, g$pos + CIS_WINDOW)
  pqtl <- fread(cmd = filter_cmd)
  if (nrow(pqtl) == 0) { cat("  No SNPs in cis window\n"); next }

  # Columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
  pqtl[, pval := 10^(-LOG10P)]
  setnames(pqtl,
    old = c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","BETA","SE","LOG10P","N"),
    new = c("chr","bp","snp_id","oa","ea","eaf","beta","se","log10p","n"))

  pqtl_n <- pqtl$n[1]
  cat(sprintf("  pQTL SNPs in cis: %d, N=%d\n", nrow(pqtl), pqtl_n))

  # ---- 2. Load GWAS for same chromosome + cis region ----
  gw_filter <- sprintf("zcat %s | awk 'NR==1 || ($1==%d && $2>=%d && $2<=%d)'",
                       GWAS_FILE, g$chr, g$pos - CIS_WINDOW, g$pos + CIS_WINDOW)
  gwas <- fread(cmd = gw_filter)
  if (nrow(gwas) == 0) { cat("  No GWAS SNPs in cis window\n"); next }

  setnames(gwas,
    old = c("chromosome","base_pair_location","effect_allele","other_allele",
            "beta","standard_error","effect_allele_frequency","p_value","rsid"),
    new = c("chr","bp","ea","oa","beta","se","eaf","pval","snp"))

  # GWAS .h.tsv.gz format has no N column → hardcode total sample size
  gwas_n <- 185616L
  cat(sprintf("  GWAS SNPs in cis: %d, N=%d\n", nrow(gwas), gwas_n))

  # ---- 3. Merge by chr:pos ----
  setkey(pqtl, chr, bp)
  setkey(gwas, chr, bp)
  merged <- merge(pqtl, gwas, by=c("chr","bp"), suffixes=c("_pqtl","_gwas"))
  cat(sprintf("  Overlap by chr:pos: %d\n", nrow(merged)))

  if (nrow(merged) < 10) { cat("  SKIP: <10 overlapping SNPs\n"); next }

  # ---- 4. Allele harmonization ----
  p_ea <- merged$ea_pqtl; p_oa <- merged$oa_pqtl
  g_ea <- merged$ea_gwas; g_oa <- merged$oa_gwas

  direct  <- (p_ea == g_ea & p_oa == g_oa)
  flip    <- (p_ea == g_oa & p_oa == g_ea)
  swgdam  <- (p_ea == complement_allele(g_ea) & p_oa == complement_allele(g_oa))
  swgdam_flip <- (p_ea == complement_allele(g_oa) & p_oa == complement_allele(g_ea))

  # For coloc: need consistent effect allele direction
  keep <- direct | flip | swgdam | swgdam_flip
  merged <- merged[keep]
  # Recompute flip vectors on subsetted table (indices changed)
  p_ea2 <- merged$ea_pqtl; p_oa2 <- merged$oa_pqtl
  g_ea2 <- merged$ea_gwas; g_oa2 <- merged$oa_gwas
  flip2 <- (p_ea2 == g_oa2 & p_oa2 == g_ea2)
  swgdam2 <- (p_ea2 == complement_allele(g_ea2) & p_oa2 == complement_allele(g_oa2))
  swgdam_flip2 <- (p_ea2 == complement_allele(g_oa2) & p_oa2 == complement_allele(g_ea2))
  direct2 <- (p_ea2 == g_ea2 & p_oa2 == g_oa2)

  # Flip pQTL beta where alleles are swapped vs GWAS
  flip_beta <- which(flip2 | swgdam_flip2)
  merged$beta_pqtl[flip_beta] <- -merged$beta_pqtl[flip_beta]
  # swap ea/oa for flipped rows
  tmp_ea <- merged$ea_pqtl[flip_beta]
  merged$ea_pqtl[flip_beta] <- merged$oa_pqtl[flip_beta]
  merged$oa_pqtl[flip_beta] <- tmp_ea

  cat(sprintf("  After harmonization: %d (direct=%d, flipped=%d, swgdam=%d)\n",
    nrow(merged), sum(direct2), sum(flip2), sum(swgdam2 | swgdam_flip2)))

  if (nrow(merged) < 10) { cat("  SKIP: <10 SNPs after harmonization\n"); next }

  # ---- 5. Prepare coloc inputs ----
  snp_ids <- merged$snp
  p_gwas  <- pmax(merged$pval_gwas, 1e-300)
  p_pqtl  <- pmax(merged$pval_pqtl, 1e-300)

  # MAF from UKB-PPP (A1FREQ is freq of ALLELE1/effect allele)
  maf <- pmin(merged$eaf_pqtl, 1 - merged$eaf_pqtl, na.rm = TRUE)
  maf[is.na(maf) | maf <= 0] <- 0.1

  # D1: GWAS (case-control)
  D1 <- list(
    snp = snp_ids,
    pvalues = p_gwas,
    N = gwas_n,
    type = "cc",
    s = 78573/185616,  # CRC cases / total ≈ 0.423
    MAF = maf
  )

  # D2: pQTL (quantitative)
  D2 <- list(
    snp = snp_ids,
    pvalues = p_pqtl,
    N = pqtl_n,
    type = "quant",
    MAF = maf
  )

  # ---- 6. Run coloc ----
  tryCatch({
    col_res <- coloc.abf(dataset1 = D1, dataset2 = D2)

    summary <- col_res$summary
    pp_h4 <- as.numeric(summary[["PP.H4.abf"]])
    pp_h3 <- as.numeric(summary[["PP.H3.abf"]])
    pp_h0 <- as.numeric(summary[["PP.H0.abf"]])
    pp_h1 <- as.numeric(summary[["PP.H1.abf"]])
    pp_h2 <- as.numeric(summary[["PP.H2.abf"]])
    nsnps <- nrow(merged)

    res_row <- data.frame(
      gene = gene_name,
      ensg = g$ensg,
      chr = g$chr,
      pos = g$pos,
      nsnps = nsnps,
      PPH0 = pp_h0,
      PPH1 = pp_h1,
      PPH2 = pp_h2,
      PPH3 = pp_h3,
      PPH4 = pp_h4,
      pqtl_N = pqtl_n,
      gwas_N = gwas_n,
      stringsAsFactors = FALSE
    )

    all_coloc[[gene_name]] <- res_row

    flag <- ifelse(pp_h4 > 0.8, "COLOC",
            ifelse(pp_h4 > 0.5, "Marginal",
              ifelse(pp_h3 > 0.8, "Distinct signals",
                "No colocalization")))
    cat(sprintf("  PPH4=%.4f PPH3=%.4f nsnp=%d  → %s\n", pp_h4, pp_h3, nsnps, flag))

  }, error = function(e) {
    cat(sprintf("  coloc error: %s\n", e$message))
  })
}

# ===== Output =====
cat("\n\n=== Final Results: coloc @ pQTL level ===\n")

if (length(all_coloc) > 0) {
  out <- rbindlist(all_coloc)
  out <- out[order(-PPH4)]

  # Merge eQTL coloc results for comparison
  eqtl_coloc_file <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_coloc_top10.csv"
  if (file.exists(eqtl_coloc_file)) {
    eqtl_coloc <- fread(eqtl_coloc_file)
    eqtl_coloc <- eqtl_coloc[gene %in% names(genes), .(gene, PPH4_eqtl = PPH4, PPH3_eqtl = PPH3)]
    out <- merge(out, eqtl_coloc, by="gene", all.x=TRUE)
  }

  fwrite(out, OUT_CSV)
  saveRDS(all_coloc, OUT_RDS)
  print(out)
  cat(sprintf("\nSaved: %s\n", OUT_CSV))
} else {
  cat("\n  No genes passed coloc filters\n")
}

cat("\nDone.\n")
