#!/usr/bin/env Rscript
# Phase 3: Bayesian colocalization (coloc) — eQTL level
# Pilot: top 10 SMR genes
# Input: eQTLGen text, GWAS .ma, 1000G BIM + afreq reference (MAF from afreq, not GWAS freq)
# Output: results/phase3_coloc_top10.csv

library(data.table)
library(coloc)
library(dplyr)

setDTthreads(8)

# === 1. Load BIM reference (SNP → chr:pos mapping) ===
cat("Loading BIM reference...\n")
bim_file <- "data/1kg_phase3/ref_eur/1000G_EUR.bim"
bim <- fread(bim_file, header=FALSE, select=1:6,
             col.names=c("chr","rsid","cm","pos","a1","a2"))
bim[, chr := as.integer(chr)]
setkey(bim, rsid)
cat(sprintf("BIM: %d SNPs\n", nrow(bim)))

# === 2. Load top 10 genes ===
cat("Loading top genes...\n")
top_genes <- fread("results/phase1_bonferroni_significant_annotated.tsv",
                   select=c("ENSG_clean","SYMBOL","ProbeChr","Probe_bp","p_SMR")) %>%
  filter(SYMBOL != "" & !is.na(SYMBOL)) %>%
  arrange(p_SMR) %>%
  head(10) %>%
  rename(gene_chr=ProbeChr, gene_pos=Probe_bp, ensg=ENSG_clean, gene=SYMBOL) %>%
  mutate(gene_chr=as.integer(gene_chr), gene_pos=as.integer(gene_pos))

cat(sprintf("Top 10 genes for coloc:\n"))
print(top_genes %>% select(gene, gene_chr, gene_pos, p_SMR))

# === 3. Pre-load eQTL data (filter to relevant chromosomes only) ===
eqtl_file <- "data/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
top_chrs <- unique(top_genes$gene_chr)
cat(sprintf("Loading eQTL data for chr: %s...\n", paste(top_chrs, collapse=",")))

eqtl_all <- fread(eqtl_file, header=TRUE,
                  select=c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                           "Zscore","Gene","GeneSymbol","NrSamples"),
                  tmpdir="/ifs1/User/zhouman/tmp")
eqtl_all <- eqtl_all[SNPChr %in% top_chrs]
eqtl_all[, SNPChr := as.integer(SNPChr)]
cat(sprintf("eQTL records loaded: %d\n", nrow(eqtl_all)))

# === 4. Pre-load GWAS data ===
gwas_file <- "data/gwas/GCST90255675.ma"
cat("Loading GWAS data...\n")
gwas_all <- fread(gwas_file, header=TRUE,
                  select=c("SNP","A1","A2","freq","b","se","p","N"))
cat(sprintf("GWAS SNPs loaded: %d\n", nrow(gwas_all)))

# Map GWAS SNPs to chr:pos via BIM
gwas_all <- merge(gwas_all, bim[, .(rsid, chr, pos)], by.x="SNP", by.y="rsid", all.x=FALSE)
gwas_all <- gwas_all[chr %in% top_chrs]
setnames(gwas_all, "chr", "gwas_chr")
setnames(gwas_all, "pos", "gwas_pos")
cat(sprintf("GWAS SNPs with chr:pos: %d\n", nrow(gwas_all)))

# Reference allele frequencies (1000G EUR) for coloc MAF.
# GWAS .ma `freq` is 100% NA, so imputing 0.1 is wrong; use real 1000G EUR MAF.
afreq <- fread("data/1kg_phase3/ref_eur/1000G_EUR.afreq",
               select=c(2,6), col.names=c("SNP","af"))
gwas_all <- merge(gwas_all, afreq, by.x="SNP", by.y="SNP", all.x=FALSE)
cat(sprintf("GWAS SNPs with reference MAF: %d\n", nrow(gwas_all)))

# === 5. Run coloc per gene ===
CIS_WINDOW <- 1e6  # ±1Mb

results <- list()

for(i in 1:nrow(top_genes)) {
  gene <- top_genes$gene[i]
  ensg <- top_genes$ensg[i]
  gchr <- top_genes$gene_chr[i]
  gpos <- top_genes$gene_pos[i]

  cat(sprintf("\n=== %s (chr%d:%d) ===\n", gene, gchr, gpos))

  # 5a. Extract eQTL SNPs for this gene in cis window
  eqtl_gene <- eqtl_all[Gene == ensg &
                         SNPChr == gchr &
                         SNPPos >= gpos - CIS_WINDOW &
                         SNPPos <= gpos + CIS_WINDOW]

  if(nrow(eqtl_gene) < 10) {
    cat(sprintf("  Skipping: only %d eQTL SNPs in cis window\n", nrow(eqtl_gene)))
    next
  }

  # 5b. Extract GWAS SNPs in same region
  gwas_region <- gwas_all[gwas_chr == gchr &
                          gwas_pos >= gpos - CIS_WINDOW &
                          gwas_pos <= gpos + CIS_WINDOW]

  # 5c. Merge by SNP (rsID)
  merged <- merge(eqtl_gene, gwas_region, by.x="SNP", by.y="SNP",
                  suffixes=c("_eqtl","_gwas"))

  if(nrow(merged) < 10) {
    cat(sprintf("  Skipping: only %d overlapping SNPs\n", nrow(merged)))
    next
  }

  cat(sprintf("  eQTL SNPs in cis: %d, GWAS SNPs in cis: %d, Overlap: %d\n",
              nrow(eqtl_gene), nrow(gwas_region), nrow(merged)))

  # 5d. Allele harmonization
  merged <- merged %>%
    mutate(
      # Calculate eQTL beta from Z-score: beta = Z / sqrt(2p(1-p)*(N+Z^2))
      # Simplified: use Z directly, coloc can handle with sdY
      eqtl_beta = Zscore / sqrt(NrSamples),  # approximate
      eqtl_se = 1 / sqrt(NrSamples),
      eqtl_varbeta = eqtl_se^2,

      # GWAS variance
      gwas_varbeta = se^2,

      # Check allele match
      allele_match = (AssessedAllele == A1 & OtherAllele == A2),
      allele_flip = (AssessedAllele == A2 & OtherAllele == A1),

      # Flip eQTL if needed
      eqtl_beta_final = ifelse(allele_flip, -eqtl_beta, eqtl_beta),
      maf = pmin(af, 1 - af)
    ) %>%
    filter((allele_match | allele_flip) & maf > 0 & maf < 1)

  if(nrow(merged) < 10) {
    cat(sprintf("  Skipping: only %d SNPs after allele harmonization\n", nrow(merged)))
    next
  }

  cat(sprintf("  After harmonization: %d SNPs\n", nrow(merged)))

  # 5e. Prepare coloc datasets
  snp_ids <- merged$SNP
  # Cap exact-zero p-values (underflow)
  p_gwas <- pmax(merged$p, 1e-300)
  p_eqtl <- pmax(2 * pnorm(-abs(merged$Zscore)), 1e-300)
  # GWAS: case-control (cc)
  D1 <- list(
    snp = snp_ids,
    pvalues = p_gwas,
    N = merged$N[1],
    type = "cc",
    s = 0.42,  # CRC prevalence ~42% (78573/185616)
    MAF = merged$maf
  )

  # eQTL: quantitative
  D2 <- list(
    snp = snp_ids,
    pvalues = p_eqtl,
    N = merged$NrSamples[1],
    type = "quant",
    MAF = D1$MAF  # same SNPs
  )

  # 5f. Run coloc
  tryCatch({
    col_res <- coloc.abf(dataset1=D1, dataset2=D2)

    # Extract key results
    summary <- col_res$summary
    pp_h4 <- as.numeric(summary[["PP.H4.abf"]])
    pp_h3 <- as.numeric(summary[["PP.H3.abf"]])
    nsnps <- length(snp_ids)

    res_row <- data.frame(
      gene = gene,
      ensg = ensg,
      chr = gchr,
      pos = gpos,
      nsnps = nsnps,
      PPH4 = pp_h4,
      PPH3 = pp_h3,
      p_SMR = top_genes$p_SMR[i],
      stringsAsFactors = FALSE
    )
    results[[length(results)+1]] <- res_row

    flag <- ifelse(pp_h4 > 0.8, "✅ COLOC", ifelse(pp_h4 > 0.5, "⚠ Marginal", "❌ No colocalization"))
    cat(sprintf("  PPH4=%.4f, PPH3=%.4f, nsnp=%d  %s\n", pp_h4, pp_h3, nsnps, flag))

  }, error=function(e) {
    cat(sprintf("  coloc error: %s\n", e$message))
  })
}

# === 6. Output ===
if(length(results) > 0) {
  out <- rbindlist(results)
  out <- out[order(-PPH4)]
  fwrite(out, "results/phase3_coloc_top10.csv")
  cat("\n\n=== Final Results ===\n")
  print(out)
  cat("\n✅ Saved: results/phase3_coloc_top10.csv\n")
} else {
  cat("\n⚠ No genes passed coloc filters\n")
}
