#!/usr/bin/env Rscript
# Phase 3b: Bayesian colocalization (coloc) — eQTL level — ALL 75 genes
# Extension of scripts/11_coloc_eqtl.R from Top 10 → full 75
# Input: eQTLGen text, GWAS .ma, 1000G BIM + afreq reference (MAF from afreq, not GWAS freq)
# Output: results/phase3_coloc_all75.csv

library(data.table)
library(coloc)
library(dplyr)

setDTthreads(16)

proj <- "/ifs1/User/zhouman/project9-v5-crc-atlas"

# === 1. Load BIM reference (SNP → chr:pos mapping) ===
cat("Loading BIM reference...\n")
bim <- fread(file.path(proj, "data/1kg_phase3/ref_eur/1000G_EUR.bim"),
             header=FALSE, select=1:6,
             col.names=c("chr","rsid","cm","pos","a1","a2"))
bim[, chr := as.integer(chr)]
setkey(bim, rsid)
cat(sprintf("BIM: %d SNPs\n", nrow(bim)))

# === 2. Load ALL 75 Bonferroni-significant genes ===
cat("Loading significant genes...\n")
sig <- fread(file.path(proj, "results/phase1_bonferroni_significant_annotated.tsv"),
             select=c("ENSG_clean","SYMBOL","ProbeChr","Probe_bp","p_SMR")) %>%
  rename(gene_chr=ProbeChr, gene_pos=Probe_bp, ensg=ENSG_clean, gene=SYMBOL) %>%
  mutate(gene_chr=as.integer(gene_chr), gene_pos=as.integer(gene_pos))

# For genes without SYMBOL, use ENSG as label
sig[is.na(gene) | gene == "", gene := paste0(substr(ensg, 1, 15), "...")]

cat(sprintf("All %d genes for coloc:\n", nrow(sig)))
print(sig %>% select(gene, gene_chr, gene_pos, p_SMR))

# === 3. Pre-load eQTL data (filter to relevant chromosomes) ===
eqtl_file <- file.path(proj,
  "data/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
target_chrs <- sort(unique(sig$gene_chr))
cat(sprintf("Loading eQTL data for %d chromosomes: %s...\n",
            length(target_chrs), paste(target_chrs, collapse=",")))

eqtl_all <- fread(eqtl_file, header=TRUE,
                  select=c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                           "Zscore","Gene","GeneSymbol","NrSamples"),
                  tmpdir="/ifs1/User/zhouman/tmp")
eqtl_all <- eqtl_all[SNPChr %in% target_chrs]
eqtl_all[, SNPChr := as.integer(SNPChr)]
cat(sprintf("eQTL records loaded: %d\n", nrow(eqtl_all)))

# === 4. Pre-load GWAS data ===
gwas_file <- file.path(proj, "data/gwas/GCST90255675.ma")
cat("Loading GWAS data...\n")
gwas_all <- fread(gwas_file, header=TRUE,
                  select=c("SNP","A1","A2","freq","b","se","p","N"))
cat(sprintf("GWAS SNPs loaded: %d\n", nrow(gwas_all)))

# Map GWAS SNPs to chr:pos via BIM
gwas_all <- merge(gwas_all, bim[, .(rsid, chr, pos)],
                  by.x="SNP", by.y="rsid", all.x=FALSE)
gwas_all <- gwas_all[chr %in% target_chrs]
setnames(gwas_all, "chr", "gwas_chr")
setnames(gwas_all, "pos", "gwas_pos")
cat(sprintf("GWAS SNPs with chr:pos: %d\n", nrow(gwas_all)))

# Reference allele frequencies (1000G EUR) for coloc MAF.
# GWAS .ma `freq` is 100% NA, so imputing 0.1 is wrong; use real 1000G EUR MAF.
afreq <- fread(file.path(proj, "data/1kg_phase3/ref_eur/1000G_EUR.afreq"),
               select=c(2,6), col.names=c("SNP","af"))
gwas_all <- merge(gwas_all, afreq, by.x="SNP", by.y="SNP", all.x=FALSE)
cat(sprintf("GWAS SNPs with reference MAF: %d\n", nrow(gwas_all)))
rm(bim, afreq); gc()  # free memory

# === 5. Run coloc per gene ===
CIS_WINDOW <- 1e6  # ±1Mb

results <- list()
gene_count <- nrow(sig)

for(i in 1:gene_count) {
  gene_label <- sig$gene[i]
  ensg      <- sig$ensg[i]
  gchr      <- sig$gene_chr[i]
  gpos      <- sig$gene_pos[i]

  cat(sprintf("\n[%d/%d] %s (chr%d:%d)", i, gene_count, gene_label, gchr, gpos))

  # 5a. Extract eQTL SNPs for this gene in cis window
  eqtl_gene <- eqtl_all[Gene == ensg &
                         SNPChr == gchr &
                         SNPPos >= gpos - CIS_WINDOW &
                         SNPPos <= gpos + CIS_WINDOW]

  if(nrow(eqtl_gene) < 10) {
    cat(sprintf(" → SKIP: only %d eQTL SNPs\n", nrow(eqtl_gene)))
    results[[length(results)+1]] <- data.frame(
      gene=gene_label, ensg=ensg, chr=gchr, pos=gpos,
      nsnps=0, PPH4=NA, PPH3=NA, p_SMR=sig$p_SMR[i],
      status="insufficient_eqtl", stringsAsFactors=FALSE
    )
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
    cat(sprintf(" → SKIP: only %d overlapping SNPs\n", nrow(merged)))
    results[[length(results)+1]] <- data.frame(
      gene=gene_label, ensg=ensg, chr=gchr, pos=gpos,
      nsnps=nrow(merged), PPH4=NA, PPH3=NA, p_SMR=sig$p_SMR[i],
      status="insufficient_overlap", stringsAsFactors=FALSE
    )
    next
  }

  # 5d. Allele harmonization
  merged <- merged %>%
    mutate(
      eqtl_beta    = Zscore / sqrt(NrSamples),
      eqtl_se      = 1 / sqrt(NrSamples),
      eqtl_varbeta = eqtl_se^2,
      gwas_varbeta = se^2,
      allele_match = (AssessedAllele == A1 & OtherAllele == A2),
      allele_flip  = (AssessedAllele == A2 & OtherAllele == A1),
      eqtl_beta_final = ifelse(allele_flip, -eqtl_beta, eqtl_beta),
      maf = pmin(af, 1 - af)
    ) %>%
    filter((allele_match | allele_flip) & maf > 0 & maf < 1)

  if(nrow(merged) < 10) {
    cat(sprintf(" → SKIP: only %d SNPs after harmonization\n", nrow(merged)))
    results[[length(results)+1]] <- data.frame(
      gene=gene_label, ensg=ensg, chr=gchr, pos=gpos,
      nsnps=nrow(merged), PPH4=NA, PPH3=NA, p_SMR=sig$p_SMR[i],
      status="insufficient_harmonized", stringsAsFactors=FALSE
    )
    next
  }

  nsnps_harm <- nrow(merged)
  cat(sprintf(" → %d harmonized SNPs", nsnps_harm))

  # 5e. Prepare coloc datasets
  snp_ids <- merged$SNP
  p_gwas  <- pmax(merged$p, 1e-300)
  p_eqtl  <- pmax(2 * pnorm(-abs(merged$Zscore)), 1e-300)

  D1 <- list(
    snp     = snp_ids,
    pvalues = p_gwas,
    N       = merged$N[1],
    type    = "cc",
    s       = 0.42,
    MAF     = merged$maf
  )

  D2 <- list(
    snp     = snp_ids,
    pvalues = p_eqtl,
    N       = merged$NrSamples[1],
    type    = "quant",
    MAF     = D1$MAF
  )

  # 5f. Run coloc
  tryCatch({
    col_res <- coloc.abf(dataset1=D1, dataset2=D2)
    summary <- col_res$summary
    pp_h4   <- as.numeric(summary[["PP.H4.abf"]])
    pp_h3   <- as.numeric(summary[["PP.H3.abf"]])

    flag <- ifelse(pp_h4 > 0.8, "PASS",
            ifelse(pp_h4 > 0.5, "MARGINAL",
            ifelse(pp_h4 > 0.01, "WEAK", "FAIL")))

    cat(sprintf(" → PPH4=%.4f %s", pp_h4, flag))

    results[[length(results)+1]] <- data.frame(
      gene   = gene_label,
      ensg   = ensg,
      chr    = gchr,
      pos    = gpos,
      nsnps  = nsnps_harm,
      PPH4   = pp_h4,
      PPH3   = pp_h3,
      p_SMR  = sig$p_SMR[i],
      status = flag,
      stringsAsFactors = FALSE
    )
  }, error=function(e) {
    cat(sprintf(" → ERROR: %s", e$message))
    results[[length(results)+1]] <- data.frame(
      gene=gene_label, ensg=ensg, chr=gchr, pos=gpos,
      nsnps=nsnps_harm, PPH4=NA, PPH3=NA, p_SMR=sig$p_SMR[i],
      status=paste0("error:", substr(e$message,1,60)),
      stringsAsFactors=FALSE
    )
  })
}

# === 6. Output ===
cat("\n\n=== Final Results ===\n")
out <- rbindlist(results)
out <- out[order(-PPH4, na.last=TRUE)]

# Summary stats
n_pass <- sum(out$status == "PASS" & !is.na(out$status), na.rm=TRUE)
n_marg <- sum(out$status == "MARGINAL" & !is.na(out$status), na.rm=TRUE)
n_weak <- sum(out$status == "WEAK" & !is.na(out$status), na.rm=TRUE)
n_fail <- sum(out$status == "FAIL" & !is.na(out$status), na.rm=TRUE)
n_skip <- sum(grepl("insufficient|error", out$status, ignore.case=TRUE))
n_total <- nrow(out)

cat(sprintf("\nTotal genes tested: %d\n", n_total))
cat(sprintf("  PPH4 > 0.8:  %d (%.1f%%)\n", n_pass, 100*n_pass/n_total))
cat(sprintf("  PPH4 0.5-0.8: %d (%.1f%%)\n", n_marg, 100*n_marg/n_total))
cat(sprintf("  PPH4 0.01-0.5: %d (%.1f%%)\n", n_weak, 100*n_weak/n_total))
cat(sprintf("  PPH4 < 0.01:  %d (%.1f%%)\n", n_fail, 100*n_fail/n_total))
cat(sprintf("  Failed/Skipped: %d\n", n_skip))
cat(sprintf("  Valid coloc calls: %d\n", n_total - n_skip))

print(out)

out_file <- file.path(proj, "results/phase3_coloc_all75.csv")
fwrite(out, out_file)
cat(sprintf("\n✅ Saved: %s\n", out_file))
