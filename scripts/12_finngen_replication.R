#!/usr/bin/env Rscript
# Phase 4: FinnGen R13 C3_COLORECTAL dual-cohort replication
# Uses eQTL instruments from top SMR genes → MR in FinnGen CRC
# Compares direction consistency with discovery GWAS (GCST90255675)

library(data.table)
library(dplyr)

setDTthreads(8)

# === 1. Load BIM reference ===
cat("Loading BIM reference...\n")
bim <- fread("data/1kg_phase3/ref_eur/1000G_EUR.bim", header=FALSE, select=1:6,
             col.names=c("chr","rsid","cm","pos","a1","a2"))
bim[, chr := as.integer(chr)]
setkey(bim, rsid)

# === 2. Load top genes (all 10 from coloc analysis) ===
cat("Loading top genes...\n")
top_genes <- fread("results/phase3_coloc_top10.csv") %>%
  mutate(gene_chr = chr, gene_pos = pos)

cat(sprintf("Genes for FinnGen replication: %d\n", nrow(top_genes)))
print(top_genes %>% select(gene, gene_chr, gene_pos, PPH4, p_SMR))

# === 3. Load eQTL data ===
top_chrs <- unique(top_genes$gene_chr)
cat(sprintf("Loading eQTL data for chr: %s...\n", paste(top_chrs, collapse=",")))

eqtl_file <- "data/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
eqtl_all <- fread(eqtl_file, header=TRUE,
                  select=c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                           "Zscore","Gene","GeneSymbol","NrSamples","Pvalue"),
                  tmpdir="/ifs1/User/zhouman/tmp")
eqtl_all <- eqtl_all[SNPChr %in% top_chrs]
eqtl_all[, SNPChr := as.integer(SNPChr)]
cat(sprintf("eQTL records loaded: %d\n", nrow(eqtl_all)))

# === 4. Load FinnGen GWAS ===
cat("Loading FinnGen R13 C3_COLORECTAL...\n")
finngen_file <- "data/gwas/finngen_R13_C3_COLORECTAL.gz"
finngen <- fread(finngen_file, header=TRUE,
                 select=c("#chrom","pos","ref","alt","rsids","beta","sebeta","pval","af_alt"),
                 tmpdir="/ifs1/User/zhouman/tmp")
setnames(finngen, "#chrom", "chr")
finngen[, chr := as.integer(chr)]

# Handle multi-allelic rsids (comma-separated) — take the first
finngen[, rsid := gsub(",.*", "", rsids)]
cat(sprintf("FinnGen SNPs loaded: %d\n", nrow(finngen)))

# === 5. MR per gene ===
CIS_WINDOW <- 1e6

results <- list()

for(i in 1:nrow(top_genes)) {
  gene <- top_genes$gene[i]
  ensg <- top_genes$ensg[i]
  gchr <- top_genes$gene_chr[i]
  gpos <- top_genes$gene_pos[i]
  pp_h4 <- top_genes$PPH4[i]
  p_smr <- top_genes$p_SMR[i]

  cat(sprintf("\n=== %s (chr%d:%d, coloc PPH4=%.4f) ===\n", gene, gchr, gpos, pp_h4))

  # 5a. Extract eQTL SNPs as instruments in cis window
  eqtl_gene <- eqtl_all[Gene == ensg &
                         SNPChr == gchr &
                         SNPPos >= gpos - CIS_WINDOW &
                         SNPPos <= gpos + CIS_WINDOW]

  if(nrow(eqtl_gene) < 3) {
    cat(sprintf("  Skipping: only %d eQTL SNPs in cis window\n", nrow(eqtl_gene)))
    next
  }

  # Select eQTL instruments with strong signal (p < 5e-4 in cis)
  # This gives us reasonable instrument candidates
  eqtl_inst <- eqtl_gene[Pvalue < 5e-4]

  if(nrow(eqtl_inst) < 1) {
    cat(sprintf("  Skipping: no SNPs with p < 5e-4\n"))
    next
  }

  cat(sprintf("  eQTL instruments (p<5e-4): %d\n", nrow(eqtl_inst)))

  # 5b. Match with FinnGen by chr:pos (hg38 → FinnGen is hg38)
  finngen_region <- finngen[chr == gchr &
                            pos >= gpos - CIS_WINDOW &
                            pos <= gpos + CIS_WINDOW]

  # Match by chr:pos
  merged <- merge(eqtl_inst, finngen_region,
                  by.x=c("SNPChr","SNPPos"), by.y=c("chr","pos"),
                  suffixes=c("_eqtl","_finngen"))

  if(nrow(merged) < 1) {
    cat(sprintf("  Skipping: no FinnGen matches by chr:pos\n"))
    next
  }

  cat(sprintf("  FinnGen matches: %d\n", nrow(merged)))

  # 5c. Allele harmonization (with strand-flip handling)
  # Complement mapping: A↔T, C↔G
  complement <- function(x) {
    case_when(
      x == "A" ~ "T", x == "T" ~ "A",
      x == "C" ~ "G", x == "G" ~ "C",
      TRUE ~ NA_character_
    )
  }

  merged <- merged %>%
    mutate(
      # Scenario 1: direct match (AssessedAllele=ref, OtherAllele=alt)
      allele_match = (AssessedAllele == ref & OtherAllele == alt),
      # Scenario 2: allele flip (AssessedAllele=alt, OtherAllele=ref)
      allele_flip  = (AssessedAllele == alt & OtherAllele == ref),
      # Scenario 3: strand flip - complement match
      allele_strand = (complement(AssessedAllele) == ref & complement(OtherAllele) == alt),
      # Scenario 4: strand flip + allele flip
      allele_strand_flip = (complement(AssessedAllele) == alt & complement(OtherAllele) == ref),

      # Any valid harmonization
      harmonized = allele_match | allele_flip | allele_strand | allele_strand_flip,

      # FinnGen beta is for the alt allele; flip if allele order reversed
      finngen_beta_harm = case_when(
        allele_match | allele_strand ~ beta,          # no flip needed
        allele_flip | allele_strand_flip ~ -beta,      # flip needed
        TRUE ~ NA_real_
      ),

      # eQTL beta approximation
      eqtl_beta = Zscore / sqrt(NrSamples),
      eqtl_beta_harm = case_when(
        allele_match | allele_strand ~ eqtl_beta,
        allele_flip | allele_strand_flip ~ -eqtl_beta,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(harmonized)

  if(nrow(merged) < 1) {
    cat(sprintf("  Skipping: no SNPs after allele harmonization\n"))
    next
  }

  cat(sprintf("  After harmonization: %d SNPs\n", nrow(merged)))

  # 5d. MR analysis
  # We got r in discovery SMR (b_SMR), and r_fg in FinnGen
  # Direction comparison: sign(b_SMR) vs sign(b_finngen)

  # First, check if SMR direction data is available
  # We'll read it from the SMR output file
  n_instruments <- nrow(merged)

  # Wald ratio for each SNP
  merged <- merged %>%
    mutate(
      wald_b = finngen_beta_harm / eqtl_beta_harm,
      wald_se = abs(sebeta / eqtl_beta_harm),
      wald_p = 2 * pnorm(-abs(wald_b / wald_se)),
      # Direction concordance: same sign?
      dir_concordant = sign(eqtl_beta_harm) == sign(finngen_beta_harm)
    )

  # IVW meta-analysis of individual Wald ratios
  if(n_instruments >= 2) {
    # IVW: weighted average
    weights <- 1 / (merged$wald_se^2)
    ivw_b <- sum(merged$wald_b * weights) / sum(weights)
    ivw_se <- sqrt(1 / sum(weights))
    ivw_p <- 2 * pnorm(-abs(ivw_b / ivw_se))

    # Cochran's Q for heterogeneity
    Q <- sum(weights * (merged$wald_b - ivw_b)^2)
    Q_pval <- pchisq(Q, df = n_instruments - 1, lower.tail = FALSE)
  } else {
    ivw_b <- merged$wald_b[1]
    ivw_se <- merged$wald_se[1]
    ivw_p <- merged$wald_p[1]
    Q <- NA_real_
    Q_pval <- NA_real_
  }

  # Median Wald ratio
  median_b <- median(merged$wald_b)
  median_p <- median(merged$wald_p)

  # Concordance summary
  n_concordant <- sum(merged$dir_concordant)
  n_discordant <- sum(!merged$dir_concordant)

  # We also need the discovery SMR b. Let me try to read it.
  # For now, we compare FinnGen direction internally.

  res_row <- data.frame(
    gene = gene,
    ensg = ensg,
    chr = gchr,
    pos = gpos,
    coloc_PPH4 = pp_h4,
    p_SMR = p_smr,
    n_instruments = n_instruments,
    n_concordant = n_concordant,
    n_discordant = n_discordant,
    ivw_b = ivw_b,
    ivw_se = ivw_se,
    ivw_p = ivw_p,
    cochran_q_pval = Q_pval,
    median_wald_b = median_b,
    median_wald_p = median_p,
    stringsAsFactors = FALSE
  )

  results[[length(results)+1]] <- res_row

  flag <- ifelse(ivw_p < 0.05, "✅ REPLICATED",
                 ifelse(ivw_p < 0.1, "⚠ Nominal", "❌ Not replicated"))
  cat(sprintf("  IVW: b=%.4f, se=%.4f, p=%.2e  %s\n",
              ivw_b, ivw_se, ivw_p, flag))
  cat(sprintf("  Direction: %d concordant / %d discordant\n",
              n_concordant, n_discordant))
  cat(sprintf("  Cochran Q p=%.2e\n", Q_pval))
}

# === 6. Output ===
if(length(results) > 0) {
  out <- rbindlist(results)

  # Now add SMR direction comparison
  # Read SMR results for discovery direction
  smr_file <- "results/phase1_bonferroni_significant_annotated.tsv"
  smr <- fread(smr_file, select=c("SYMBOL","b_SMR","p_SMR")) %>%
    filter(SYMBOL != "" & !is.na(SYMBOL))

  out <- out %>%
    left_join(smr %>% select(SYMBOL, b_SMR_discovery = b_SMR),
              by=c("gene"="SYMBOL")) %>%
    mutate(
      discovery_dir = sign(b_SMR_discovery),
      finngen_dir = sign(ivw_b),
      dir_match = case_when(
        is.na(discovery_dir) | is.na(finngen_dir) ~ "N/A",
        discovery_dir == finngen_dir ~ "✅ Concordant",
        TRUE ~ "⚠ Discordant"
      ),
      replication = case_when(
        ivw_p < 0.05 & dir_match == "✅ Concordant" ~ "Full",
        ivw_p < 0.05 & dir_match == "⚠ Discordant" ~ "Direction mismatch",
        ivw_p < 0.1 ~ "Nominal",
        TRUE ~ "Not replicated"
      )
    )

  out <- out[order(-coloc_PPH4)]
  fwrite(out, "results/phase4_finngen_replication.csv")

  cat("\n\n=== FinnGen Replication Results ===\n")
  cat(sprintf("%-10s %6s %10s %10s %10s %12s %18s %14s\n",
              "Gene", "nsnp", "PPH4", "IVW_b", "IVW_p", "CochranQ_p",
              "Direction", "Verdict"))
  for(i in 1:nrow(out)) {
    row <- out[i,]
    cat(sprintf("%-10s %6d %10.4f %10.4f %10.2e %12.2e %18s %14s\n",
                row$gene, row$n_instruments, row$coloc_PPH4,
                row$ivw_b, row$ivw_p, row$cochran_q_pval,
                row$dir_match, row$replication))
  }

  # Summary
  n_replicated <- sum(out$replication == "Full", na.rm=TRUE)
  n_total <- nrow(out)
  cat(sprintf("\n=== Summary ===\n"))
  cat(sprintf("Full replication (p<0.05 + concordant): %d/%d (%.0f%%)\n",
              n_replicated, n_total, 100*n_replicated/n_total))
  cat(sprintf("Nominal (p<0.1): %d\n", sum(out$replication == "Nominal", na.rm=TRUE)))
  cat(sprintf("Direction mismatch: %d\n", sum(out$replication == "Direction mismatch", na.rm=TRUE)))
  cat(sprintf("Not replicated: %d\n", sum(out$replication == "Not replicated", na.rm=TRUE)))

  cat("\n✅ Saved: results/phase4_finngen_replication.csv\n")
} else {
  cat("\n⚠ No genes passed replication filters\n")
}
