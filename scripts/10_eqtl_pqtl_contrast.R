#!/usr/bin/env Rscript
# eQTL vs pQTL 对比矩阵
# 当前: UKB-PPP 5 基因（Phase 1 eQTL SMR + Phase 2 pQTL MR）
# 后续: 加入 deCODE 12 基因

library(dplyr)
library(readr)

# === 1. eQTL SMR+HEIDI results ===
# HEIDI file already contains all SMR columns plus p_HEIDI, nsnp_HEIDI
heidi <- read_tsv("results/phase1_heidi_75genes.smr", show_col_types=FALSE)

# === 2. pQTL MR results (UKB-PPP 5 genes) ===
mr <- read_csv("results/phase2_pqtl/phase2_pqtl_mr_combined.csv", show_col_types=FALSE)

# === 3. Gene name mapping ===
gene_map <- read_tsv("results/phase1_bonferroni_significant_annotated.tsv",
                     show_col_types=FALSE) %>%
  select(ENSG_clean, SYMBOL) %>%
  filter(SYMBOL != "" & !is.na(SYMBOL)) %>%
  distinct()

# === 4. Merge: HEIDI + gene symbols ===
smr_annot <- heidi %>%
  left_join(gene_map, by=c("probeID"="ENSG_clean")) %>%
  filter(!is.na(SYMBOL))

# pQTL MR results (already has gene column)
cat("pQTL MR results:\n")
print(mr %>% select(gene, n_instruments, b, pval, ivw_pval, cochran_q_pval))

# === 5. Build contrast matrix ===
contrast <- mr %>%
  select(gene, n_instruments, pqtl_b=b, pqtl_p=pval,
         ivw_b, ivw_pval, cochran_q_pval) %>%
  left_join(
    smr_annot %>% select(SYMBOL, b_SMR, p_SMR, p_HEIDI, nsnp_HEIDI),
    by=c("gene"="SYMBOL")
  )

cat("\n\n=== eQTL vs pQTL 对比矩阵 (5 UKB-PPP genes) ===\n")
cat(sprintf("%-10s %12s %12s %10s %10s %12s %10s %15s %15s\n",
            "Gene", "eQTL_SMR_b", "eQTL_SMR_p", "HEIDI_p", "HEIDI_pass",
            "pQTL_Wald_b", "pQTL_p", "Direction", "Classification"))

for(i in 1:nrow(contrast)) {
  row <- contrast[i, ]
  smr_b <- row$b_SMR
  pqtl_b <- row$pqtl_b
  pqtl_p <- row$pqtl_p

  # Direction
  if(is.na(smr_b) || is.na(pqtl_b)) {
    dir_label <- "N/A"
    dir_flag <- "?"
  } else if(sign(smr_b) == sign(pqtl_b)) {
    dir_label <- "Concordant"
    dir_flag <- "✓"
  } else {
    dir_label <- "Discordant"
    dir_flag <- "⚠"
  }

  # HEIDI pass
  heidi_p <- row$p_HEIDI
  heidi_pass <- ifelse(is.na(heidi_p), "N/A",
                       ifelse(heidi_p >= 0.01, "PASS", "FAIL"))

  # Classification
  if(row$n_instruments == 0) {
    cls <- "C (无pQTL)"
  } else if(pqtl_p < 0.05 & heidi_pass == "PASS" & dir_flag == "✓") {
    cls <- "A (pQTL独立+)"
  } else if(pqtl_p < 0.05) {
    cls <- "B (confirmatory)"
  } else {
    cls <- "C (pQTL不显著)"
  }

  cat(sprintf("%-10s %+12.4f %12.2e %10.4f %10s %+12.4f %10.4f %15s %15s\n",
              row$gene, smr_b, row$p_SMR, heidi_p, heidi_pass,
              pqtl_b, pqtl_p, paste(dir_flag, dir_label), cls))
}

# === 6. Summary stats ===
cat("\n\n=== 汇总 ===\n")
cat(sprintf("5 基因中:\n"))
cat(sprintf("  pQTL显著(Wald p<0.05): %d (CDKN1A, TNF)\n",
            sum(contrast$pqtl_p < 0.05, na.rm=TRUE)))
cat(sprintf("  HEIDI通过(p>=0.01): %d\n",
            sum(contrast$p_HEIDI >= 0.01, na.rm=TRUE)))
cat(sprintf("  方向一致: %d\n",
            sum(sign(contrast$b_SMR) == sign(contrast$pqtl_b), na.rm=TRUE)))
cat(sprintf("  ⚠ 方向相反: %d (CDKN1A)\n",
            sum(sign(contrast$b_SMR) != sign(contrast$pqtl_b) & contrast$n_instruments > 0, na.rm=TRUE)))
cat(sprintf("  Cochran's Q 显著异质性(p<0.05): %d (需要进一步分析)\n",
            sum(contrast$cochran_q_pval < 0.05, na.rm=TRUE)))

# Save
write_csv(contrast, "results/phase2_eqtl_pqtl_contrast_5genes.csv")
cat("\n✅ 已保存: results/phase2_eqtl_pqtl_contrast_5genes.csv\n")
