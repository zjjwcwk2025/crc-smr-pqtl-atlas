#!/usr/bin/env Rscript
# GTEx Colon Transverse SMR vs eQTLGen Blood SMR 对比
# 回应: 审稿人对"血液eQTL局限性"的质疑
library(dplyr)
library(readr)

# === 1. 读取数据 ===
eqtlgen <- read_tsv("results/phase1_eqtl_smr.smr", show_col_types=FALSE)
gtex <- read_tsv("results/gtex_colon_transverse_smr.smr", show_col_types=FALSE)
annot <- read_tsv("results/phase1_bonferroni_significant_annotated.tsv", show_col_types=FALSE)

# Bonferroni threshold (same as Phase 1)
bonf_thresh <- 3.21e-06

cat("=== GTEx Colon Transverse SMR 结果摘要 ===\n")
cat(sprintf("Total probes tested: %d\n", nrow(gtex)))
cat(sprintf("Bonferroni significant (p < %.2e): %d\n",
            bonf_thresh, sum(gtex$p_SMR < bonf_thresh, na.rm=TRUE)))

# === 2. 方向一致性 ===
eqtlgen_sub <- eqtlgen[eqtlgen$probeID %in% gtex$probeID, ]
merged <- merge(gtex[, c("probeID","Gene","b_SMR","p_SMR","p_HEIDI","nsnp_HEIDI")],
                eqtlgen_sub[, c("probeID","b_SMR","p_SMR","p_HEIDI")],
                by="probeID", suffixes=c("_gtex","_blood"))

merged$direction <- ifelse(merged$b_SMR_gtex * merged$b_SMR_blood > 0,
                           "✓ Concordant", "✗ Discordant")

cat(sprintf("\n配对基因数: %d\n", nrow(merged)))
cat(sprintf("方向一致: %d (%.0f%%)\n",
            sum(merged$direction == "✓ Concordant"),
            100 * sum(merged$direction == "✓ Concordant") / nrow(merged)))

# === 3. 逐个基因对比 ===
cat("\n=== 逐个基因 SMR 对比 (GTEx Colon vs eQTLGen Blood) ===\n")
cat(sprintf("%-15s | %12s %12s | %12s %12s | %s\n",
            "Gene", "b_SMR(colon)", "p_SMR(colon)", "b_SMR(blood)", "p_SMR(blood)", "Direction"))

# Sort by GTEx p-value
merged <- merged[order(merged$p_SMR_gtex), ]

for(i in seq_len(min(nrow(merged), 30))) {
  r <- merged[i, ]
  # Determine significance flags
  gtex_sig <- ifelse(r$p_SMR_gtex < bonf_thresh, "***",
                     ifelse(r$p_SMR_gtex < 0.05, "*", ""))
  blood_sig <- ifelse(r$p_SMR_blood < bonf_thresh, "***",
                      ifelse(r$p_SMR_blood < 0.05, "*", ""))
  dir_sym <- ifelse(r$direction == "✓ Concordant", "✅", "⚠️")

  cat(sprintf("%-15s | %+10.3f %10.2e | %+10.3f %10.2e | %s %s %s\n",
              r$Gene,
              r$b_SMR_gtex, r$p_SMR_gtex,
              r$b_SMR_blood, r$p_SMR_blood,
              r$direction, gtex_sig, blood_sig))
}

# === 4. 重点基因: Tier 1 genes ===
cat("\n\n=== Tier 1 (coloc+SuSiE) 基因在 GTEx Colon 中的表现 ===\n")
tier1_genes <- c("UTP11","BMP2","PNKD","LAMC1","TMBIM1","GPBAR1")
tier1_in_gtex <- merged[merged$Gene %in% tier1_genes, ]

cat(sprintf("%-12s | %12s %12s | %12s %12s | %s | %s\n",
            "Gene", "b_SMR_Colon", "p_SMR_Colon", "b_SMR_Blood", "p_SMR_Blood", "Direction", "HEIDI_Colon"))

for(i in seq_len(nrow(tier1_in_gtex))) {
  r <- tier1_in_gtex[i, ]
  cat(sprintf("%-12s | %+10.3f %10.2e | %+10.3f %10.2e | %s | p=%.2e\n",
              r$Gene,
              r$b_SMR_gtex, r$p_SMR_gtex,
              r$b_SMR_blood, r$p_SMR_blood,
              r$direction, r$p_HEIDI_gtex))
}

# Check which tier 1 genes are MISSING from GTEx Colon
missing_tier1 <- setdiff(tier1_genes, gtex$Gene)
if(length(missing_tier1) > 0) {
  cat(sprintf("\n⚠️ Tier 1 基因未在 GTEx Colon 检测到 cis-eQTL: %s\n",
              paste(missing_tier1, collapse=", ")))
}

# === 5. 统计总结 ===
cat("\n\n=== 统计总结 ===\n")
cat(sprintf("GTEx Colon SMR 显著基因 (p < %.2e): %d\n",
            bonf_thresh, sum(gtex$p_SMR < bonf_thresh, na.rm=TRUE)))
cat(sprintf("eQTLGen Blood SMR 显著基因: %d\n",
            sum(eqtlgen$p_SMR < bonf_thresh, na.rm=TRUE)))

# Both significant
both_sig <- merged$p_SMR_gtex < bonf_thresh & merged$p_SMR_blood < bonf_thresh
cat(sprintf("双队列同时显著: %d\n", sum(both_sig, na.rm=TRUE)))

# GTEx only
gtex_only <- merged$p_SMR_gtex < bonf_thresh & merged$p_SMR_blood >= bonf_thresh
cat(sprintf("仅 GTEx 显著: %d genes: %s\n",
            sum(gtex_only, na.rm=TRUE),
            if(sum(gtex_only)>0) paste(merged$Gene[gtex_only], collapse=", ") else "—"))

# === 6. 效应量相关性（铁律 6: 相关用 Spearman，报 rho + exact p） ===
cor_test <- cor.test(merged$b_SMR_gtex, merged$b_SMR_blood, method = "spearman", use = "complete.obs")
cor_b <- cor_test$estimate
cor_p <- cor_test$p.value
cat(sprintf("\neQTLGen Blood vs GTEx Colon 效应量相关系数 (Spearman rho): %.3f\n", cor_b))
cat(sprintf("Spearman p = %.4g (two-sided)\n", cor_p))

# === 7. HEIDI 对比 ===
cat("\n=== HEIDI 对比 ===\n")
# SMR standard HEIDI threshold: p_HEIDI > 0.01 = pass (no evidence of pleiotropy)
cat(sprintf("GTEx Colon HEIDI pass (p > 0.01): %d/%d\n",
            sum(merged$p_HEIDI_gtex > 0.01, na.rm=TRUE), nrow(merged)))
cat(sprintf("eQTLGen HEIDI pass (p > 0.01): %d/%d\n",
            sum(merged$p_HEIDI_blood > 0.01, na.rm=TRUE), nrow(merged)))

# === 8. 对审稿人的回应 ===
cat("\n\n=== 审稿人回应要点 (Major Concern #2: 血液 eQTL 局限性) ===\n")
cat("Summary for Discussion:\n")
cat(sprintf("  \"To address potential tissue specificity concerns, we performed a sensitivity\n"))
cat(sprintf("  analysis using GTEx Colon Transverse (n=%d eQTL probes) as eQTL reference.\n", nrow(gtex)))
cat(sprintf("  Of %d paired genes tested, %d%% showed concordant effect direction\n",
            nrow(merged), round(100*sum(merged$direction=="✓ Concordant")/nrow(merged))))
cat(sprintf("  between eQTLGen blood and GTEx colon tissue (Spearman rho=%.2f).\n", cor_b))
cat(sprintf("  %d genes reached Bonferroni significance in both tissues,\n",
            sum(both_sig, na.rm=TRUE)))
cat(sprintf("  suggesting that the eQTLGen blood-based SMR findings are broadly\n"))
cat(sprintf("  consistent with colon tissue eQTL effects.\"\n"))

cat("\nLimitations to acknowledge:\n")
cat(sprintf("  - GTEx Colon Transverse sample size (%d genes tested) is smaller than eQTLGen (15,582)\n",
            nrow(gtex)))
cat(sprintf("  - Some Tier 1 genes were not detected in GTEx Colon due to limited cis-eQTL power\n"))
cat(sprintf("  - GTEx colon tissue samples are from non-CRC individuals, so CRC-specific eQTL effects may differ\n"))

# === 9. 落盘输出（消除软支撑：方向一致 21/25 + Spearman rho=0.58 等） ===
out_dir <- "results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n_concordant <- sum(merged$direction == "✓ Concordant")

write_csv(merged, file.path(out_dir, "gtex_colon_eqtlgen_compare_per_gene.csv"))

summary_df <- data.frame(
  metric = c("n_paired_genes",
             "n_direction_concordant",
             "direction_concordant_pct",
             "spearman_rho",
             "spearman_p",
             "n_both_bonferroni_sig",
             "n_gtex_only_bonferroni_sig",
             "n_heidi_pass_gtex",
             "n_heidi_pass_blood"),
  value = c(nrow(merged),
            n_concordant,
            round(100 * n_concordant / nrow(merged), 2),
            round(cor_b, 4),
            signif(cor_p, 4),
            sum(both_sig, na.rm = TRUE),
            sum(gtex_only, na.rm = TRUE),
            sum(merged$p_HEIDI_gtex > 0.01, na.rm = TRUE),
            sum(merged$p_HEIDI_blood > 0.01, na.rm = TRUE))
)
write_csv(summary_df, file.path(out_dir, "gtex_colon_eqtlgen_compare_summary.csv"))

cat(sprintf("\n落盘完成: %s\n",
            "results/gtex_colon_eqtlgen_compare_summary.csv"))

cat("\n=== 完成 ===\n")
