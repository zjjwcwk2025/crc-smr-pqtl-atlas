#!/usr/bin/env Rscript
# 60b_spearman_correlation.R
# Compute Spearman correlations for manuscript Fig 2c and Fig 3 (replacing Pearson).
# 铁律 6: 相关用 Spearman（报 rho + exact p）。
# 2026-08-31: switch Pearson -> Spearman for the two effect-size correlations.

setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")

cat("===== Fig 2c: eQTL SMR b vs pQTL Wald b (Spearman) =====\n")
contrast <- read.csv("results/phase2_decode/eqtl_pqtl_contrast_matrix.csv",
                     stringsAsFactors = FALSE)
contrast$b_SMR       <- suppressWarnings(as.numeric(contrast$b_SMR))
contrast$pqtl_wald_b <- suppressWarnings(as.numeric(contrast$pqtl_wald_b))
contrast$pqtl_wald_p <- suppressWarnings(as.numeric(contrast$pqtl_wald_p))
both <- contrast[!is.na(contrast$b_SMR) & !is.na(contrast$pqtl_wald_b) &
                 !is.na(contrast$pqtl_wald_p), ]
cat(sprintf("n genes with both eQTL and pQTL estimates: %d\n", nrow(both)))
ct2 <- cor.test(both$b_SMR, both$pqtl_wald_b, method = "spearman")
cat(sprintf("Spearman rho = %.6f, p = %.6g (two-sided)\n",
            ct2$estimate, ct2$p.value))
cat("genes:", paste(both$SYMBOL, collapse = ", "), "\n\n")

cat("===== Fig 3: GTEx Colon b_SMR vs eQTLGen Blood b_SMR (Spearman) =====\n")
gtex <- read.delim("results/gtex_colon_transverse_smr.smr", stringsAsFactors = FALSE)
e1   <- read.delim("results/phase1_eqtl_smr.smr", stringsAsFactors = FALSE)
e1c  <- data.frame(probeID = e1$probeID, b_SMR_blood = e1$b_SMR,
                   stringsAsFactors = FALSE)
merged <- merge(gtex[, c("probeID", "b_SMR")], e1c, by = "probeID")
merged <- merged[is.finite(merged$b_SMR) & is.finite(merged$b_SMR_blood), ]
cat(sprintf("n genes with SMR estimates in both tissues: %d\n", nrow(merged)))
ct3 <- cor.test(merged$b_SMR, merged$b_SMR_blood, method = "spearman")
cat(sprintf("Spearman rho = %.6f, p = %.6g (two-sided)\n",
            ct3$estimate, ct3$p.value))
concordant_n <- sum(sign(merged$b_SMR) == sign(merged$b_SMR_blood))
cat(sprintf("direction concordant: %d/%d (%.1f%%)\n",
            concordant_n, nrow(merged), 100 * concordant_n / nrow(merged)))

cat("\n===== TCGA survival log-rank p (6 non-MHC Tier 1 genes) =====\n")
surv <- read.csv("results/phase6_scrna/TableS_tcga_survival.csv",
                 stringsAsFactors = FALSE)
tier1 <- c("UTP11", "BMP2", "PNKD", "LAMC1", "TMBIM1", "GPBAR1")
sub <- surv[surv$gene %in% tier1, c("gene", "logrank_p")]
sub <- sub[order(sub$logrank_p), ]
print(sub, row.names = FALSE)
cat(sprintf("min p = %.4f (%s), max p = %.4f (%s)\n",
            min(sub$logrank_p), sub$gene[1],
            max(sub$logrank_p), sub$gene[nrow(sub)]))
