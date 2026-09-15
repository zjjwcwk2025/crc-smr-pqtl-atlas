#!/usr/bin/env Rscript
# Phase 6c: TCGA COAD survival analysis for key SMR/coloc genes
# Uses pre-existing TCGA COAD TPM + clinical data
# Output: KM curves, Cox HR table

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

# ===== Config =====
EXPR_RDS   <- "/ifs1/User/zhouman/CRC-Analysis-Final/TCGA_COAD_TPM.rds"
CLIN_RDS   <- "/ifs1/User/zhouman/CRC-Analysis-Final/TCGA_COAD_Clinical.rds"
COLOC_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_coloc_top10.csv"
OUT_DIR    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase6_scrna"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===== 1. Load data =====
cat("Loading TCGA COAD data...\n")
expr <- readRDS(EXPR_RDS)   # 59427 genes x 483 samples
clin <- readRDS(CLIN_RDS)   # 437 patients

cat(sprintf("Expression: %d genes x %d samples\n", nrow(expr), ncol(expr)))
cat(sprintf("Clinical: %d patients\n", nrow(clin)))

# Keep only primary tumor samples (01A)
tcga_samples <- colnames(expr)
primary_samples <- tcga_samples[substr(tcga_samples, 14, 15) == "01"]
cat(sprintf("Primary tumor samples: %d\n", length(primary_samples)))

expr <- expr[, primary_samples, drop = FALSE]

# Map sample IDs to patient IDs (first 12 chars)
sample_to_patient <- substr(primary_samples, 1, 12)
colnames(expr) <- sample_to_patient

# Average duplicates (some patients have multiple samples)
dup_patients <- names(which(table(colnames(expr)) > 1))
if (length(dup_patients) > 0) {
  cat(sprintf("Averaging %d patients with multiple samples\n", length(dup_patients)))
  # Keep unique per gene, average
  for (pat in dup_patients) {
    cols <- which(colnames(expr) == pat)
    avg <- rowMeans(expr[, cols, drop = FALSE], na.rm = TRUE)
    # Replace first column, mark others NA
    expr[, cols[1]] <- avg
    expr[, cols[-1]] <- NA
  }
  # Remove NA columns
  expr <- expr[, !is.na(colnames(expr)), drop = FALSE]
  # Remove duplicated colnames
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
}

cat(sprintf("After dedup: %d patients\n", ncol(expr)))

# ===== 2. Define genes to analyze =====
coloc <- fread(COLOC_FILE)
coloc_pass <- coloc[PPH4 > 0.8, gene]

# Key genes: coloc-passing + known drug targets + CAF-expressed
target_genes <- unique(c(
  coloc_pass,                    # 6 coloc-passing
  "CDKN1A", "TNF", "LTA", "NCF2", # pQTL genes
  "MMP11", "SMAD7", "CASC8", "POLD3", "PNKD", # additional key genes
  "BMP2", "GPBAR1", "TMBIM1", "LAMC1", "UTP11", "LIMA1" # coloc genes
))

# Check which are in TCGA
found <- intersect(target_genes, rownames(expr))
missing <- setdiff(target_genes, rownames(expr))
cat(sprintf("\nGenes found in TCGA: %d / %d\n", length(found), length(target_genes)))
if (length(missing) > 0) {
  cat(sprintf("Missing: %s\n", paste(missing, collapse = ", ")))
}

# ===== 3. Expression to log2(TPM+1) =====
expr_log <- log2(as.matrix(expr[found, , drop = FALSE]) + 1)

# ===== 4. Merge with clinical data =====
# Match patient IDs
patients <- colnames(expr_log)
match_idx <- match(patients, clin$Patient_ID)
cat(sprintf("Matched clinical: %d / %d patients\n", sum(!is.na(match_idx)), length(patients)))

clin_matched <- clin[match_idx[!is.na(match_idx)], ]
expr_matched <- expr_log[, !is.na(match_idx), drop = FALSE]

# ===== 5. Survival analysis per gene =====
cat("\n--- Survival Analysis ---\n")

surv_results <- list()

for (g in found) {
  gene_expr <- expr_matched[g, ]

  # Dichotomize by median
  cutoff <- median(gene_expr, na.rm = TRUE)
  group <- ifelse(gene_expr > cutoff, "High", "Low")

  # Create survival data
  surv_df <- data.frame(
    time   = clin_matched$OS.time,
    status = clin_matched$OS,
    group  = factor(group, levels = c("Low", "High")),
    expr   = gene_expr,
    stringsAsFactors = FALSE
  )

  # Remove NAs
  surv_df <- surv_df[complete.cases(surv_df), ]
  if (nrow(surv_df) < 20) next

  # Log-rank test
  fit <- survfit(Surv(time, status) ~ group, data = surv_df)
  lr_test <- survdiff(Surv(time, status) ~ group, data = surv_df)
  lr_p <- 1 - pchisq(lr_test$chisq, df = 1)

  # Cox regression for HR
  cox <- coxph(Surv(time, status) ~ group, data = surv_df)
  cox_sum <- summary(cox)
  hr <- cox_sum$conf.int[1, "exp(coef)"]
  hr_lower <- cox_sum$conf.int[1, "lower .95"]
  hr_upper <- cox_sum$conf.int[1, "upper .95"]
  cox_p <- cox_sum$coefficients[1, "Pr(>|z|)"]

  surv_results[[g]] <- data.frame(
    gene = g,
    n_low = sum(group == "Low"),
    n_high = sum(group == "High"),
    n_total = nrow(surv_df),
    mean_expr = mean(gene_expr),
    hr = hr,
    hr_lower = hr_lower,
    hr_upper = hr_upper,
    cox_p = cox_p,
    logrank_p = lr_p,
    significant = ifelse(lr_p < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %-8s: HR=%.2f (%.2f-%.2f), logrank_p=%.4f, cox_p=%.4f\n",
              g, hr, hr_lower, hr_upper, lr_p, cox_p))
}

# ===== 6. Output survival table =====
surv_df <- rbindlist(surv_results)
surv_df <- surv_df[order(logrank_p)]

cat("\n=== TCGA Survival Results (OS) ===\n")
print(surv_df[, .(gene, n_total, hr, hr_lower, hr_upper, logrank_p, significant)])

fwrite(surv_df, file.path(OUT_DIR, "TableS_tcga_survival.csv"))
cat(sprintf("\nSaved: %s\n", file.path(OUT_DIR, "TableS_tcga_survival.csv")))

# ===== 7. KM Plots for significant genes =====
cat("\nGenerating Kaplan-Meier plots...\n")

sig_genes <- surv_df[logrank_p < 0.05, gene]
if (length(sig_genes) == 0) {
  # Fall back to top 6 by significance
  sig_genes <- head(surv_df$gene, 6)
}

for (g in sig_genes) {
  gene_expr <- expr_matched[g, ]
  cutoff <- median(gene_expr, na.rm = TRUE)
  group <- ifelse(gene_expr > cutoff, "High", "Low")

  surv_input <- data.frame(
    time   = clin_matched$OS.time,
    status = clin_matched$OS,
    group  = factor(group, levels = c("Low", "High"))
  )
  surv_input <- surv_input[complete.cases(surv_input), ]

  fit <- survfit(Surv(time, status) ~ group, data = surv_input)

  p <- ggsurvplot(fit, data = surv_input,
                   pval = TRUE, risk.table = TRUE,
                   palette = c("#2E86C1", "#E74C3C"),
                   title = sprintf("%s — TCGA COAD Overall Survival", g),
                   xlab = "Time (days)", ylab = "Overall Survival",
                   legend.title = sprintf("%s Expression", g),
                   legend.labs = c("Low", "High"),
                   ggtheme = theme_minimal(),
                   risk.table.height = 0.2)

  ggsave(file.path(OUT_DIR, sprintf("FigX_KM_%s.pdf", g)),
         p$plot, width = 6, height = 5)
  cat(sprintf("  KM plot saved: %s\n", g))
}

# ===== 8. Combined KM grid for coloc-passing genes =====
coloc_found <- intersect(coloc_pass, found)
if (length(coloc_found) >= 2) {
  km_plots <- list()
  for (g in coloc_found) {
    gene_expr <- expr_matched[g, ]
    cutoff <- median(gene_expr, na.rm = TRUE)
    group <- ifelse(gene_expr > cutoff, "High", "Low")

    surv_input <- data.frame(
      time   = clin_matched$OS.time,
      status = clin_matched$OS,
      group  = factor(group, levels = c("Low", "High"))
    )
    surv_input <- surv_input[complete.cases(surv_input), ]

    if (length(unique(surv_input$group)) < 2) next

    fit <- survfit(Surv(time, status) ~ group, data = surv_input)
    lr_p <- surv_df[gene == g, logrank_p]
    hr_v <- surv_df[gene == g, hr]

    km_plots[[g]] <- ggsurvplot(fit, data = surv_input,
                                 pval = TRUE, risk.table = FALSE,
                                 palette = c("#2E86C1", "#E74C3C"),
                                 title = sprintf("%s (HR=%.2f, p=%.3f)", g, hr_v, lr_p),
                                 legend = "none",
                                 ggtheme = theme_minimal())$plot
  }

  if (length(km_plots) > 0) {
    combined <- cowplot::plot_grid(plotlist = km_plots, ncol = min(3, length(km_plots)))
    ggsave(file.path(OUT_DIR, "FigX_KM_coloc_panel.pdf"),
           combined, width = 10, height = 3.5 * ceiling(length(km_plots)/3))
    cat("Combined KM panel saved.\n")
  }
}

cat("\nDone.\n")
