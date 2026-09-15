#!/usr/bin/env Rscript
# Phase 6d: TCGA COAD+READ paired tumor-normal differential expression
# Data: UCSC Xena HiSeqV2 (log2(RSEM+1) TPM)
# Output: DE boxplots, paired line plots, volcano-style summary

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})

# ===== Config =====
EXPR_FILE  <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/tcga_coad/TCGA_COADREAD_HiSeqV2.gz"
COLOC_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_coloc_top10.csv"
OUT_DIR    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase6_scrna"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===== 1. Load and prepare data =====
cat("Loading TCGA COADREAD data from UCSC Xena...\n")
expr_raw <- fread(EXPR_FILE)
# First column is gene symbol, rest are TCGA samples
genes_vec <- expr_raw[[1]]
expr_mat <- as.matrix(expr_raw[, -1, with = FALSE])
rownames(expr_mat) <- genes_vec

cat(sprintf("Expression: %d genes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))

# Parse sample types from TCGA barcode (chars 14-15)
samples <- colnames(expr_mat)
sample_types <- substr(samples, 14, 15)
patient_ids <- substr(samples, 1, 12)

cat(sprintf("Primary Tumor (01): %d\n", sum(sample_types == "01")))
cat(sprintf("Solid Normal (11): %d\n", sum(sample_types == "11")))

# Identify paired patients
tumor_patients <- patient_ids[sample_types == "01"]
normal_patients <- patient_ids[sample_types == "11"]
paired_patients <- intersect(tumor_patients, normal_patients)
cat(sprintf("Paired patients: %d\n", length(paired_patients)))

# Select tumor and normal columns (all, not just paired)
tumor_cols <- which(sample_types == "01")
normal_cols <- which(sample_types == "11")

# ===== 2. Define genes =====
coloc <- fread(COLOC_FILE)
coloc_pass <- coloc[PPH4 > 0.8, gene]

target_genes <- unique(c(
  coloc_pass,                          # 6 coloc-passing
  "CDKN1A", "TNF", "LTA", "NCF2",     # pQTL genes
  "MMP11", "SMAD7", "CASC8", "POLD3"  # key SMR genes
))

found <- intersect(target_genes, genes_vec)
cat(sprintf("\nGenes in TCGA data: %d / %d\n", length(found), length(target_genes)))
missing <- setdiff(target_genes, genes_vec)
if (length(missing) > 0) cat(sprintf("Missing: %s\n", paste(missing, collapse = ", ")))

# ===== 3. Differential expression (all tumor vs all normal) =====
cat("\n=== Differential Expression (Tumor vs Normal) ===\n")

de_results <- list()

for (g in found) {
  t_expr <- expr_mat[g, tumor_cols]
  n_expr <- expr_mat[g, normal_cols]

  log2fc <- mean(t_expr) - mean(n_expr)  # Xena data is already log2(RSEM+1)
  test <- wilcox.test(t_expr, n_expr, paired = FALSE)
  p_val <- test$p.value

  de_results[[g]] <- data.frame(
    gene = g,
    mean_tumor = mean(t_expr),
    mean_normal = mean(n_expr),
    log2FC = log2fc,
    p_value = p_val,
    direction = ifelse(log2fc > 0, "Up", "Down"),
    significant = ifelse(p_val < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %-8s: Tumor=%.2f Normal=%.2f log2FC=%+.2f p=%.2e %s\n",
              g, mean(t_expr), mean(n_expr), log2fc, p_val,
              ifelse(p_val < 0.05, "*", "")))
}

de_df <- rbindlist(de_results)
de_df <- de_df[order(p_value)]

cat("\n=== DE Results Summary ===\n")
print(de_df[, .(gene, log2FC, p_value, direction, significant)])

fwrite(de_df, file.path(OUT_DIR, "TableS_tcga_tumor_vs_normal.csv"))

# ===== 4. Paired DE (for patients with matched tumor-normal) =====
if (length(paired_patients) >= 5) {
  cat(sprintf("\n=== Paired DE (%d pairs) ===\n", length(paired_patients)))

  # Align paired samples
  paired_tumor_cols <- which(sample_types == "01" & patient_ids %in% paired_patients)
  paired_normal_cols <- which(sample_types == "11" & patient_ids %in% paired_patients)

  # Order by patient ID for proper pairing
  t_order <- order(patient_ids[paired_tumor_cols])
  n_order <- order(patient_ids[paired_normal_cols])
  paired_tumor_cols <- paired_tumor_cols[t_order]
  paired_normal_cols <- paired_normal_cols[n_order]

  paired_de <- list()
  for (g in found) {
    t_expr <- expr_mat[g, paired_tumor_cols]
    n_expr <- expr_mat[g, paired_normal_cols]
    fc <- mean(t_expr - n_expr)
    test <- wilcox.test(t_expr, n_expr, paired = TRUE)
    paired_de[[g]] <- data.frame(
      gene = g, paired_log2FC = fc, paired_p = test$p.value,
      stringsAsFactors = FALSE
    )
  }
  paired_df <- rbindlist(paired_de)
  paired_df <- paired_df[order(paired_p)]
  cat("\nPaired DE (top):\n")
  print(head(paired_df, 10))

  # Merge with main DE
  de_df <- merge(de_df, paired_df, by = "gene", all.x = TRUE)
  fwrite(de_df, file.path(OUT_DIR, "TableS_tcga_tumor_vs_normal.csv"))
}

# ===== 5. Boxplots for significant genes =====
cat("\nGenerating boxplots...\n")

boxplot_genes <- de_df[significant == "Yes", gene]
if (length(boxplot_genes) > 8) boxplot_genes <- head(boxplot_genes, 8)
if (length(boxplot_genes) == 0) boxplot_genes <- head(de_df[order(p_value), gene], 6)

for (g in boxplot_genes) {
  t_expr <- expr_mat[g, tumor_cols]
  n_expr <- expr_mat[g, normal_cols]

  plot_df <- rbind(
    data.frame(Expression = t_expr, Group = "Tumor", stringsAsFactors = FALSE),
    data.frame(Expression = n_expr, Group = "Normal", stringsAsFactors = FALSE)
  )

  pv <- de_df[gene == g, p_value]
  fc <- de_df[gene == g, log2FC]

  p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.size = 0.5, width = 0.5, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
    scale_fill_manual(values = c("Normal" = "#3498DB", "Tumor" = "#E74C3C")) +
    labs(title = g, subtitle = sprintf("log2FC=%+.2f, p=%.2e", fc, pv),
         y = "log2(RSEM+1)") +
    theme_minimal(base_size = 11) + theme(legend.position = "none")

  ggsave(file.path(OUT_DIR, sprintf("FigX_boxplot_TvsN_%s.pdf", g)),
         p, width = 3.5, height = 4)
}

# ===== 6. Combined boxplot panel (coloc-passing genes) =====
coloc_found <- intersect(coloc_pass, found)
if (length(coloc_found) >= 2) {
  boxplots <- list()
  for (g in coloc_found) {
    t_expr <- expr_mat[g, tumor_cols]
    n_expr <- expr_mat[g, normal_cols]
    plot_df <- rbind(
      data.frame(Expression = t_expr, Group = "Tumor"),
      data.frame(Expression = n_expr, Group = "Normal")
    )
    pv <- de_df[gene == g, p_value]
    fc <- de_df[gene == g, log2FC]

    boxplots[[g]] <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(outlier.size = 0.3, width = 0.5, alpha = 0.85) +
      scale_fill_manual(values = c("Normal" = "#3498DB", "Tumor" = "#E74C3C")) +
      labs(title = sprintf("%s\nFC=%+.2f p=%.1e", g, fc, pv), y = "log2(RSEM+1)") +
      theme_minimal(base_size = 9) + theme(legend.position = "none")
  }

  combo <- plot_grid(plotlist = boxplots, ncol = min(3, length(boxplots)))
  ggsave(file.path(OUT_DIR, "FigX_boxplot_coloc_panel.pdf"),
         combo, width = 9, height = 3.2 * ceiling(length(boxplots)/3))
  cat("Combined coloc boxplot panel saved.\n")
}

# ===== 7. Volcano-style plot =====
cat("Generating volcano plot...\n")

volc <- copy(de_df)
volc[, neg_log10p := -log10(p_value)]
volc[, label := ifelse(significant == "Yes" | abs(log2FC) > 1 | neg_log10p > 2, gene, "")]
# Add coloc indicator
volc[, is_coloc := gene %in% coloc_pass]

p_volc <- ggplot(volc, aes(x = log2FC, y = neg_log10p)) +
  geom_point(aes(color = significant == "Yes", shape = is_coloc), size = 2.5, alpha = 0.8) +
  geom_text_repel(aes(label = label), size = 3.2, max.overlaps = 15, box.padding = 0.3) +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.4) +
  labs(title = "Key SMR Genes: Tumor vs Normal (TCGA COAD+READ)",
       subtitle = sprintf("%d tumors vs %d normals; triangle = coloc-passing",
                          length(tumor_cols), length(normal_cols)),
       x = "log2 Fold Change (Tumor/Normal)",
       y = "-log10(p-value)",
       shape = "Coloc-passing", color = "Significant (p<0.05)") +
  theme_minimal(base_size = 11)

ggsave(file.path(OUT_DIR, "FigX_volcano_TvsN.pdf"), p_volc, width = 8.5, height = 6)
cat("Volcano plot saved.\n")

# ===== 8. Paired line plots for coloc genes =====
if (length(paired_patients) >= 5) {
  cat("Generating paired line plots...\n")

  paired_plots <- list()
  for (g in head(coloc_found, 6)) {
    t_expr <- expr_mat[g, paired_tumor_cols]
    n_expr <- expr_mat[g, paired_normal_cols]
    pair_df <- data.frame(
      Patient = rep(1:length(t_expr), 2),
      Expression = c(t_expr, n_expr),
      Group = rep(c("Tumor", "Normal"), each = length(t_expr))
    )
    fc <- de_df[gene == g, log2FC]
    pv <- de_df[gene == g, p_value]

    paired_plots[[g]] <- ggplot(pair_df, aes(x = Group, y = Expression, group = Patient)) +
      geom_line(color = "gray70", alpha = 0.4, linewidth = 0.3) +
      geom_point(aes(color = Group), size = 1) +
      scale_color_manual(values = c("Normal" = "#3498DB", "Tumor" = "#E74C3C")) +
      labs(title = sprintf("%s (FC=%+.2f, p=%.3f)", g, fc, pv), y = "log2(RSEM+1)") +
      theme_minimal(base_size = 9) + theme(legend.position = "none")
  }

  combo <- plot_grid(plotlist = paired_plots,
                      ncol = min(3, length(paired_plots)))
  ggsave(file.path(OUT_DIR, "FigX_paired_TvsN_coloc.pdf"),
         combo, width = 9, height = 3 * ceiling(length(paired_plots)/3))
  cat("Paired line plots saved.\n")
}

cat(sprintf("\nAll outputs saved to: %s\n", OUT_DIR))
cat("Done.\n")
