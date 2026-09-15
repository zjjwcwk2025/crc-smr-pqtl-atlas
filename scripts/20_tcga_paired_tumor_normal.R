#!/usr/bin/env Rscript
# Phase 6d: TCGA COAD paired tumor-normal differential expression
# Downloads via TCGAbiolinks, analyses key SMR/coloc genes
# Output: DE boxplots, volcano data, summary table

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# ===== Config =====
OUT_DIR    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase6_scrna"
TCGA_DIR   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/data/tcga_coad"
COLOC_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_coloc_top10.csv"
GDC_CACHE  <- "/ifs1/User/zhouman/tmp/GDCdata"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TCGA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GDC_CACHE, showWarnings = FALSE, recursive = TRUE)

# ===== 1. Download TCGA-COAD RNA-seq (STAR counts) =====
cat("=== Step 1: Query GDC for TCGA-COAD RNA-seq ===\n")

query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open"
)

cat(sprintf("Query returned %d cases, %d files\n",
            length(unique(query$results[[1]]$cases)),
            nrow(query$results[[1]])))

# Download if not already cached
cat("Downloading...\n")
GDCdownload(query, directory = GDC_CACHE, method = "api")

# ===== 2. Prepare expression matrix =====
cat("\n=== Step 2: Prepare expression data ===\n")

exp_data <- GDCprepare(query, directory = GDC_CACHE, summarizedExperiment = TRUE)
cat(sprintf("Expression: %d genes x %d samples\n", nrow(exp_data), ncol(exp_data)))

# Extract counts
counts <- assay(exp_data, "unstranded")
gene_info <- rowData(exp_data)
sample_info <- colData(exp_data)

cat(sprintf("Sample types:\n"))
print(table(sample_info$sample_type))

# ===== 3. Identify paired tumor-normal samples =====
cat("\n=== Step 3: Identify paired samples ===\n")

# Patient ID = first 12 chars of TCGA barcode
sample_info$patient <- substr(sample_info$barcode, 1, 12)

# Find patients with both tumor (01) and normal (11)
tumor_samples <- sample_info[sample_info$sample_type == "Primary Tumor", ]
normal_samples <- sample_info[sample_info$sample_type == "Solid Tissue Normal", ]

paired_patients <- intersect(tumor_samples$patient, normal_samples$patient)
cat(sprintf("Paired tumor-normal patients: %d\n", length(paired_patients)))

if (length(paired_patients) < 5) {
  cat("Too few paired samples, using all available samples instead.\n")

  # Fallback: use all tumor + normal samples (unpaired comparison)
  tumor_idx <- which(sample_info$sample_type == "Primary Tumor")
  normal_idx <- which(sample_info$sample_type == "Solid Tissue Normal")

  cat(sprintf("Tumor samples: %d, Normal samples: %d\n",
              length(tumor_idx), length(normal_idx)))

  if (length(normal_idx) < 3) {
    cat("ERROR: Too few normal samples. Aborting.\n")
    quit(save = "no", status = 1)
  }
} else {
  # Use only paired samples
  tumor_idx <- which(sample_info$patient %in% paired_patients &
                      sample_info$sample_type == "Primary Tumor")
  normal_idx <- which(sample_info$patient %in% paired_patients &
                       sample_info$sample_type == "Solid Tissue Normal")

  # Ensure same order for paired analysis
  tumor_ordered <- sample_info$patient[tumor_idx]
  normal_ordered <- sample_info$patient[normal_idx]
  matched <- match(tumor_ordered, normal_ordered)
  normal_idx <- normal_idx[matched]
}

# ===== 4. Define genes =====
coloc <- fread(COLOC_FILE)
coloc_pass <- coloc[PPH4 > 0.8, gene]

target_genes <- unique(c(
  coloc_pass, "CDKN1A", "TNF", "LTA", "NCF2",
  "MMP11", "SMAD7", "CASC8", "POLD3", "BMP2",
  "GPBAR1", "TMBIM1", "LAMC1", "UTP11", "LIMA1", "PNKD"
))

# Map gene symbols to rowData
gene_names <- gene_info$gene_name
found <- intersect(target_genes, gene_names)
cat(sprintf("\nGenes found: %d / %d\n", length(found), length(target_genes)))

# ===== 5. Differential expression for each gene =====
cat("\n=== Step 4: Differential expression ===\n")

de_results <- list()

for (g in found) {
  gene_row <- which(gene_names == g)[1]
  if (is.na(gene_row)) next

  tumor_expr <- counts[gene_row, tumor_idx]
  normal_expr <- counts[gene_row, normal_idx]

  # Handle paired or unpaired
  if (length(paired_patients) >= 5) {
    # Paired: compute fold change per pair then test
    fc_values <- log2((tumor_expr + 1) / (normal_expr + 1))
    log2fc <- mean(fc_values, na.rm = TRUE)
    test <- wilcox.test(tumor_expr, normal_expr, paired = TRUE)
  } else {
    log2fc <- log2((mean(tumor_expr) + 1) / (mean(normal_expr) + 1))
    test <- wilcox.test(tumor_expr, normal_expr, paired = FALSE)
  }

  de_results[[g]] <- data.frame(
    gene = g,
    mean_tumor = mean(tumor_expr),
    mean_normal = mean(normal_expr),
    log2FC = log2fc,
    p_value = test$p.value,
    direction = ifelse(log2fc > 0, "Up in Tumor", "Down in Tumor"),
    significant = ifelse(test$p.value < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %-8s: TPM_tumor=%.1f TPM_normal=%.1f log2FC=%.2f p=%.2e %s\n",
              g, mean(tumor_expr), mean(normal_expr),
              log2fc, test$p.value,
              ifelse(test$p.value < 0.05, "*", "")))
}

# ===== 6. Output DE table =====
de_df <- rbindlist(de_results)
de_df <- de_df[order(-abs(log2FC))]

cat("\n=== Tumor vs Normal Differential Expression ===\n")
print(de_df[, .(gene, log2FC, p_value, direction, significant)])

fwrite(de_df, file.path(OUT_DIR, "TableS_tcga_tumor_vs_normal.csv"))

# ===== 7. Boxplots for top DEGs =====
cat("\nGenerating boxplots...\n")

top_de_genes <- de_df[significant == "Yes"]
if (nrow(top_de_genes) > 8) top_de_genes <- head(top_de_genes, 8)
if (nrow(top_de_genes) == 0) top_de_genes <- head(de_df[order(p_value)], 6)

for (g in top_de_genes$gene) {
  gene_row <- which(gene_names == g)[1]
  if (is.na(gene_row)) next

  tumor_expr <- counts[gene_row, tumor_idx]
  normal_expr <- counts[gene_row, normal_idx]

  plot_df <- rbind(
    data.frame(Expression = log2(tumor_expr + 1), Group = "Tumor", stringsAsFactors = FALSE),
    data.frame(Expression = log2(normal_expr + 1), Group = "Normal", stringsAsFactors = FALSE)
  )

  p_val <- de_df[gene == g, p_value]
  logfc <- de_df[gene == g, log2FC]

  p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.size = 0.8, width = 0.5) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
    scale_fill_manual(values = c("Normal" = "#3498DB", "Tumor" = "#E74C3C")) +
    labs(title = sprintf("%s: Tumor vs Normal", g),
         subtitle = sprintf("log2FC=%.2f, Wilcoxon p=%.2e", logfc, p_val),
         y = "log2(Count + 1)") +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave(file.path(OUT_DIR, sprintf("FigX_boxplot_%s_TvsN.pdf", g)),
         p, width = 4, height = 4.5)
}

# ===== 8. Volcano-style summary =====
cat("\nGenerating volcano-style summary plot...\n")

de_plot <- copy(de_df)
de_plot[, neg_log10p := -log10(p_value)]
de_plot[, label := ifelse(significant == "Yes" | abs(log2FC) > 1 | neg_log10p > 2, gene, "")]

p_volcano <- ggplot(de_plot, aes(x = log2FC, y = neg_log10p)) +
  geom_point(aes(color = significant == "Yes"), size = 2.5, alpha = 0.8) +
  geom_text_repel(aes(label = label), size = 3.5, max.overlaps = 15) +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  labs(title = "Key SMR Genes: Tumor vs Normal (TCGA COAD)",
       x = "log2 Fold Change (Tumor/Normal)",
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "FigX_volcano_TvsN.pdf"), p_volcano, width = 8, height = 6)

# ===== 9. Paired line plot for coloc-passing genes =====
if (length(paired_patients) >= 5) {
  cat("\nGenerating paired plots for coloc genes...\n")

  coloc_found <- intersect(coloc_pass, found)
  if (length(coloc_found) >= 2) {
    paired_plots <- list()
    for (g in head(coloc_found, 6)) {
      gene_row <- which(gene_names == g)[1]
      t_expr <- log2(counts[gene_row, tumor_idx] + 1)
      n_expr <- log2(counts[gene_row, normal_idx] + 1)

      pair_df <- data.frame(
        Patient = rep(1:length(t_expr), 2),
        Expression = c(t_expr, n_expr),
        Group = rep(c("Tumor", "Normal"), each = length(t_expr))
      )

      p_val <- de_df[gene == g, p_value]
      logfc <- de_df[gene == g, log2FC]

      paired_plots[[g]] <- ggplot(pair_df, aes(x = Group, y = Expression, group = Patient)) +
        geom_line(color = "gray70", alpha = 0.4) +
        geom_point(aes(color = Group), size = 1.2) +
        scale_color_manual(values = c("Normal" = "#3498DB", "Tumor" = "#E74C3C")) +
        labs(title = sprintf("%s (FC=%.2f, p=%.3f)", g, logfc, p_val),
             y = "log2(Count+1)") +
        theme_minimal() + theme(legend.position = "none")
    }

    combo <- cowplot::plot_grid(plotlist = paired_plots,
                                 ncol = min(3, length(paired_plots)))
    ggsave(file.path(OUT_DIR, "FigX_paired_TvsN_coloc.pdf"),
           combo, width = 10, height = 3.5 * ceiling(length(paired_plots)/3))
  }
}

cat(sprintf("\nAll outputs saved to: %s\n", OUT_DIR))
cat("Done.\n")
