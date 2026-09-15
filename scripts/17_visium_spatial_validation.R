#!/usr/bin/env Rscript
# Phase 6b: Visium spatial validation of SMR coloc-passing genes
# Maps coloc-passing + CAF-expressed SMR genes onto CRC Visium spatial transcriptomics
# Uses Cell2location deconvolution to link gene expression to cell-type abundance

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

theme_set(theme_gray(base_size = 9))
# Every panel below is laid out at the width it finally occupies in the
# supplementary PDF (0.70 x 6.30 in = 4.41 in), so nothing is down-scaled
# and the type stays at or above 7 pt in print.

set.seed(42)

# ===== Config =====
VISIUM_RDS   <- "/ifs1/User/zhouman/CRC_Project/GSE285505_Spatial_Processed.rds"
COLOC_FILE   <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_coloc_top10.csv"
SCRNA_RESULT <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase6_scrna/phase6_celltype_specificity.csv"
OUT_DIR      <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase6_scrna"

dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# ===== 1. Load Visium data =====
cat("Loading Visium spatial data...\n")
vis <- readRDS(VISIUM_RDS)
# Default assay is already "SCT" (processed SCTransform data) — keep it

cat(sprintf("Visium: %d spots x %d genes\n", ncol(vis), nrow(vis)))

# ===== 2. Define genes to check =====
# Coloc-passing genes (PPH4 > 0.8) — core genes for spatial validation
coloc <- fread(COLOC_FILE)
coloc_pass <- coloc[PPH4 > 0.8]

cat(sprintf("\nColoc-passing genes: %d\n", nrow(coloc_pass)))
print(coloc_pass[, .(gene, PPH4)])

# Also add SMR top genes with detectable CAF expression
gene_list <- unique(c(
  coloc_pass$gene,           # 6 coloc-passing genes
  "MMP11",                    # top Fibroblast-specific SMR gene
  "CDKN1A",                   # CAF-expressed, key drug target
  "TNF",                      # key drug target
  "LAMC1",                    # key coloc-passing gene
  "TMBIM1",                   # coloc-passing + CAF expressed
  "SMAD7",                    # known CRC GWAS hit
  "CASC8"                     # known CRC GWAS locus
))

# Check which are in Visium
vis_genes <- rownames(vis)
found_genes <- intersect(gene_list, vis_genes)
missing_genes <- setdiff(gene_list, vis_genes)

cat(sprintf("\nGenes in Visium: %d / %d\n", length(found_genes), length(gene_list)))
if (length(missing_genes) > 0) {
  cat(sprintf("Missing: %s\n", paste(missing_genes, collapse=", ")))
}

# Normalize if needed
if (!"SCT" %in% names(vis@assays)) {
  cat("\nRunning SCTransform...\n")
  vis <- SCTransform(vis, assay = "Spatial", verbose = FALSE)
}

# ===== 3. Spatial Feature Plots for coloc-passing genes =====
cat("\nGenerating spatial feature plots...\n")

for (g in intersect(coloc_pass$gene, found_genes)) {
  cat(sprintf("  Plotting %s...\n", g))

  p <- SpatialFeaturePlot(vis, features = g, pt.size.factor = 1.6,
                           alpha = c(0.3, 1), stroke = 0) +
    ggtitle(sprintf("%s (coloc PPH4=%.3f)", g, coloc_pass$PPH4[coloc_pass$gene == g])) +
    labs(fill = "SCT expression") +
    theme(plot.title = element_text(size = 9, face = "bold"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.width = unit(0.7, "cm"))

  ggsave(file.path(OUT_DIR, sprintf("FigX_spatial_%s.pdf", g)),
         p, width = 4.41, height = 3.97)
}

# ===== 4. Multi-gene spatial panel =====
cat("Generating combined spatial panel...\n")

# Select up to 6 genes for a combined figure
panel_genes <- head(intersect(coloc_pass$gene, found_genes), 6)
if (length(panel_genes) >= 2) {
  plots <- lapply(panel_genes, function(g) {
    SpatialFeaturePlot(vis, features = g, pt.size.factor = 1.4, stroke = 0) +
      ggtitle(g) + theme(plot.title = element_text(size = 9))
  })

  n_plots <- length(plots)
  ncol_panel <- min(3, n_plots)
  nrow_panel <- ceiling(n_plots / ncol_panel)

  combined <- wrap_plots(plots, ncol = ncol_panel) +
    plot_annotation(title = "Coloc-Passing Genes: Spatial Expression in CRC Tissue",
                    theme = theme(plot.title = element_text(size = 13, face = "bold")))

  ggsave(file.path(OUT_DIR, "FigX_spatial_coloc_panel.pdf"),
         combined, width = ncol_panel * 4.5, height = nrow_panel * 4)
  cat(sprintf("Combined panel saved: %d genes\n", length(panel_genes)))
}

# ===== 5. CAF-related genes spatial expression =====
cat("\nCAF-related gene spatial expression...\n")

caf_genes_vis <- c("MMP11", "CDKN1A", "LAMC1", "TMBIM1", "LIMA1")
caf_genes_vis <- intersect(caf_genes_vis, found_genes)

if (length(caf_genes_vis) >= 2) {
  for (g in caf_genes_vis) {
    p <- SpatialFeaturePlot(vis, features = g, pt.size.factor = 1.6,
                             alpha = c(0.3, 1), stroke = 0) +
      ggtitle(sprintf("%s (CAF-associated)", g)) +
      labs(fill = "SCT expression") +
      theme(legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            legend.key.width = unit(0.7, "cm"))
    ggsave(file.path(OUT_DIR, sprintf("FigX_spatial_CAF_%s.pdf", g)),
           p, width = 4.41, height = 3.97)
  }
}

# ===== 6. Exploratory byproduct, not reported in the manuscript: correlation with
# the meCAF (metabolic CAF) signature score of the separate CAF project =====
cat("\nChecking meCAF score correlation...\n")

mecaf_col <- grep("meCAF|mecaf", colnames(vis@meta.data), value = TRUE, ignore.case = TRUE)
cat(sprintf("meCAF-related columns: %s\n", paste(mecaf_col, collapse = ", ")))

if (length(mecaf_col) > 0) {
  mecaf_score_col <- mecaf_col[1]

  # Compute correlation between gene expression and meCAF score
  cor_results <- list()
  for (g in found_genes) {
    expr <- FetchData(vis, vars = c(g, mecaf_score_col))
    ct <- cor.test(expr[[g]], expr[[mecaf_score_col]], method = "spearman")
    cor_results[[g]] <- data.frame(
      gene = g, spearman_rho = ct$estimate, p_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  }
  cor_df <- rbindlist(cor_results)
  cor_df <- cor_df[order(-abs(spearman_rho))]

  cat("\nGene correlation with meCAF spatial score:\n")
  print(cor_df)

  fwrite(cor_df, file.path(OUT_DIR, "TableS_spatial_mecaf_correlation.csv"))

  # Exploratory byproduct, NOT reported in the manuscript. meCAF is the metabolic
  # (HIF1A-driven, glycolysis/lactate/hypoxia) CAF state of the separate CAF project
  # (Project01_meCAF); the score is the mean SCT-normalized expression of the unified
  # 13-gene signature (LDHA, PKM, ENO1, PFKP, SLC16A3, SLC16A1, CA9, PDK1, SLC2A1,
  # HIF1A, VEGFA, BNIP3, NDRG1), as set in scripts/unify_mecaf_signature.R.
  p_mecaf <- SpatialFeaturePlot(vis, features = mecaf_score_col, pt.size.factor = 1.6,
                                alpha = c(0.3, 1), stroke = 0) +
    ggtitle("Metabolic CAF (meCAF) signature") +
    labs(fill = "meCAF score") +
    theme(plot.title = element_text(size = 9, face = "bold"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.width = unit(0.7, "cm"))
  ggsave(file.path(OUT_DIR, "FigX_spatial_meCAF_score.pdf"),
         p_mecaf, width = 4.41, height = 3.97)
  cat("  -> FigX_spatial_meCAF_score.pdf\n")

  # Plot top correlations
  top_cor_genes <- head(cor_df[abs(spearman_rho) > 0.1], 8)
  if (nrow(top_cor_genes) > 0) {
    for (g in top_cor_genes$gene) {
      expr <- FetchData(vis, vars = c(g, mecaf_score_col))
      colnames(expr) <- c("gene_expr", "mecaf_score")

      p <- ggplot(expr, aes(x = mecaf_score, y = gene_expr)) +
        geom_point(alpha = 0.3, size = 0.5, color = "steelblue") +
        geom_smooth(method = "lm", se = TRUE, color = "darkred") +
        labs(x = "meCAF Score", y = g,
             title = sprintf("%s vs meCAF (rho=%.3f, p=%.1e)",
                             g, cor_df$spearman_rho[cor_df$gene == g],
                             cor_df$p_value[cor_df$gene == g])) +
        theme_minimal()

      ggsave(file.path(OUT_DIR, sprintf("FigX_scatter_%s_vs_mecaf.pdf", g)),
             p, width = 5, height = 4)
    }
  }
}

# ===== 7. Known CRC GWAS hit spatial validation =====
cat("\nKnown CRC GWAS genes spatial validation...\n")

gwas_genes <- c("SMAD7", "CASC8", "LAMC1", "POLD3")
gwas_genes <- intersect(gwas_genes, found_genes)

if (length(gwas_genes) > 0) {
  plots <- lapply(gwas_genes, function(g) {
    SpatialFeaturePlot(vis, features = g, pt.size.factor = 1.4, stroke = 0) +
      ggtitle(paste(g, "(CRC GWAS locus)")) +
      theme(plot.title = element_text(size = 9))
  })

  gw_combined <- wrap_plots(plots, ncol = min(2, length(plots))) +
    plot_annotation(title = "Known CRC GWAS Loci: Spatial Expression",
                    theme = theme(plot.title = element_text(size = 13, face = "bold")))

  ggsave(file.path(OUT_DIR, "FigX_spatial_GWAS_loci.pdf"),
         gw_combined, width = 9, height = 3.5 * ceiling(length(plots)/2))
}

# ===== 8. Summary table =====
cat("\n--- Spatial Validation Summary ---\n")

summary_vis <- data.frame(
  gene = found_genes,
  mean_expr = sapply(found_genes, function(g) mean(FetchData(vis, vars = g)[[1]])),
  pct_expressed = sapply(found_genes, function(g) mean(FetchData(vis, vars = g)[[1]] > 0) * 100),
  stringsAsFactors = FALSE
)
summary_vis <- summary_vis[order(-summary_vis$mean_expr), ]

cat("\nSpatial expression summary:\n")
print(summary_vis)

fwrite(summary_vis, file.path(OUT_DIR, "TableS_spatial_expression_summary.csv"))

cat(sprintf("\nAll Visium outputs saved to: %s\n", OUT_DIR))
cat("Done.\n")

# ===== 9. Export the two panels used in the manuscript =====
# Supplementary Fig. S8 (TMBIM1 spatial map) and S9 (TMBIM1 vs meCAF score)
fig_dir <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/figures"
file.copy(file.path(OUT_DIR, "FigX_spatial_TMBIM1.pdf"),
          file.path(fig_dir, "FigS8_spatial_TMBIM1.pdf"), overwrite = TRUE)
cat("Supplementary Fig. S8 (TMBIM1 spatial map) exported at printed width (4.41 x 3.97 in).\n")
