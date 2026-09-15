#!/usr/bin/env Rscript
# Phase 6: scRNA cell-type localization of 75 SMR genes
# Maps SMR-significant genes to CRC single-cell atlas cell types
# Focus on CAF subpopulations (meCAF, iCAF, myCAF) as drug-target niches

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

set.seed(42)

# ===== Config =====
SMR_FILE     <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase1_bonferroni_significant.tsv"
SCRNA_RDS    <- "/ifs1/User/zhouman/CRC-Analysis-Final/run_genesis_v13_Global_Integration_20260515_134302/13_CRC_Global_Integrated_Fine_20260515_134302.rds"
CAF_RDS      <- "/ifs1/User/zhouman/CRC-Analysis-Final/run_genesis_v11_CAF_Annotated_20260515_125949/11_CRC_CAF_Final_Annotated_20260515_125949.rds"
OUT_DIR      <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase6_scrna"

dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# ===== 1. Map Ensembl → gene symbols =====
smr <- fread(SMR_FILE)
ensg_ids <- unique(smr$Gene)
cat(sprintf("Total unique Ensembl IDs: %d\n", length(ensg_ids)))

# Map Ensembl to gene symbol
map_res <- select(org.Hs.eg.db, keys = ensg_ids,
                  columns = c("SYMBOL"), keytype = "ENSEMBL")
# Keep one symbol per Ensembl ID (first match)
map_res <- map_res[!duplicated(map_res$ENSEMBL), ]

cat(sprintf("Mapped to symbols: %d / %d\n", sum(!is.na(map_res$SYMBOL)), length(ensg_ids)))

# Merge with SMR data
smr_mapped <- merge(smr, map_res, by.x = "Gene", by.y = "ENSEMBL", all.x = TRUE)
# Use probeID as fallback for unmapped
smr_mapped[, gene_symbol := ifelse(is.na(SYMBOL), probeID, SYMBOL)]

cat(sprintf("Genes with no symbol mapping: %d\n", sum(is.na(smr_mapped$SYMBOL))))
cat("Unmapped:\n")
print(smr_mapped[is.na(SYMBOL), unique(gene_symbol)])

# ===== 2. Load scRNA data =====
cat("\nLoading scRNA atlas...\n")
scrna <- readRDS(SCRNA_RDS)
DefaultAssay(scrna) <- "RNA"

# Normalize if not already
if (!"data" %in% names(scrna@assays$RNA)) {
  scrna <- NormalizeData(scrna, verbose = FALSE)
}

# Get gene overlap
sc_genes <- rownames(scrna)
smr_symbols <- unique(smr_mapped$gene_symbol)
overlap <- intersect(smr_symbols, sc_genes)
cat(sprintf("\nSMR genes found in scRNA: %d / %d (%.1f%%)\n",
            length(overlap), length(smr_symbols), 100*length(overlap)/length(smr_symbols)))

missing <- setdiff(smr_symbols, sc_genes)
if (length(missing) > 0) {
  cat(sprintf("Missing genes (%d):\n", length(missing)))
  print(missing)
}

# Keep only top 50 by SMR significance for visualization
smr_sig <- smr_mapped[order(p_SMR)]
smr_sig <- smr_sig[gene_symbol %in% overlap]
top_genes <- head(smr_sig$gene_symbol, 50)
cat(sprintf("\nTop 50 SMR genes for visualization: %d in scRNA\n", length(top_genes)))

# ===== 3. Cell-type average expression =====
cat("\nComputing cell-type expression profiles...\n")

# Use CellType_Major for the main analysis
Idents(scrna) <- "CellType_Major"
cell_types <- levels(Idents(scrna))

# Average expression per cell type
avg_expr <- AverageExpression(scrna, features = top_genes, assays = "RNA",
                               group.by = "CellType_Major", slot = "data")$RNA
avg_expr <- as.data.frame(avg_expr)
avg_expr$gene <- rownames(avg_expr)

# Percent expressed
pct_expr <- matrix(NA, nrow = length(top_genes), ncol = length(cell_types))
rownames(pct_expr) <- top_genes
colnames(pct_expr) <- cell_types
for (ct in cell_types) {
  cells_ct <- WhichCells(scrna, idents = ct)
  expr_ct <- GetAssayData(scrna, assay = "RNA", layer = "data")[top_genes, cells_ct, drop = FALSE]
  pct_expr[, ct] <- rowMeans(expr_ct > 0) * 100
}

cat("Average expression computed.\n")

# ===== 4. Cell-type specificity score =====
# For each gene, compute: max_expression / mean_expression across cell types
# Higher score = more cell-type specific
spec_score <- apply(avg_expr[, cell_types, drop = FALSE], 1, function(x) {
  if (max(x) == 0) return(0)
  max(x) / (mean(x[x > 0]) + 1e-10)
})
spec_df <- data.frame(
  gene = names(spec_score),
  specificity = spec_score,
  top_celltype = apply(avg_expr[, cell_types, drop = FALSE], 1, function(x) names(which.max(x))),
  top_expr = apply(avg_expr[, cell_types, drop = FALSE], 1, max),
  stringsAsFactors = FALSE
)
spec_df <- spec_df[order(-spec_df$specificity), ]
cat("\nTop 15 cell-type-specific genes:\n")
print(head(spec_df, 15))

# ===== 5. Dot Plot =====
cat("\nGenerating dot plot...\n")

# Focus on genes enriched in Fibroblasts/CAFs
fibro_genes <- spec_df$gene[spec_df$top_celltype == "Fibroblasts"]
other_genes <- setdiff(top_genes, fibro_genes)

# Select ~30 genes for readable dot plot: top 15 fibro-specific + top 15 others
plot_genes <- c(
  head(fibro_genes, 15),
  head(setdiff(spec_df$gene, fibro_genes), 15)
)
plot_genes <- intersect(plot_genes, top_genes)

# Generate dot plot
p_dot <- DotPlot(scrna, features = plot_genes, group.by = "CellType_Major",
                  assay = "RNA", dot.scale = 8) +
  RotatedAxis() +
  # Short title: the gene / cell / cell-type counts live in the figure caption.
  ggtitle("Cell-type expression of SMR genes") +
  labs(x = "Gene", y = "Cell type") +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9))

ggsave(file.path(OUT_DIR, "FigX_dotplot_celltype.pdf"), p_dot, width = 6.30, height = 2.70)
dir.create("/ifs1/User/zhouman/project9-v5-crc-atlas/results/figures", showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(OUT_DIR, "FigX_dotplot_celltype.pdf"),
          "/ifs1/User/zhouman/project9-v5-crc-atlas/results/figures/FigS7_celltype_dotplot.pdf",
          overwrite = TRUE)
cat("Dot plot saved.\n")

# ===== 6. CAF-specific Expression Analysis =====
cat("\n--- CAF Subpopulation Analysis ---\n")

# Load CAF object
caf <- readRDS(CAF_RDS)
DefaultAssay(caf) <- "RNA"

# Check CAF subtype column
caf_cols <- colnames(caf@meta.data)
cat("CAF metadata columns:", paste(grep("subtype|Subtype|cluster|Cluster|caf|CAF", caf_cols, value = TRUE), collapse = ", "), "\n")

# Identify CAF subtype column
caf_subtype_col <- grep("CAF_Subtype|caf_subtype", caf_cols, value = TRUE, ignore.case = TRUE)
if (length(caf_subtype_col) == 0) {
  # Try other columns
  caf_subtype_col <- grep("subtype", caf_cols, value = TRUE, ignore.case = TRUE)[1]
}
cat(sprintf("Using CAF subtype column: %s\n", caf_subtype_col))

caf_subtypes <- unique(caf@meta.data[[caf_subtype_col]])
cat(sprintf("CAF subtypes: %s\n", paste(caf_subtypes, collapse = ", ")))

# Which SMR genes are expressed in CAFs?
caf_genes <- intersect(top_genes, rownames(caf))
cat(sprintf("SMR genes in CAF data: %d\n", length(caf_genes)))

if (length(caf_genes) >= 3) {
  # Average expression in CAF subtypes
  Idents(caf) <- caf_subtype_col
  caf_avg <- AverageExpression(caf, features = caf_genes, assays = "RNA",
                                group.by = caf_subtype_col, slot = "data")$RNA

  # Heatmap of CAF-specific expression
  # Scale within genes
  caf_avg_scaled <- t(scale(t(as.matrix(caf_avg))))

  # Save CAF expression data
  caf_out <- as.data.frame(caf_avg)
  caf_out$gene <- rownames(caf_out)
  caf_out <- caf_out[order(-caf_out[["meCAF"]]), ]  # sort by meCAF expression
  fwrite(caf_out, file.path(OUT_DIR, "TableS_caf_subtype_expression.csv"))
  cat("CAF subtype expression table saved.\n")

  # Print top genes by CAF subtype
  cat("\nTop genes in meCAF:\n")
  meCAF_genes <- caf_out$gene[caf_out$meCAF > 0]
  if (length(meCAF_genes) > 0) print(head(caf_out[caf_out$gene %in% meCAF_genes, c("gene", "meCAF", "iCAF", "myCAF")], 10))
  else cat("  None with detectable expression\n")

  cat("\nTop genes in iCAF:\n")
  iCAF_genes <- caf_out$gene[caf_out$iCAF > 0]
  if (length(iCAF_genes) > 0) print(head(caf_out[caf_out$gene %in% iCAF_genes, c("gene", "meCAF", "iCAF", "myCAF")], 10))
  else cat("  None with detectable expression\n")

  cat("\nTop genes in myCAF:\n")
  myCAF_genes <- caf_out$gene[caf_out$myCAF > 0]
  if (length(myCAF_genes) > 0) print(head(caf_out[caf_out$gene %in% myCAF_genes, c("gene", "meCAF", "iCAF", "myCAF")], 10))
  else cat("  None with detectable expression\n")
}

# ===== 7. Heatmap: Top 30 SMR genes across cell types =====
cat("\nGenerating heatmap...\n")

top30 <- head(top_genes, 30)

# Get average expression
ht_data <- avg_expr[top30, cell_types, drop = FALSE]
ht_scaled <- t(scale(t(as.matrix(ht_data))))
ht_scaled[is.na(ht_scaled)] <- 0

# Clip for visualization
ht_scaled[ht_scaled > 2.5] <- 2.5
ht_scaled[ht_scaled < -2.5] <- -2.5

# Simple heatmap using pheatmap or base
if (requireNamespace("pheatmap", quietly = TRUE)) {
  library(pheatmap)

  # Add SMR p-value annotation
  smr_pvals <- smr_mapped[match(rownames(ht_scaled), smr_mapped$gene_symbol), ]
  anno_row <- data.frame(
    SMR_log10P = -log10(pmin(smr_pvals$p_SMR, 1e-50)),
    row.names = rownames(ht_scaled)
  )

  pdf(file.path(OUT_DIR, "FigX_heatmap_top30_celltype.pdf"), width = 8, height = 10)
  pheatmap(ht_scaled,
           main = "Top 30 SMR Significant Genes: Cell-Type Expression",
           cluster_rows = TRUE, cluster_cols = TRUE,
           annotation_row = anno_row,
           show_rownames = TRUE, fontsize_row = 7,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()
  cat("Heatmap saved.\n")
} else {
  cat("pheatmap not available, skipping heatmap.\n")
}

# ===== 8. Summary table =====
cat("\nGenerating summary...\n")

summary_tab <- smr_mapped[gene_symbol %in% overlap]
summary_tab <- summary_tab[order(p_SMR)]

# Add cell type info
summary_tab$top_celltype <- spec_df$top_celltype[match(summary_tab$gene_symbol, spec_df$gene)]
summary_tab$specificity <- spec_df$specificity[match(summary_tab$gene_symbol, spec_df$gene)]

# Add average expression in Fibroblasts
if ("Fibroblasts" %in% cell_types) {
  summary_tab$expr_fibro <- avg_expr$Fibroblasts[match(summary_tab$gene_symbol, avg_expr$gene)]
  summary_tab$expr_epithelial <- avg_expr$`Epithelial cells`[match(summary_tab$gene_symbol, avg_expr$gene)]
  summary_tab$expr_immune <- rowMeans(cbind(
    avg_expr$`T/NK cells`[match(summary_tab$gene_symbol, avg_expr$gene)],
    avg_expr$`B/Plasma cells`[match(summary_tab$gene_symbol, avg_expr$gene)],
    avg_expr$`Myeloid cells`[match(summary_tab$gene_symbol, avg_expr$gene)]
  ), na.rm = TRUE)
} else {
  summary_tab$expr_fibro <- NA
  summary_tab$expr_epithelial <- NA
  summary_tab$expr_immune <- NA
}

# Select columns
out_cols <- intersect(c("gene_symbol","Gene","p_SMR","b_SMR","top_celltype","specificity",
                         "expr_fibro","expr_epithelial","expr_immune"), names(summary_tab))
fwrite(summary_tab[, out_cols, with = FALSE], file.path(OUT_DIR, "phase6_celltype_specificity.csv"))
cat(sprintf("Summary table saved: %d genes\n", nrow(summary_tab)))

# ===== 9. Print key findings =====
cat("\n=== Phase 6 Key Findings ===\n")

# Genes with highest fibroblast expression
fibro_specific <- summary_tab[top_celltype == "Fibroblasts"][order(-expr_fibro)]
if (nrow(fibro_specific) > 0) {
  cat(sprintf("\n%d genes with highest Fibroblast expression:\n", min(10, nrow(fibro_specific))))
  print(head(fibro_specific[, .(gene_symbol, expr_fibro, specificity, p_SMR)], 10))
} else {
  cat("\nNo SMR genes show top expression in Fibroblasts.\n")
  # Show genes with any fibroblast expression
  has_fibro <- summary_tab[expr_fibro > 0.1][order(-expr_fibro)]
  if (nrow(has_fibro) > 0) {
    cat(sprintf("Genes with detectable fibroblast expression (>0.1): %d\n", nrow(has_fibro)))
    print(head(has_fibro[, .(gene_symbol, expr_fibro, top_celltype, specificity, p_SMR)], 10))
  }
}

# Top 5 cell-type specific genes overall
cat("\nTop 10 most cell-type-specific SMR genes:\n")
print(spec_df[spec_df$gene %in% overlap, ][1:min(10, sum(spec_df$gene %in% overlap)), ])

# Proportion by top cell type
cat("\nGenes by top cell type:\n")
print(table(summary_tab$top_celltype))

cat(sprintf("\nAll outputs saved to: %s\n", OUT_DIR))
cat("Done.\n")
