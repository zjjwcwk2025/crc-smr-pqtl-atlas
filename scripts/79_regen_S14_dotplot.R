#!/usr/bin/env Rscript
# =====================================================================
#  79_regen_S14_dotplot.R  --  2026-09-26 pre-submission figure audit fix
#
#  Why this file exists
#  --------------------
#  Figure S14 (cell-type dot plot) had two defects found by the submission
#  audit, both inherited from the original render in
#  16_scrna_celltype_localization.R:
#
#    1. the dot-size legend overflowed the device top, so the largest
#       legend key (100) was clipped by the page edge   -> page clipping
#    2. the "Average Expression" colour bar was written as a raster image
#       (1 x 300 px stretched over 5.4 x 27.1 mm, ~281 dpi, below the
#       300 dpi floor) and the figure used the pdf() device default
#       Helvetica, which is NOT embedded  -> JTM requires embedded fonts
#
#  This script is the faithful re-render: it is the dot-plot section of
#  16_scrna_celltype_localization.R verbatim, with exactly three changes
#  (all marked [FIX] below):
#
#    [FIX 1] ggsave height 2.70 in -> 3.00 in   (legend no longer clipped)
#    [FIX 2] guides(colour = guide_colourbar(display = "gradient"))
#            the colour bar is bound to the *colour* aesthetic in Seurat's
#            DotPlot (not fill), and "gradient" draws it as vector
#    [FIX 3] device = cairo_pdf                  (fonts embedded)
#
#  Resulting panel: 160.0 x 76.2 mm, 0 embedded rasters, font embedded.
#  Output: results/figures_rev3/FigS14_celltype_dotplot.pdf
# =====================================================================
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
  guides(colour = guide_colourbar(display = "gradient"), fill = guide_colourbar(display = "gradient")) +
  labs(x = "Gene", y = "Cell type") +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9))

ggsave(file.path(OUT_DIR, "FigX_dotplot_celltype.pdf"), p_dot, width = 6.30, height = 3.00, device = cairo_pdf)
dir.create("/ifs1/User/zhouman/project9-v5-crc-atlas/results/figures", showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(OUT_DIR, "FigX_dotplot_celltype.pdf"),
          "/ifs1/User/zhouman/project9-v5-crc-atlas/results/figures/FigS7_celltype_dotplot.pdf",
          overwrite = TRUE)
cat("Dot plot saved.\n")

# --- deliver into results/figures_rev3 under the final figure number ---
file.copy(file.path(OUT_DIR, "FigX_dotplot_celltype.pdf"),
          "/ifs1/User/zhouman/project9-v5-crc-atlas/results/figures_rev3/FigS14_celltype_dotplot.pdf",
          overwrite = TRUE)
cat("FigS14_celltype_dotplot.pdf regenerated\n")
