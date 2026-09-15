#!/usr/bin/env Rscript
# 思路 1: 空间微环境分区分析
# Uses scRNA-derived cell-type signature scoring (Seurat AddModuleScore)
# to classify Visium spots into microenvironment types
# Then computes Tier 1 gene expression per microenvironment zone

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

# ===== Config =====
SCRNA_RDS  <- "/ifs1/User/zhouman/CRC-Analysis-Final/run_genesis_v13_Global_Integration_20260515_134302/13_CRC_Global_Integrated_Fine_20260515_134302.rds"
VISIUM_RDS <- "/ifs1/User/zhouman/CRC_Project/GSE285505_Spatial_Processed.rds"
OUT_DIR    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase6_scrna"

# Tier 1 genes defined by coloc (PPH4 > 0.8) + SuSiE fine-mapping
MIN_MARKERS <- 5L
TOP_MARKERS <- 30L

# Tier 1 genes — read from coloc results (data-driven)
COLOC_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_coloc_all75.csv"
if (file.exists(COLOC_FILE)) {
  coloc <- fread(COLOC_FILE)
  # Coloc file uses 'gene' column
  tier1_raw <- coloc[PPH4 > 0.8 & status == "PASS", unique(gene)]
  # Exclude non-coding/lncRNA
  coding_regex <- "^(RP|AC|AL|AP|LINC|MIR|SNOR|[A-Z]+[0-9]+\\.[0-9]+)"
  tier1_genes <- tier1_raw[!grepl(coding_regex, tier1_raw)]
} else {
  # Fallback: known 6 Tier 1 from prior analysis
  tier1_genes <- c("UTP11", "BMP2", "PNKD", "LAMC1", "TMBIM1", "GPBAR1")
}
cat(sprintf("Tier 1 genes: %s\n", paste(tier1_genes, collapse=", ")))

# ===== 1. Load scRNA reference and extract cell-type markers =====
cat("Loading scRNA atlas...\n")
scrna <- readRDS(SCRNA_RDS)
cat(sprintf("scRNA: %d cells x %d genes\n", ncol(scrna), nrow(scrna)))

# Find cell-type column
celltype_col <- grep("CellType|celltype|cell_type|cluster", colnames(scrna@meta.data),
                     value=TRUE, ignore.case=TRUE)
if (length(celltype_col) == 0) {
  cell_types <- Idents(scrna)
  cat("Using Seurat active ident as cell-type labels\n")
  type_labels <- levels(cell_types)
} else {
  celltype_col <- celltype_col[1]
  type_labels <- unique(scrna@meta.data[[celltype_col]])
  cat(sprintf("Using metadata column '%s' with %d cell types\n", celltype_col, length(type_labels)))
  Idents(scrna) <- celltype_col
}

cat(sprintf("Cell types detected: %s\n", paste(type_labels, collapse=", ")))

# ===== 2. Get cell-type marker genes from scRNA =====
cat("\nComputing cell-type marker genes...\n")

assay_use <- if ("SCT" %in% names(scrna@assays)) "SCT" else "RNA"
DefaultAssay(scrna) <- assay_use

marker_list <- list()

for (ct in type_labels) {
  tryCatch({
    markers <- FindMarkers(scrna, ident.1=ct, group.by=celltype_col,
                           only.pos=TRUE, logfc.threshold=0.5, min.pct=0.2,
                           max.cells.per.ident=500)
    top_markers <- head(rownames(markers[order(-markers$avg_log2FC), ]), TOP_MARKERS)
    if (length(top_markers) >= MIN_MARKERS) {
      marker_list[[ct]] <- top_markers
      cat(sprintf("  %s: %d markers (top: %s)\n", ct, length(top_markers),
                  paste(head(top_markers, 3), collapse=", ")))
    } else {
      cat(sprintf("  %s: only %d markers (<%d), skipping\n", ct, length(top_markers), MIN_MARKERS))
    }
  }, error=function(e) {
    cat(sprintf("  %s: ERROR - %s\n", ct, e$message))
  })
}

safe_names <- make.names(names(marker_list))
cat(sprintf("\nTotal cell types with sufficient markers: %d\n", length(marker_list)))

# ===== 3. Load Visium and score each spot for each cell type =====
cat("\nLoading Visium data...\n")
vis <- readRDS(VISIUM_RDS)
cat(sprintf("Visium: %d spots x %d genes\n", ncol(vis), nrow(vis)))

if (!"SCT" %in% names(vis@assays)) {
  cat("Running SCTransform...\n")
  vis <- SCTransform(vis, assay="Spatial", verbose=FALSE)
}

cat("\nScoring Visium spots for cell-type signatures...\n")

for (i in seq_along(marker_list)) {
  ct_name <- names(marker_list)[i]
  safe_name <- safe_names[i]
  ct_markers <- intersect(marker_list[[ct_name]], rownames(vis))
  if (length(ct_markers) < MIN_MARKERS) {
    cat(sprintf("  %s: only %d markers in Visium, skipping\n", ct_name, length(ct_markers)))
    next
  }
  cat(sprintf("  %s: scoring %d markers...\n", ct_name, length(ct_markers)))
  vis <- AddModuleScore(vis, features=list(ct_markers), name=safe_name, ctrl=100)
}

# ===== 4. Collect cell-type scores =====
cat("\nCollecting cell-type scores...\n")

score_cols <- grep("^[A-Za-z].*[0-9]$", colnames(vis@meta.data), value=TRUE)
score_cols <- score_cols[grepl("1$", score_cols)]

ct_score_cols <- list()
for (i in seq_along(marker_list)) {
  safe_name <- safe_names[i]
  matching <- grep(paste0("^", safe_name), score_cols, value=TRUE)
  if (length(matching) > 0) {
    ct_score_cols[[names(marker_list)[i]]] <- matching[1]
  }
}

cat(sprintf("Collected scores for %d / %d cell types\n", length(ct_score_cols), length(marker_list)))

# ===== 5. Data-driven selection of k (elbow method) =====
cat("\nSelecting optimal number of microenvironment clusters...\n")

score_mat <- as.data.frame(vis@meta.data[, unlist(ct_score_cols), drop=FALSE])
colnames(score_mat) <- names(ct_score_cols)
score_mat <- scale(score_mat)

# Scan k=2 to min(8, n_celltypes*2) using within-cluster sum of squares
n_features <- length(ct_score_cols)
k_max <- min(8, n_features * 2)
k_range <- 2:k_max

wss <- sapply(k_range, function(k) {
  km_try <- kmeans(score_mat, centers=k, nstart=25)
  km_try$tot.withinss
})

# Find elbow: point of maximum curvature (second derivative)
d1 <- diff(wss)
d2 <- diff(d1)
# elbow = k with max second derivative + 2 (offset for k=2 start + diff)
elbow_idx <- which.min(d2) + 2
opt_k <- k_range[elbow_idx]

cat(sprintf("WSS by k: %s\n", paste(sprintf("k=%d:%.1f", k_range, wss), collapse=", ")))
cat(sprintf("Optimal k (elbow method) = %d\n", opt_k))

# Generate elbow plot
elbow_df <- data.frame(k=k_range, wss=wss)
# Rendered at the final printed width (0.56 x 6.30 in = 3.53 in) with a 9 pt
# base size, so the type is not down-scaled below 7 pt in print. Panel titles
# are omitted because the figure caption carries that information.
elbow_plot <- ggplot(elbow_df, aes(x=k, y=wss)) +
  geom_point(size=1.6) +
  geom_line(linewidth=0.5) +
  geom_vline(xintercept=opt_k, linetype="dashed", color="darkred", linewidth=0.5) +
  annotate("text", x=opt_k, y=max(wss), label=sprintf("k=%d", opt_k),
           vjust=-0.5, color="darkred", fontface="bold", size=3) +
  labs(x="Number of clusters (k)", y="Total within-cluster sum of squares") +
  theme_minimal(base_size=9)

ggsave(file.path(OUT_DIR, "FigS14_microenvironment_elbow_plot.pdf"),
       elbow_plot, width=3.28, height=2.65)

# ===== 6. Cluster spots into microenvironment types =====
cat(sprintf("\nClustering spots into %d microenvironment types...\n", opt_k))

set.seed(42)
km <- kmeans(score_mat, centers=opt_k, nstart=25)
vis$microenvironment <- factor(km$cluster,
                               labels=paste0("ME", 1:opt_k))

cat("Microenvironment cluster sizes:\n")
print(table(vis$microenvironment))

# ===== 7. Characterize each microenvironment by cell-type scores =====
cat("\nCharacterizing microenvironments...\n")

me_char <- data.frame()
for (me in levels(vis$microenvironment)) {
  me_spots <- which(vis$microenvironment == me)
  me_scores <- colMeans(score_mat[me_spots, , drop=FALSE])
  top_cts <- names(sort(me_scores, decreasing=TRUE))

  me_char <- rbind(me_char, data.frame(
    microenvironment = me,
    n_spots = length(me_spots),
    pct_spots = round(length(me_spots) / ncol(vis) * 100, 1),
    top_ct1 = top_cts[1], top_ct1_score = round(me_scores[top_cts[1]], 3),
    top_ct2 = top_cts[2], top_ct2_score = round(me_scores[top_cts[2]], 3),
    top_ct3 = top_cts[3], top_ct3_score = round(me_scores[top_cts[3]], 3),
    stringsAsFactors = FALSE
  ))
}

cat("\nMicroenvironment characterization:\n")
print(me_char)
fwrite(me_char, file.path(OUT_DIR, "TableS_microenvironment_characterization.csv"))

# ===== 8. Tier 1 gene expression per microenvironment =====
cat("\nComputing Tier 1 gene expression per microenvironment...\n")

tier1_found <- intersect(tier1_genes, rownames(vis))
cat(sprintf("Tier 1 genes in Visium: %d / %d\n", length(tier1_found), length(tier1_genes)))

me_expr <- data.frame()
for (g in tier1_found) {
  for (me in levels(vis$microenvironment)) {
    me_spots <- which(vis$microenvironment == me)
    g_expr <- FetchData(vis, cells=colnames(vis)[me_spots], vars=g)
    me_expr <- rbind(me_expr, data.frame(
      gene = g,
      microenvironment = me,
      n_spots = length(me_spots),
      mean_expr = mean(g_expr[[1]]),
      sd_expr = sd(g_expr[[1]]),
      pct_expressed = mean(g_expr[[1]] > 0) * 100,
      stringsAsFactors = FALSE
    ))
  }
}

cat("\nTier 1 gene expression by microenvironment:\n")
print(me_expr[, c("gene", "microenvironment", "mean_expr", "pct_expressed")])

fwrite(me_expr, file.path(OUT_DIR, "TableS_tier1_microenvironment_expression.csv"))

# ===== 9. Visualization: microenvironment map on tissue =====
cat("\nGenerating spatial microenvironment map...\n")

me_spatial <- SpatialDimPlot(vis, group.by="microenvironment", pt.size.factor=1.6,
                              stroke=0, alpha=0.8) +
  scale_fill_brewer(palette="Set2", name="Microenvironment") +
  ggtitle("CRC Tissue Microenvironment Zones",
          subtitle=sprintf("%d microenvironment types (k selected by elbow method)", opt_k)) +
  theme(plot.title=element_text(size=13, face="bold"))

ggsave(file.path(OUT_DIR, "FigX_spatial_microenvironment_map.pdf"),
       me_spatial, width=7, height=6)

# ===== 10. Tier 1 genes x microenvironment heatmap =====
cat("Generating Tier 1 x microenvironment heatmap...\n")

me_expr_norm <- me_expr %>%
  group_by(gene) %>%
  mutate(scaled_expr = scale(mean_expr)) %>%
  ungroup()

hm2 <- ggplot(me_expr_norm, aes(x=microenvironment, y=gene, fill=scaled_expr)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=sprintf("%.2f", mean_expr)), size=3) +
  scale_fill_gradient2(low="steelblue", mid="white", high="darkred",
                       name="Expression\n(Z-score)") +
  labs(title="Tier 1 Genes: Expression by Microenvironment Zone",
       x="Microenvironment Type", y="") +
  theme_minimal(base_size=12) +
  theme(axis.text.y=element_text(face="bold"),
        panel.grid=element_blank())

ggsave(file.path(OUT_DIR, "FigX_spatial_tier1_microenvironment_heatmap.pdf"),
       hm2, width=8, height=5)

# ===== 11. Cell-type score heatmap per microenvironment =====
cat("Generating cell-type score heatmap...\n")

ct_scores_long <- data.frame()
for (me in levels(vis$microenvironment)) {
  me_spots <- which(vis$microenvironment == me)
  me_means <- colMeans(score_mat[me_spots, , drop=FALSE])
  for (ct in names(ct_score_cols)) {
    ct_scores_long <- rbind(ct_scores_long, data.frame(
      microenvironment = me,
      celltype = ct,
      score = me_means[ct],
      stringsAsFactors = FALSE
    ))
  }
}

hm3 <- ggplot(ct_scores_long, aes(x=microenvironment, y=celltype, fill=score)) +
  geom_tile(color="white", linewidth=0.4) +
  scale_fill_gradient2(low="steelblue", mid="white", high="darkred",
                       name="Z-score") +
  labs(title="Cell-Type Signature Scores per Microenvironment Zone",
       x="Microenvironment Type", y="Cell Type") +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(face="bold"),
        axis.text.y=element_text(face="bold"),
        panel.grid=element_blank())

ggsave(file.path(OUT_DIR, "FigX_spatial_celltype_scores_microenvironment.pdf"),
       hm3, width=9, height=6)

# ===== 12. GPBAR1 spatial overlay on microenvironment =====
if ("GPBAR1" %in% tier1_found) {
  cat("\nGenerating GPBAR1 overlay on microenvironment map...\n")

  gpbar1_expr <- FetchData(vis, vars="GPBAR1")
  vis$GPBAR1_expr <- gpbar1_expr$GPBAR1

  gpbar1_spots <- rownames(gpbar1_expr)[gpbar1_expr$GPBAR1 > 0]
  cat(sprintf("GPBAR1 positive spots: %d / %d (%.2f%%)\n",
              length(gpbar1_spots), ncol(vis),
              length(gpbar1_spots)/ncol(vis)*100))

  if (length(gpbar1_spots) > 0) {
    gpbar1_me <- table(vis$microenvironment[gpbar1_spots])
    cat("\nGPBAR1+ spots by microenvironment:\n")
    print(gpbar1_me)

    gpbar1_me_df <- data.frame(
      microenvironment = names(gpbar1_me),
      gpbar1_spots = as.integer(gpbar1_me),
      stringsAsFactors = FALSE
    )
    fwrite(gpbar1_me_df, file.path(OUT_DIR, "TableS_gpbar1_microenvironment.csv"))
  }
}

cat(sprintf("\nAll outputs saved to: %s\n", OUT_DIR))
cat("Done.\n")
