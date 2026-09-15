#!/usr/bin/env Rscript
# 思路 5: 6个Tier 1基因之间的空间共定位分析
# Pairwise spatial correlation + cross-correlogram for chr2:219Mb genes
# Output: 6×6 correlation heatmap + cross-correlogram

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
})

set.seed(42)

# ===== Config =====
VISIUM_RDS <- "/ifs1/User/zhouman/CRC_Project/GSE285505_Spatial_Processed.rds"
OUT_DIR    <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase6_scrna"

# Tier 1 genes — read from coloc results (data-driven)
COLOC_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase3_coloc_all75.csv"
GENE_POS_FILE <- "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase1_bonferroni_significant_annotated.tsv"
if (file.exists(COLOC_FILE)) {
  coloc <- fread(COLOC_FILE)
  # Coloc file uses 'gene' column
  tier1_raw <- coloc[PPH4 > 0.8 & status == "PASS", unique(gene)]

  coding_regex <- "^(RP|AC|AL|AP|LINC|MIR|SNOR|[A-Z]+[0-9]+\\.[0-9]+)"
  TIER1_GENES <- tier1_raw[!grepl(coding_regex, tier1_raw)]
} else {
  TIER1_GENES <- c("UTP11", "BMP2", "PNKD", "LAMC1", "TMBIM1", "GPBAR1")
}

# Chr2 genes: derived from Tier 1 on chr2 (data-driven)
# Gene position file uses SYMBOL + ProbeChr columns
if (file.exists(GENE_POS_FILE)) {
  gene_pos <- fread(GENE_POS_FILE)
  chr2_tier1 <- gene_pos[SYMBOL %in% TIER1_GENES & ProbeChr == 2, unique(SYMBOL)]
  CHR2_GENES <- intersect(chr2_tier1, TIER1_GENES)
} else {
  CHR2_GENES <- c("PNKD", "TMBIM1", "GPBAR1")  # fallback
}
cat(sprintf("Tier 1 genes: %s\n", paste(TIER1_GENES, collapse=", ")))
cat(sprintf("Chr2 genes: %s\n", paste(CHR2_GENES, collapse=", ")))

# Cross-correlogram: data-driven distance parameters
# Visium spot diameter = 55µm, center-to-center spacing ≈ 100µm
# max_dist = 30 × spacing covers spatial decay up to ~3mm tissue span
VISIUM_SPOT_SPACING_UM <- 100
CROSSCORR_MAX_DIST_UM  <- VISIUM_SPOT_SPACING_UM * 30  # 3000µm
CROSSCORR_N_BINS       <- 30
CROSSCORR_MIN_PAIRS     <- 10L

# ===== 1. Load Visium =====
cat("Loading Visium data...\n")
vis <- readRDS(VISIUM_RDS)
cat(sprintf("Visium: %d spots x %d genes\n", ncol(vis), nrow(vis)))

# Check gene availability
found <- intersect(TIER1_GENES, rownames(vis))
missing <- setdiff(TIER1_GENES, rownames(vis))
cat(sprintf("Tier 1 genes found in Visium: %d / %d\n", length(found), length(TIER1_GENES)))
if (length(missing) > 0) cat(sprintf("Missing: %s\n", paste(missing, collapse=", ")))

# ===== 2. Extract expression matrix for Tier 1 genes =====
expr_mat <- FetchData(vis, vars = found)
cat(sprintf("Expression matrix: %d spots x %d genes\n", nrow(expr_mat), ncol(expr_mat)))

# ===== 3. Pairwise Spearman correlation (all 6 Tier 1 genes) =====
cat("\nComputing pairwise Spearman correlations...\n")

cor_matrix <- matrix(NA, nrow=length(found), ncol=length(found))
rownames(cor_matrix) <- found
colnames(cor_matrix) <- found

p_matrix <- cor_matrix

for (i in seq_along(found)) {
  for (j in seq_along(found)) {
    if (i == j) {
      cor_matrix[i,j] <- 1.0
      p_matrix[i,j] <- 0
    } else if (i < j) {
      ct <- cor.test(expr_mat[[found[i]]], expr_mat[[found[j]]], method="spearman")
      cor_matrix[i,j] <- ct$estimate
      cor_matrix[j,i] <- ct$estimate
      p_matrix[i,j] <- ct$p.value
      p_matrix[j,i] <- ct$p.value
    }
  }
}

cat("Pairwise Spearman correlation matrix:\n")
print(round(cor_matrix, 4))

# Save correlation matrix
cor_df <- melt(cor_matrix, varnames=c("gene1", "gene2"), value.name="spearman_rho")
p_df   <- melt(p_matrix, varnames=c("gene1", "gene2"), value.name="p_value")
cor_full <- merge(cor_df, p_df, by=c("gene1", "gene2"))
cor_full <- cor_full[cor_full$gene1 != cor_full$gene2, ]
fwrite(cor_full, file.path(OUT_DIR, "TableS_spatial_gene_gene_correlation.csv"))

# ===== 4. Heatmap: 6×6 pairwise correlation =====
cat("\nGenerating correlation heatmap...\n")

cor_df_plot <- cor_df
cor_df_plot$label <- sprintf("%.3f", round(cor_df_plot$spearman_rho, 3))
cor_df_plot$gene1 <- factor(cor_df_plot$gene1, levels=found)
cor_df_plot$gene2 <- factor(cor_df_plot$gene2, levels=rev(found))

hm <- ggplot(cor_df_plot, aes(x=gene1, y=gene2, fill=spearman_rho)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=label), size=3.5, color="grey20") +
  scale_fill_gradient2(low="steelblue", mid="white", high="darkred",
                       midpoint=0, limits=c(-1,1), name="Spearman rho") +
  labs(title="Tier 1 Genes: Pairwise Spatial Expression Correlation (Visium)",
       subtitle=sprintf("n = %d spots, SCT-normalized expression", nrow(expr_mat)),
       x="", y="") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1, face="bold"),
        axis.text.y=element_text(face="bold"),
        panel.grid=element_blank())

ggsave(file.path(OUT_DIR, "FigX_spatial_tier1_correlation_heatmap.pdf"),
       hm, width=7, height=6)

# ===== 5. Spatial cross-correlogram for chr2:219Mb three genes =====
cat("\nComputing spatial cross-correlograms for chr2:219Mb genes...\n")

# Get spatial coordinates (Seurat v5 / VisiumV2 compatible)
coords <- GetTissueCoordinates(vis)
coords <- coords[, c("x", "y")]

chr2_expr <- as.data.frame(expr_mat[, CHR2_GENES, drop=FALSE])
chr2_expr$spot_id <- rownames(chr2_expr)

# Compute pairwise Euclidean distances between all spot pairs
coords_mat <- as.matrix(coords)
spot_dist <- as.matrix(dist(coords_mat))

# Helper: cross-correlogram for a gene pair
compute_crosscorrelogram <- function(g1, g2, dist_mat, expr_data,
                                      max_dist=CROSSCORR_MAX_DIST_UM,
                                      n_bins=CROSSCORR_N_BINS,
                                      min_pairs=CROSSCORR_MIN_PAIRS) {
  g1_expr <- expr_data[[g1]]
  g2_expr <- expr_data[[g2]]

  bin_edges <- seq(0, max_dist, length.out=n_bins+1)
  bin_centers <- (bin_edges[-1] + bin_edges[-(n_bins+1)]) / 2

  cors <- numeric(n_bins)
  ns   <- integer(n_bins)

  for (b in seq_len(n_bins)) {
    in_bin <- dist_mat > bin_edges[b] & dist_mat <= bin_edges[b+1]
    # Upper triangle only to avoid double-counting
    in_bin[lower.tri(in_bin, diag=TRUE)] <- FALSE

    idx <- which(in_bin, arr.ind=TRUE)
    ns[b] <- nrow(idx)

    if (ns[b] >= min_pairs) {
      cors[b] <- cor(g1_expr[idx[,1]], g2_expr[idx[,2]], method="spearman")
    } else {
      cors[b] <- NA
    }
  }

  data.frame(
    pair = paste(g1, g2, sep="_"),
    gene1 = g1,
    gene2 = g2,
    distance_um = bin_centers,
    spearman_rho = cors,
    n_pairs = ns,
    stringsAsFactors = FALSE
  )
}

# Compute for all chr2 gene pairs
chr2_pairs <- combn(CHR2_GENES, 2, simplify=FALSE)
crosscorr_list <- list()

for (pair in chr2_pairs) {
  cat(sprintf("  Cross-correlogram: %s vs %s...\n", pair[1], pair[2]))
  crosscorr_list[[paste(pair, collapse="_")]] <- compute_crosscorrelogram(
    pair[1], pair[2], spot_dist, chr2_expr
  )
}

crosscorr_df <- rbindlist(crosscorr_list)
fwrite(crosscorr_df, file.path(OUT_DIR, "TableS_spatial_crosscorrelogram_chr2.csv"))

# ===== 6. Plot cross-correlogram =====
cat("Plotting cross-correlogram...\n")

cc_plot <- ggplot(crosscorr_df, aes(x=distance_um/1000, y=spearman_rho, color=pair, group=pair)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line(linewidth=1) +
  geom_point(size=1.5) +
  scale_color_manual(values=c("PNKD_TMBIM1"="darkred",
                               "PNKD_GPBAR1"="steelblue",
                               "TMBIM1_GPBAR1"="darkgreen"),
                      labels=c("PNKD_TMBIM1"="PNKD vs TMBIM1",
                               "PNKD_GPBAR1"="PNKD vs GPBAR1",
                               "TMBIM1_GPBAR1"="TMBIM1 vs GPBAR1")) +
  labs(title="Cross-correlogram: chr2:219Mb Tier 1 Genes",
       subtitle="Spatial correlation decays with distance",
       x="Distance between spots (mm)", y="Spearman rho",
       color="Gene Pair") +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom")

ggsave(file.path(OUT_DIR, "FigS17_spatial_chr2_crosscorrelogram.pdf"),
       cc_plot, width=8, height=5)

# ===== 7. Summary statistics =====
cat("\n===== Summary =====\n")

# Mean pairwise correlation
pairwise <- cor_full[!duplicated(paste(pmin(as.character(cor_full$gene1), as.character(cor_full$gene2)),
                                       pmax(as.character(cor_full$gene1), as.character(cor_full$gene2)))), ]
cat(sprintf("Total unique gene pairs: %d\n", nrow(pairwise)))
cat(sprintf("Significant pairs (p<0.05): %d\n", sum(pairwise$p_value < 0.05, na.rm=TRUE)))

# chr2 cluster analysis
cat("\n--- chr2:219Mb intra-cluster correlations ---\n")
for (g1 in CHR2_GENES) {
  for (g2 in CHR2_GENES) {
    if (g1 < g2) {
      r <- cor_matrix[g1, g2]
      p <- p_matrix[g1, g2]
      cat(sprintf("  %s vs %s: ρ=%.4f, p=%.2e\n", g1, g2, r, p))
    }
  }
}

# Cross-cluster: chr2 genes vs non-chr2 Tier 1
non_chr2 <- setdiff(found, CHR2_GENES)
cat("\n--- chr2 vs non-chr2 correlations ---\n")
for (g1 in CHR2_GENES) {
  for (g2 in non_chr2) {
    r <- cor_matrix[g1, g2]
    p <- p_matrix[g1, g2]
    cat(sprintf("  %s vs %s: ρ=%.4f, p=%.2e\n", g1, g2, r, p))
  }
}

cat(sprintf("\nAll outputs saved to: %s\n", OUT_DIR))
cat("Done.\n")
