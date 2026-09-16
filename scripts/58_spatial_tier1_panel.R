#!/usr/bin/env Rscript
# 58_spatial_tier1_panel.R  (rev7: NATIVE hires background; geometry as rev4)
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(png)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
VISIUM_RDS <- "/ifs1/User/zhouman/CRC_Project/GSE285505_Spatial_Processed.rds"
HIRES_PNG  <- "/ifs1/User/zhouman/CRC-Analysis-Final/GSE285505_Sample1/spatial/tissue_hires_image.png"
OUT <- "results/figures"
genes <- c("UTP11", "BMP2", "PNKD", "LAMC1", "TMBIM1", "GPBAR1")

vis <- readRDS(VISIUM_RDS)
im <- Images(vis)[1]
DefaultAssay(vis) <- if ("SCT" %in% names(vis@assays)) "SCT" else "Spatial"
img <- readPNG(HIRES_PNG)
if (length(dim(img)) == 3 && dim(img)[3] == 4) img <- img[, , 1:3, drop = FALSE]
vis@images[[im]]@image <- img
cat(sprintf("spots=%d assay=%s bg=%s\n", ncol(vis), DefaultAssay(vis), paste(dim(img), collapse=" x ")))

pct <- sapply(genes, function(g) 100 * mean(FetchData(vis, vars = g)[[1]] > 0))
print(round(pct, 2))

plots <- lapply(genes, function(g) {
  lab <- if (pct[[g]] < 1) sprintf("%.2f%% of spots", pct[[g]]) else sprintf("%.1f%% of spots", pct[[g]])
  SpatialFeaturePlot(vis, features = g, pt.size.factor = 1.6, alpha = c(0.3, 1), stroke = 0,
                     image.scale = "hires") +
    labs(title = g, subtitle = lab) +
    theme(plot.title    = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 0)),
          plot.subtitle = element_text(size = 7, hjust = 0.5, margin = margin(t = 0, b = 1)),
          legend.title = element_text(size = 7),
          legend.text  = element_text(size = 7),
          legend.key.width  = unit(0.40, "cm"),
          legend.key.height = unit(0.20, "cm"),
          plot.margin = margin(2, 2, 2, 2))
})
p <- wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave(file.path(OUT, "FigS8_spatial_tier1_genes_raw.pdf"), p, width = 6.30, height = 4.65, device = cairo_pdf)
cat("-> FigS8_spatial_tier1_genes_raw.pdf\n")
cat("=== done ===\n")