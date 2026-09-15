#!/usr/bin/env Rscript
#'
#' Phase 5c-4: SuSiE Locus Plots for 6 Tier 1 genes
#' Regional association plots showing GWAS signal, eQTL signal, and SuSiE credible sets
#' Uses harmonized GWAS file (chr/pos/rsid) for regional extraction
#'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
})

PROJ <- "/ifs1/User/zhouman/project9-v5-crc-atlas"
setwd(PROJ)
dir.create("results/phase5c_locus_plots", showWarnings = FALSE, recursive = TRUE)

# ── 1. Load Tier 1 genes ──
susie <- fread("results/phase3_susie_finemap.csv")
tier1 <- susie[tier == 1]
cat("Tier 1 genes:", paste(tier1$gene, collapse = ", "), "\n")

# Parse credible set SNPs
parse_cs_snps <- function(cs_text) {
  parts <- strsplit(cs_text, "|", fixed = TRUE)[[1]]
  snps <- gsub(".*(rs\\d+).*", "\\1", parts)
  snps <- snps[grepl("^rs\\d+$", snps)]
  unique(snps)
}

cs_snps_list <- lapply(tier1$cs_summary, parse_cs_snps)
names(cs_snps_list) <- tier1$gene

# ── 2. Load GWAS data (harmonized, has chr/pos/rsid) ──
cat("\nLoading GWAS summary stats (harmonized)...\n")
gwas <- fread("data/gwas/GCST90255675.h.tsv.gz")
setnames(gwas, old = c("chromosome", "base_pair_location", "p_value", "beta", "standard_error", "effect_allele", "other_allele", "rsid"),
         new = c("chr", "pos", "p", "b", "se", "A1", "A2", "SNP"))
gwas <- gwas[grepl("^rs", SNP)]
gwas[, log10p := -log10(p)]
setkey(gwas, chr, pos)
cat(sprintf("  %d GWAS SNPs loaded\n", nrow(gwas)))

# ── 3. Load eQTL SMR results ──
cat("Loading eQTL SMR results...\n")
smr <- fread("results/phase1_eqtl_smr.smr")
# Get gene-level data
smr_genes <- smr[Gene %in% tier1$ensg]
cat(sprintf("  %d SMR gene entries for Tier 1 genes\n", nrow(smr_genes)))

# ── 4. Generate locus plots ──
window_bp <- 500000  # +/-500kb

for (i in seq_len(nrow(tier1))) {
  gene <- tier1$gene[i]
  ensg_i <- tier1$ensg[i]
  chr_i <- tier1$chr[i]
  pos_i <- tier1$pos[i]
  cs_snps <- cs_snps_list[[gene]]
  top_snp_i <- tier1$top_pip_snp[i]
  n_cs_i <- tier1$n_cs[i]
  pph4_i <- tier1$PPH4[i]

  cat(sprintf("\n[%s] chr%s:%s (PPH4=%.3f, %d CSs, %d CS SNPs, top=%s)\n",
              gene, chr_i, pos_i, pph4_i, n_cs_i, length(cs_snps), top_snp_i))

  # GWAS data in region
  chr_n <- as.integer(chr_i)
  gwas_region <- gwas[chr == chr_n & pos >= pos_i - window_bp & pos <= pos_i + window_bp]

  if (nrow(gwas_region) < 10) {
    # Try expanding window
    gwas_region <- gwas[chr == chr_n & pos >= pos_i - 1e6 & pos <= pos_i + 1e6]
  }

  if (nrow(gwas_region) < 10) {
    cat(sprintf("  SKIP: only %d GWAS SNPs in region\n", nrow(gwas_region)))
    next
  }

  cat(sprintf("  %d GWAS SNPs in region\n", nrow(gwas_region)))

  # Annotate credible set SNPs — match by rsID
  gwas_region[, is_cs := SNP %in% cs_snps]
  cs_in_gwas <- gwas_region[is_cs == TRUE]
  cat(sprintf("  %d CS SNPs matched in GWAS data\n", nrow(cs_in_gwas)))

  # Top SNP annotation
  gwas_region[, is_top := SNP == top_snp_i]

  # ── Panel A: Full regional plot ──
  x_range <- range(gwas_region$pos) / 1e6
  x_breaks <- pretty(x_range, n = 5)

  p_full <- ggplot(gwas_region, aes(x = pos / 1e6, y = log10p)) +
    geom_point(aes(color = is_cs, size = is_cs), alpha = 0.6) +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                       labels = c("TRUE" = "In credible set", "FALSE" = "Other"),
                       name = "") +
    scale_size_manual(values = c("TRUE" = 1.5, "FALSE" = 0.5), guide = "none") +
    labs(x = paste0("Chromosome ", chr_i, " position (Mb)"),
         y = expression(-log[10](italic(p))),
         title = paste0(gene, " — GWAS CRC (78K cases)"),
         subtitle = sprintf("PPH4=%.3f, %d CSs, %d CS SNPs, top=%s",
                            pph4_i, n_cs_i, length(cs_snps), top_snp_i)) +
    theme_bw(base_size = 11) +
    theme(legend.position = c(0.85, 0.9),
          legend.background = element_rect(fill = "white", color = "grey80"),
          panel.grid.minor = element_blank())

  # Label CS SNPs
  if (nrow(cs_in_gwas) > 0) {
    p_full <- p_full +
      geom_text_repel(
        data = cs_in_gwas[1:min(10, .N)],
        aes(label = SNP), size = 2.5, color = "#E41A1C",
        max.overlaps = 15, nudge_y = 0.5, segment.size = 0.2)
  }

  # Gene position marker
  p_full <- p_full +
    geom_vline(xintercept = pos_i / 1e6, linetype = "dashed",
               color = "darkgreen", linewidth = 0.8) +
    annotate("text", x = pos_i / 1e6,
             y = min(gwas_region$log10p) + diff(range(gwas_region$log10p)) * 0.95,
             label = gene, color = "darkgreen", fontface = "bold", size = 4, hjust = -0.1)

  # ── Panel B: Zoom on top CS SNP (±100kb) ──
  if (length(top_snp_i) > 0 && top_snp_i %in% gwas_region$SNP) {
    top_pos <- gwas_region[SNP == top_snp_i, pos[1]]
  } else if (nrow(cs_in_gwas) > 0) {
    top_pos <- cs_in_gwas[order(-log10p)][1, pos]
  } else {
    top_pos <- pos_i
  }

  zoom_range <- 100000
  gwas_zoom <- gwas_region[pos >= top_pos - zoom_range & pos <= top_pos + zoom_range]

  if (nrow(gwas_zoom) > 5) {
    p_zoom <- ggplot(gwas_zoom, aes(x = pos / 1e6, y = log10p)) +
      geom_point(aes(color = is_cs, size = is_cs), alpha = 0.7) +
      scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
                         guide = "none") +
      scale_size_manual(values = c("TRUE" = 2, "FALSE" = 1), guide = "none") +
      geom_text_repel(
        data = gwas_zoom[is_cs == TRUE],
        aes(label = SNP), size = 3, color = "#E41A1C", max.overlaps = 20) +
      geom_vline(xintercept = pos_i / 1e6, linetype = "dashed",
                 color = "darkgreen", linewidth = 0.8) +
      labs(x = paste0("Chr", chr_i, " position (Mb)"),
           y = expression(-log[10](italic(p))),
           title = sprintf("Credible set zoom (±100kb around %s)",
                          ifelse(top_snp_i %in% gwas_zoom$SNP, top_snp_i, "top CS SNP"))) +
      theme_bw(base_size = 9) +
      theme(panel.grid.minor = element_blank())
  } else {
    p_zoom <- ggplot() + theme_void() +
      ggtitle("Zoom region: insufficient SNPs")
  }

  # Combine and save
  p_combined <- p_full / p_zoom + plot_layout(heights = c(2, 1))

  ggsave(sprintf("results/phase5c_locus_plots/%s_locus.pdf", gene),
         p_combined, width = 10, height = 8, dpi = 150)
  ggsave(sprintf("results/phase5c_locus_plots/%s_locus.png", gene),
         p_combined, width = 10, height = 8, dpi = 150)

  cat(sprintf("  -> %s_locus.pdf/png\n", gene))
}

cat("\nDone. Plots saved to results/phase5c_locus_plots/\n")
