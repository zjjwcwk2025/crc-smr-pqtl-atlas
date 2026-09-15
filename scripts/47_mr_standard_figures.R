#!/usr/bin/env Rscript
# 47_mr_standard_figures.R
# Fill the MR-standard figure gaps identified in the standards audit:
#   FigS23_tier1_regional_plots.pdf  — 6-panel regional (locus) association plots
#                                      for the six non-MHC Tier 1 genes (GWAS signal
#                                      + SuSiE credible set + gene annotation). This is
#                                      the canonical figure type for a drug-target MR /
#                                      gene-prioritisation study.
#   FigS24_pqtl_mr_scatter.pdf       — SNP-exposure vs SNP-outcome scatter for the
#                                      signal multi-instrument pQTL MR genes (slope =
#                                      IVW causal estimate).
#   FigS25_pqtl_mr_loo.pdf           — leave-one-out IVW plots for those genes
#                                      (robustness to a single instrument).
#
# Raw data (no hard-coded effects):
#   results/phase3_susie_finemap.csv            (Tier 1 genes + SuSiE credible sets)
#   data/gwas/GCST90255675.h.tsv.gz            (harmonised CRC GWAS)
#   results/phase2_decode/{gene}_instruments.csv  (deCODE pQTL + GWAS effects)
#   results/phase2_pqtl/{gene}_instruments.csv    (UKB-PPP pQTL + GWAS effects)
#   results/phase2_decode/phase2_decode_pqtl_mr_combined.csv
#   results/phase2_pqtl/phase2_pqtl_mr_combined.csv

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggrepel)
  library(patchwork)
})

setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "results/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── shared theme (mirrors generate_manuscript_figures.R) ───────────
theme_nc <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text  = element_text(size = 9, color = "black"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ══════════════════════════════════════════════════════════════════
# PART (1) — FigS23: 6-panel regional / locus association plots
# ══════════════════════════════════════════════════════════════════
cat("== FigS23: Tier-1 regional locus plots ==\n")
susie <- fread("results/phase3_susie_finemap.csv")
tier1 <- susie[tier == 1]
cat("  Tier 1 genes:", paste(tier1$gene, collapse = ", "), "\n")

# Parse credible-set SNP rsids from cs_summary ("rs...|rs...|...")
parse_cs <- function(cs_text) {
  if (is.na(cs_text) || cs_text == "") return(character(0))
  parts <- strsplit(cs_text, "|", fixed = TRUE)[[1]]
  snps <- gsub(".*(rs\\d+).*", "\\1", parts)
  snps <- snps[grepl("^rs\\d+$", snps)]
  unique(snps)
}
tier1[, cs_snps := lapply(cs_summary, parse_cs)]

cat("  Loading harmonised GWAS ...\n")
gwas <- fread("data/gwas/GCST90255675.h.tsv.gz")
setnames(gwas, old = c("chromosome", "base_pair_location", "p_value",
                       "beta", "standard_error", "effect_allele",
                       "other_allele", "rsid"),
         new = c("chr", "pos", "p", "b", "se", "A1", "A2", "SNP"))
gwas <- gwas[grepl("^rs", SNP)]
gwas[, log10p := -log10(p)]
setkey(gwas, chr, pos)

WINDOW <- 500e3
p_list <- list()
for (i in seq_len(nrow(tier1))) {
  gene <- tier1$gene[i]; chr_i <- as.integer(tier1$chr[i]); pos_i <- tier1$pos[i]
  cs <- unlist(tier1$cs_snps[[i]])
  top_snp <- tier1$top_pip_snp[i]

  region <- gwas[chr == chr_i & pos >= pos_i - WINDOW & pos <= pos_i + WINDOW]
  if (nrow(region) < 10) region <- gwas[chr == chr_i & pos >= pos_i - 1e6 & pos <= pos_i + 1e6]
  region[, is_cs := SNP %in% cs]

  cs_in_region <- region[is_cs == TRUE]
  lab <- if (nrow(cs_in_region) > 0) head(cs_in_region[order(-log10p)], 6) else region[0]

  p <- ggplot(region, aes(x = pos / 1e6, y = log10p)) +
    geom_point(aes(color = is_cs), alpha = 0.55) +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
    geom_vline(xintercept = pos_i / 1e6, linetype = "dashed",
               color = "darkgreen", linewidth = 0.7) +
    annotate("text", x = pos_i / 1e6,
             y = min(region$log10p) + diff(range(region$log10p)) * 0.93,
             label = gene, color = "darkgreen", fontface = "bold", size = 3.2, hjust = -0.1) +
    labs(x = paste0("Chr", chr_i, " pos (Mb)"), y = expression(-log[10](italic(p))),
         title = sprintf("%s  (PPH4=%.3f)", gene, tier1$PPH4[i])) +
    theme_nc

  if (nrow(lab) > 0) {
    p <- p + geom_text_repel(data = lab, aes(label = SNP), size = 2.2,
                             color = "#E41A1C", max.overlaps = 20,
                             nudge_y = 0.3, segment.size = 0.2)
  }
  p_list[[gene]] <- p
}

p23 <- wrap_plots(p_list, ncol = 3) +
  plot_annotation(title = "Regional association plots — six non-MHC Tier 1 genes",
                  subtitle = "CRC GWAS (-log10 p) with SuSiE credible-set SNPs (red); dashed line = gene position") &
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"))

ggsave(file.path(outdir, "FigS23_tier1_regional_plots.pdf"), p23,
       device = cairo_pdf, width = 12, height = 8)
cat("  -> FigS23_tier1_regional_plots.pdf\n")

# ══════════════════════════════════════════════════════════════════
# PART (2)+(3) — FigS24 scatter + FigS25 leave-one-out for pQTL MR
# ══════════════════════════════════════════════════════════════════
cat("\n== FigS24/FigS25: pQTL MR instrument diagnostics ==\n")

# Significant multi-instrument pQTL MR genes (from combined CSVs) and their source panel.
decode_sig <- c("CCM2", "LIMA1", "STAT6")       # deCODE, Wald p < 0.05
ukb_sig    <- c("CDKN1A", "TNF")                # UKB-PPP, Wald p < 0.05

load_instruments <- function(gene, source) {
  f <- if (source == "deCODE") sprintf("results/phase2_decode/%s_instruments.csv", gene)
       else sprintf("results/phase2_pqtl/%s_instruments.csv", gene)
  d <- fread(f, header = TRUE)
  # harmonise to exposure/outcome effect columns
  d[, .(gene = gene, source = source,
        e_beta = pqtl.beta, e_se = pqtl.se,
        o_beta = gwas.beta, o_se = gwas.se,
        rsid = if ("pqtl.rsid" %in% names(d)) pqtl.rsid else pqtl.snp)]
}

# Gather instrument rows + reported IVW from combined CSV
combined <- data.frame()
instruments <- rbindlist(lapply(c(decode_sig, ukb_sig), function(g) {
  src <- if (g %in% decode_sig) "deCODE" else "UKB-PPP"
  load_instruments(g, src)
}))
instruments <- instruments[e_beta != 0 & !is.na(e_beta) & is.finite(o_beta) & !is.na(o_beta)]

get_ivw <- function(gene) {
  d <- if (gene %in% decode_sig)
    fread("results/phase2_decode/phase2_decode_pqtl_mr_combined.csv")[gene == i]
  else
    fread("results/phase2_pqtl/phase2_pqtl_mr_combined.csv")[gene == i]
  list(b = d$ivw_b, se = d$ivw_se, wald_b = d$b, wald_se = d$se)
}

# Fixed-effect IVW on per-SNP Wald ratios (delta-method SE)
ratio_ivw <- function(dt) {
  r <- dt$o_beta / dt$e_beta
  se_r <- sqrt(dt$o_se^2 / dt$e_beta^2 + dt$o_beta^2 * dt$e_se^2 / dt$e_beta^4)
  w <- 1 / se_r^2
  b <- sum(w * r) / sum(w)
  se <- 1 / sqrt(sum(w))
  c(b = b, se = se)
}

## --- FigS24: SNP-exposure vs SNP-outcome scatter (one facet per gene) ---
scatter_tbl <- instruments[, .(gene, source, e_beta, o_beta)]
scatter_tbl[, gene := factor(gene, levels = c(decode_sig, ukb_sig))]

slopes <- rbindlist(lapply(c(decode_sig, ukb_sig), function(g) {
  iv <- get_ivw(g)
  data.table(gene = g, b = iv$b)
}))
scatter_tbl <- merge(scatter_tbl, slopes, by = "gene", all.x = TRUE)

p24 <- ggplot(scatter_tbl, aes(x = e_beta, y = o_beta)) +
  geom_point(aes(color = gene), alpha = 0.6, size = 1.3) +
  geom_abline(aes(intercept = 0, slope = b, color = gene), linewidth = 0.9, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey60") +
  geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
  facet_wrap(~ gene, scales = "free", ncol = 3) +
  scale_color_manual(values = c("CCM2" = "#2166AC", "LIMA1" = "#B2182B", "STAT6" = "#1B7837",
                                "CDKN1A" = "#9970AB", "TNF" = "#D6604D")) +
  labs(x = "SNP effect on protein level (exposure, pQTL beta)",
       y = "SNP effect on CRC risk (outcome, GWAS beta)",
       title = "pQTL MR — instrument scatter (slope = IVW estimate)") +
  theme_nc + theme(strip.text = element_text(face = "bold"), legend.position = "none")

ggsave(file.path(outdir, "FigS24_pqtl_mr_scatter.pdf"), p24,
       device = cairo_pdf, width = 11, height = 7)
cat("  -> FigS24_pqtl_mr_scatter.pdf\n")

## --- FigS25: leave-one-out IVW (one facet per gene) ---
loo_rows <- list()
for (g in c(decode_sig, ukb_sig)) {
  d <- instruments[gene == g]
  if (nrow(d) < 3) next
  full <- ratio_ivw(d)
  # validate full IVW against reported value
  rep <- get_ivw(g)
  cat(sprintf("  [%s] manual IVW b=%.4f vs reported=%.4f\n", g, full["b"], rep$b))
  mk <- data.table(gene = g, exclude = "ALL", b = full["b"], se = full["se"])
  for (k in seq_len(nrow(d))) {
    sub <- ratio_ivw(d[-k])
    mk <- rbind(mk, data.table(gene = g, exclude = d$rsid[k], b = sub["b"], se = sub["se"]))
  }
  loo_rows[[g]] <- mk
}
loo <- rbindlist(loo_rows)
loo[, lo := b - 1.96 * se]
loo[, hi := b + 1.96 * se]
loo[, exclude := factor(exclude, levels = rev(unique(exclude)))]
loo[, gene := factor(gene, levels = c(decode_sig, ukb_sig))]

p25 <- ggplot(loo, aes(x = b, y = exclude)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = gene), height = 0.25, linewidth = 0.55) +
  geom_point(aes(shape = exclude == "ALL", color = gene), size = 1.6) +
  facet_wrap(~ gene, scales = "free", ncol = 3) +
  scale_color_manual(values = c("CCM2" = "#2166AC", "LIMA1" = "#B2182B", "STAT6" = "#1B7837",
                                "CDKN1A" = "#9970AB", "TNF" = "#D6604D")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16), guide = "none") +
  labs(x = "IVW causal estimate (95% CI), leave-one-out", y = NULL,
       title = "pQTL MR — leave-one-out sensitivity (red 'ALL' = all instruments)") +
  theme_nc + theme(strip.text = element_text(face = "bold"))

ggsave(file.path(outdir, "FigS25_pqtl_mr_loo.pdf"), p25,
       device = cairo_pdf, width = 11, height = 7)
cat("  -> FigS25_pqtl_mr_loo.pdf\n")

# ── postcondition checks ──────────────────────────────────────────
cat("\nPostcondition checks:\n")
f23 <- file.path(outdir, "FigS23_tier1_regional_plots.pdf")
f24 <- file.path(outdir, "FigS24_pqtl_mr_scatter.pdf")
f25 <- file.path(outdir, "FigS25_pqtl_mr_loo.pdf")
stopifnot(file.exists(f23), file.exists(f24), file.exists(f25))
cat("  All 3 PDFs exist.\n")
