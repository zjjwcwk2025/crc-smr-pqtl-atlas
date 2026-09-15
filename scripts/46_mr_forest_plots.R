#!/usr/bin/env Rscript
# 46_mr_forest_plots.R
# Generate MR-related publication figures (previously missing).
# Three forest plots, one per MR group, following generate_manuscript_figures.R conventions
# (theme_nc, palette_direction, cairo_pdf).
#
#   (1) FigS20_pqtl_mr_decode_forest.pdf   — deCODE  pQTL MR (5 genes), Wald + IVW, 95% CI
#   (2) FigS21_pqtl_mr_ukbppp_forest.pdf   — UKB-PPP pQTL MR (4 genes), Wald + IVW, 95% CI
#   (3) FigS22_finngen_replication_forest.pdf — FinnGen eQTL-tool Mr (8 genes), IVW, concordance colour
#
# Reads from the canonical result CSVs (no hard-coded values):
#   results/phase2_decode/phase2_decode_pqtl_mr_combined.csv
#   results/phase2_pqtl/phase2_pqtl_mr_combined.csv
#   results/phase4_finngen_replication.csv
#   results/fstatistics_summary.csv   (per-gene mean F for annotation)

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "results/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── theme (mirror of generate_manuscript_figures.R) ───────────
theme_nc <- theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 11, color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(12, 12, 12, 12)
  )

palette_direction <- c("concordant" = "#2166AC", "discordant" = "#B2182B",
                       "not_significant" = "#999999")

# ── load data ─────────────────────────────────────────────────
decode <- read_csv("results/phase2_decode/phase2_decode_pqtl_mr_combined.csv",
                   show_col_types = FALSE)
ukb   <- read_csv("results/phase2_pqtl/phase2_pqtl_mr_combined.csv",
                  show_col_types = FALSE)
fg    <- read_csv("results/phase4_finngen_replication.csv",
                  show_col_types = FALSE)
fstat <- read_csv("results/fstatistics_summary.csv", show_col_types = FALSE)

# Build a clean annotation table: gene -> (source, n_instruments, f_mean)
fstat <- fstat %>%
  filter(!is.na(n_instruments)) %>%
  select(gene, source, n_instruments, f_mean)

# ── (1) deCODE pQTL MR forest ─────────────────────────────────
# Wald ratio (single top instrument) and IVW (all instruments), 95% CI = b ± 1.96*se
decode_plot <- decode %>%
  select(gene, wald_b, wald_se, wald_pval, wald_significant, ivw_b, ivw_se, ivw_pval) %>%
  mutate(
    wald_ci_low  = wald_b - 1.96 * wald_se,
    wald_ci_high = wald_b + 1.96 * wald_se,
    ivw_ci_low   = ivw_b - 1.96 * ivw_se,
    ivw_ci_high  = ivw_b + 1.96 * ivw_se
  ) %>%
  select(gene, wald_b, wald_ci_low, wald_ci_high, wald_significant,
         ivw_b, ivw_ci_low, ivw_ci_high, ivw_pval) %>%
  mutate(gene = factor(gene, levels = rev(unique(gene))))

subtitle_decode <- sprintf(
  "Wald ratio (single top cis-pQTL) and IVW (all instruments).\nn = %d genes; shaded = p < 0.05.",
  nrow(decode_plot))

p_decode <- ggplot(decode_plot, aes(y = gene)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_errorbarh(aes(xmin = ivw_ci_low, xmax = ivw_ci_high, color = "IVW"),
                 height = 0.22, linewidth = 0.9,
                 position = position_nudge(y = -0.17)) +
  geom_point(aes(x = ivw_b, color = "IVW"), size = 2.4, shape = 17,
             position = position_nudge(y = -0.17)) +
  geom_errorbarh(aes(xmin = wald_ci_low, xmax = wald_ci_high, color = "Wald"),
                 height = 0.22, linewidth = 0.9,
                 position = position_nudge(y = 0.17)) +
  geom_point(aes(x = wald_b, color = "Wald"), size = 2.4, shape = 16,
             position = position_nudge(y = 0.17)) +
  scale_color_manual(values = c("Wald" = "#2166AC", "IVW" = "#B2182B"),
                     name = NULL) +
  coord_cartesian(xlim = range(c(decode_plot$wald_ci_low, decode_plot$wald_ci_high,
                                 decode_plot$ivw_ci_low, decode_plot$ivw_ci_high)) * 1.15) +
  labs(x = "Effect on CRC risk per SD (b, 95% CI)", y = NULL,
       title = "pQTL MR — deCODE instruments",
       subtitle = subtitle_decode) +
  theme_nc

# ── (2) UKB-PPP pQTL MR forest ────────────────────────────────
ukb_plot <- ukb %>%
  select(gene, n_instruments, b, se, pval, ivw_b, ivw_se, ivw_pval) %>%
  mutate(
    wald_b      = b,
    wald_pval   = pval,
    wald_ci_low  = b - 1.96 * se,
    wald_ci_high = b + 1.96 * se,
    ivw_ci_low   = ivw_b - 1.96 * ivw_se,
    ivw_ci_high  = ivw_b + 1.96 * ivw_se,
    gene = factor(gene, levels = rev(unique(gene)))
  )

subtitle_ukb <- sprintf(
  "Wald ratio (single top cis-pQTL) and IVW (all instruments). n = %d genes. UKB-PPP plasma pQTL.",
  nrow(ukb_plot))

p_ukb <- ggplot(ukb_plot, aes(y = gene)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_errorbarh(aes(xmin = ivw_ci_low, xmax = ivw_ci_high, color = "IVW"),
                 height = 0.22, linewidth = 0.9, position = position_nudge(y = -0.17)) +
  geom_point(aes(x = ivw_b, color = "IVW"), size = 2.4, shape = 17,
             position = position_nudge(y = -0.17)) +
  geom_errorbarh(aes(xmin = wald_ci_low, xmax = wald_ci_high, color = "Wald"),
                 height = 0.22, linewidth = 0.9, position = position_nudge(y = 0.17)) +
  geom_point(aes(x = wald_b, color = "Wald"), size = 2.4, shape = 16,
             position = position_nudge(y = 0.17)) +
  scale_color_manual(values = c("Wald" = "#2166AC", "IVW" = "#B2182B"), name = NULL) +
  coord_cartesian(xlim = range(c(ukb_plot$wald_ci_low, ukb_plot$wald_ci_high,
                                 ukb_plot$ivw_ci_low, ukb_plot$ivw_ci_high)) * 1.15) +
  labs(x = "Effect on CRC risk per SD (b, 95% CI)", y = NULL,
       title = "pQTL MR — UKB-PPP instruments",
       subtitle = subtitle_ukb) +
  theme_nc

# ── (3) FinnGen replication forest ────────────────────────────
fg_plot <- fg %>%
  select(gene, n_instruments, ivw_b, ivw_se, ivw_p, dir_match, replication) %>%
  mutate(
    ci_low  = ivw_b - 1.96 * ivw_se,
    ci_high = ivw_b + 1.96 * ivw_se,
    concord = case_when(
      str_starts(dir_match, "✅") ~ "concordant",
      str_starts(dir_match, "⚠") ~ "discordant",
      TRUE ~ "not_significant"
    ),
    rep_lab = ifelse(str_detect(replication, "Full"), "Full", "Power-limited"),
    gene = factor(gene, levels = rev(unique(gene)))
  )

subtitle_fg <- sprintf(
  "eQTL-IVW MR of CRC risk in FinnGen R13. n = %d genes; colour = direction vs discovery SMR.",
  nrow(fg_plot))

p_fg <- ggplot(fg_plot, aes(x = ivw_b, y = gene, color = concord)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.22, linewidth = 0.9) +
  geom_point(size = 2.6, shape = 16) +
  scale_color_manual(values = palette_direction,
                     limits = c("concordant", "discordant", "not_significant"),
                     name = "Direction vs discovery", labels = c("Concordant", "Discordant", "n.s.")) +
  labs(x = "IVW effect on CRC risk (b, 95% CI)", y = NULL,
       title = "FinnGen replication — eQTL instrument MR",
       subtitle = subtitle_fg) +
  theme_nc

# ── save all ──────────────────────────────────────────────────
ggsave(file.path(outdir, "FigS20_pqtl_mr_decode_forest.pdf"), p_decode,
       device = cairo_pdf, width = 8, height = 5)
ggsave(file.path(outdir, "FigS21_pqtl_mr_ukbppp_forest.pdf"), p_ukb,
       device = cairo_pdf, width = 8, height = 4.5)
ggsave(file.path(outdir, "FigS22_finngen_replication_forest.pdf"), p_fg,
       device = cairo_pdf, width = 8, height = 6)

cat("Saved 3 MR forest plots to", outdir, "\n")
cat("  - FigS20_pqtl_mr_decode_forest.pdf\n")
cat("  - FigS21_pqtl_mr_ukbppp_forest.pdf\n")
cat("  - FigS22_finngen_replication_forest.pdf\n")

# postcondition checks (iron rule 5)
stopifnot(file.exists(file.path(outdir, "FigS20_pqtl_mr_decode_forest.pdf")))
stopifnot(file.exists(file.path(outdir, "FigS21_pqtl_mr_ukbppp_forest.pdf")))
stopifnot(file.exists(file.path(outdir, "FigS22_finngen_replication_forest.pdf")))
cat("Postcondition: all 3 PDFs exist.\n")
