#!/usr/bin/env Rscript
# =====================================================================
#  Figure 1c (eQTL-GWAS colocalization screen, 49/76 pairs) in a compact
#  transposed layout: 49 genes on the x axis (labels rotated 90 deg),
#  log10 PPH4 on the y axis.  Same data, colours, thresholds and wording
#  as the full-height panel (scripts/rev4_printsizes.R, block "Figure 1c"),
#  but ~74 mm tall instead of 178 mm, so that Figure 1 built from panels
#  (a)+(b)+(c) stacked at native size stays inside the 225 mm height limit
#  without shrinking any text below 7 pt.
#  Overrides results/figures_rev3/Fig1c_coloc_pph4.pdf (tall version kept
#  in ~/work_bi/fig1c_orig/).
# =====================================================================
suppressMessages({
  library(data.table); library(ggplot2); library(dplyr); library(readr)
  library(stringr)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "results/figures_rev3"
altdir <- "results/figures_rev6"
bak    <- "/ifs1/User/zhouman/work_bi/fig1c_orig"
dir.create(bak, showWarnings = FALSE, recursive = TRUE)
for (d in c(outdir, altdir)) {
  f <- file.path(d, "Fig1c_coloc_pph4.pdf")
  if (file.exists(f) && !file.exists(file.path(bak, basename(d), "Fig1c_coloc_pph4.pdf"))) {
    dir.create(file.path(bak, basename(d)), showWarnings = FALSE, recursive = TRUE)
    file.copy(f, file.path(bak, basename(d), "Fig1c_coloc_pph4.pdf"))
  }
}
source("scripts/theme_pub.R")
theme_nc <- theme_pub()

coloc <- read_csv("results/phase3_coloc_all75.csv", show_col_types = FALSE)
if ("SYMBOL" %in% names(coloc)) coloc <- coloc %>% rename(gene = SYMBOL)
if ("pph4" %in% names(coloc))   coloc <- coloc %>% rename(PPH4 = pph4)
biotype_map <- read_tsv("results/phase1_ensg_biotype.tsv", show_col_types = FALSE) %>%
  select(ensg = ENSG, gene_type, gene_name)
coloc <- coloc %>%
  left_join(biotype_map, by = "ensg") %>%
  mutate(gene_label = case_when(
    !grepl("^ENSG", gene) ~ gene,
    !is.na(gene_name) & gene_name != "" & !grepl("^ENSG", gene_name) ~ gene_name,
    !is.na(gene_type) & gene_type == "not_in_GENCODE_v47" ~ paste0(ensg, " (not in GENCODE)"),
    !is.na(gene_type) ~ paste0(ensg, " (", gene_type, ")"),
    TRUE ~ ensg))

n_total  <- nrow(coloc)
n_pass   <- sum(coloc$status == "PASS")
n_margin <- sum(coloc$status == "MARGINAL")
n_hidden <- sum(coloc$status %in% c("WEAK", "FAIL"))
n_show   <- n_pass + n_margin

status_lvl <- c(sprintf("PASS (PPH4 > 0.8; n = %d)", n_pass),
                sprintf("Marginal (0.5\u20130.8; n = %d)", n_margin))
status_vals <- setNames(c("#2166AC", "#92C5DE"), status_lvl)

cs <- coloc %>%
  filter(status %in% c("PASS", "MARGINAL")) %>%
  mutate(coloc_status = factor(status, levels = c("PASS", "MARGINAL"), labels = status_lvl),
         PPH4_plot = PPH4) %>%
  arrange(desc(PPH4_plot), gene_label) %>%
  mutate(gene_label = factor(gene_label, levels = rev(gene_label)))

cat(sprintf("pairs shown: %d of %d  (PASS %d + marginal %d; omitted %d)\n",
            nrow(cs), n_total, n_pass, n_margin, n_hidden))
cat(sprintf("longest y label: %d chars\n", max(nchar(as.character(cs$gene_label)))))

p <- ggplot(cs, aes(x = gene_label, y = PPH4_plot)) +
  geom_segment(aes(x = gene_label, xend = gene_label, y = 0.30, yend = PPH4_plot),
               colour = "#BDBDBD", linewidth = 0.45) +
  geom_point(aes(colour = coloc_status), size = 1.9) +
  scale_colour_manual(name = NULL, values = status_vals, drop = FALSE) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "#B2182B", linewidth = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "#6B6B6B", linewidth = 0.5) +
  scale_y_log10(breaks = c(0.3, 0.5, 0.8, 1.0), labels = c("0.3", "0.5", "0.8", "1.0"),
                limits = c(0.28, 1.03), expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = NULL, y = "Posterior Probability of Colocalization (PPH4, log10)") +
  theme_nc +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5,
                                   margin = margin(t = 2)),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(10, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        plot.margin = margin(t = 3, r = 4, b = 1, l = 4))

W <- 6.69; H <- 2.90
for (d in c(outdir, altdir)) {
  ggsave(file.path(d, "Fig1c_coloc_pph4.pdf"), p, width = W, height = H, device = cairo_pdf)
}
preview <- "/ifs1/User/zhouman/work_bi/fig1c_compact"
dir.create(preview, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(preview, "Fig1c_coloc_pph4.png"), p, width = W, height = H, dpi = 300)
cat(sprintf("saved compact Fig1c  %.2f x %.2f in  (%.1f x %.1f mm)\n", W, H, W*25.4, H*25.4))
cat("== 75_fig1c_compact done ==\n")