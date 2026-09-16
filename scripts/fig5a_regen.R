#!/usr/bin/env Rscript
# Regenerate ONLY Fig 5a (competitor overlap bar) after the "v5" -> "this study"
# label change.  Isolated from rev3_main_figures.R so no other figure is touched.
suppressMessages({
  library(data.table); library(ggplot2); library(dplyr); library(readr)
  library(tidyr); library(tibble); library(forcats); library(stringr)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "results/figures_rev3"
source("scripts/theme_pub.R")
theme_nc <- theme_pub()

comp <- read_csv("results/competitor_overlap.tsv", show_col_types = FALSE)

v5_sym <- unique(read_tsv("results/phase1_bonferroni_significant_annotated.tsv",
                          show_col_types = FALSE)$SYMBOL)
v5_sym <- v5_sym[!is.na(v5_sym) & v5_sym != ""]
chen_lines_all <- readLines("manuscript/tables/TableS8_chen2024_gene_list.tex")
chen_lines_all <- chen_lines_all[grepl("^\\s*[A-Za-z0-9][A-Za-z0-9._-]*\\s*&", chen_lines_all)]
chen_sym_all <- unique(sub("^\\s*([A-Za-z0-9._-]+)\\s*&.*$", "\\1", chen_lines_all))
chen_sym_all <- chen_sym_all[chen_sym_all != "Gene"]
n_chen_ov <- sum(v5_sym %in% chen_sym_all)
comp <- bind_rows(comp, tibble(
  study = "Chen_2024_NatCommun", journal_IF = "16",
  total_genes = length(chen_sym_all), overlap_n = n_chen_ov,
  overlap_rate = sprintf("%.1f%%", 100 * n_chen_ov / length(v5_sym)),
  overlapping_genes = ""))
cat(sprintf("  Chen 2024: %d/%d overlap\n", n_chen_ov, length(v5_sym)))

comp_plot <- comp %>%
  filter(study != "ALL_COMBINED") %>%
  mutate(study_short = case_when(
    grepl("Hong", study) ~ "Hong 2025\nFront Immunol",
    grepl("Wang_2025", study) ~ "Wang 2025\nDisc Oncol",
    grepl("Wang_suppl", study) ~ "Wang suppl\nBMC Cancer",
    grepl("Tian", study) ~ "Tian 2025\nBiomedicines",
    grepl("Chen", study) ~ "Chen 2024\nNat Commun",
    grepl("Hazelwood", study) ~ "Hazelwood 2025\nNat Commun",
  ),
  study_short = fct_reorder(study_short, overlap_n),
  journal_if_num = as.numeric(journal_IF))

p5c_bar <- ggplot(comp_plot, aes(y = study_short, x = overlap_n)) +
  geom_col(aes(fill = journal_if_num), width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(overlap_n, " (", overlap_rate, ")")),
            hjust = -0.12, size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "#FDDBC7", high = "#2166AC", name = "Journal IF") +
  scale_x_continuous(limits = c(0, max(comp_plot$overlap_n) * 1.35), expand = c(0, 0)) +
  labs(x = "Overlapping genes with this study", y = "",
       title = "Overlap with prior CRC target-discovery studies",
       subtitle = paste0("This study: 69 unique genes. ",
                         comp$overlap_n[comp$study == "ALL_COMBINED"],
                         " overlap (", comp$overlap_rate[comp$study == "ALL_COMBINED"], ")"),
       caption = "Overlap counts are gene-symbol based; journal impact factor shown by fill.") +
  theme_nc +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3))

ggsave(file.path(outdir, "Fig5a_competitor_overlap.pdf"), p5c_bar,
       width = 6.30, height = 3.5, device = cairo_pdf)
cat("  -> Fig5a_competitor_overlap.pdf\n")