#!/usr/bin/env Rscript
# ============================================================
#  rev8_fig1.R -- Figure 1, redrawn to CNS art standard.
#  Panel sizes are identical to the previous round, so the user's
#  Illustrator layout does not change. No number or claim was altered;
#  only the presentation, the hierarchy and the grouping are new.
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2); library(dplyr); library(readr)
  library(tidyr); library(tibble); library(forcats); library(stringr)
  library(ggrepel); library(scales); library(grid)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
source("scripts/theme_pub2.R")
outdir <- "results/figures_rev8_cns"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
th <- theme_pub2()

U <- "\u00d7"; NE <- "\u2260"; MINUS <- "\u2212"; ARR <- "\u2192"; LE <- "\u2264"
SUP6 <- "\u207b\u2076"; EN <- "\u2013"

# ============================================================
#  Figure 1a -- the two evidence layers (conceptual overview)
# ============================================================
cat("\n>>> Figure 1a: two evidence layers\n")

W <- 219; H <- 100                      # drawing units == panel aspect
cols <- list(c1 = c(2, 66), c2 = c(92, 156), card = c(168, 216))
cy   <- c(79, 57.5, 36, 14.5)
bh   <- 6.8

boxes <- data.frame(
  xmin = c(rep(cols$c1[1], 4), rep(cols$c2[1], 4)),
  xmax = c(rep(cols$c1[2], 4), rep(cols$c2[2], 4)),
  id = 1:8, ymin = cy - bh, ymax = cy + bh, ycen = rep(cy, 2),
  layer = rep(c("genetic", "protein"), each = 4),
  fill  = c(rep(COL$tint_blue, 4), rep(COL$tint_amb, 4)),
  accent= c(rep(COL$genetic, 4), rep(COL$protein, 4)),
  line1 = c("15,582 genes tested",
            "75 probes significant",
            "43/76 (57%) colocalised",
            "19 genes Tier 1",
            "62 protein-coding candidates",
            "15 measured in plasma",
            "9 with cis-pQTL",
            "1 colocalised"),
  line2 = c("",
            "p < 3.21e-6",
            "",
            "(SuSiE-converged)",
            "",
            "(17 with TNF pathway)",
            "p < 5e-8",
            "CCM2 (PPH4 = 0.95)"),
  expr = rep(FALSE, 8),
  stringsAsFactors = FALSE)

vseg <- data.frame(
  x    = rep(c(mean(cols$c1), mean(cols$c2)), each = 3),
  y    = rep(cy[1:3] - bh, 2),
  yend = rep(cy[2:4] + bh, 2))

p1a <- ggplot() +
  # vertical chain inside each layer
  geom_segment(data = vseg, aes(x = x, y = y, xend = x, yend = yend),
               colour = COL$ink3, linewidth = 0.45,
               arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  # boxes
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = boxes$fill, colour = NA) +
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmin + 1.7, ymin = ymin, ymax = ymax),
            fill = boxes$accent, colour = NA) +
  geom_text(data = boxes, aes(x = (xmin + xmax) / 2 + 0.9, y = ifelse(line2 == "", ycen, ymax - 3.7), label = line1),
            size = sz(6.9), fontface = "bold", colour = COL$ink) +
  geom_text(data = boxes, aes(x = (xmin + xmax) / 2 + 0.9, y = ymin + 3.7, label = line2),
            size = sz(6.3), colour = COL$ink2) +
  # layer headers
  annotate("text", x = mean(cols$c1), y = 97.5, label = "Genetic layer",
           size = sz(8.2), fontface = "bold", colour = COL$genetic) +
  annotate("text", x = mean(cols$c1), y = 91.5, label = "blood cis-eQTL",
           size = sz(6.4), colour = COL$ink2) +
  annotate("text", x = mean(cols$c2), y = 97.5, label = "Protein layer",
           size = sz(8.2), fontface = "bold", colour = COL$protein) +
  annotate("text", x = mean(cols$c2), y = 91.5, label = "plasma pQTL",
           size = sz(6.4), colour = COL$ink2) +
  # cross-layer arrows
  geom_segment(aes(x = cols$c1[2], y = cy[1], xend = cols$c2[1], yend = cy[1]),
               colour = COL$ink3, linewidth = 0.45,
               arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  annotate("text", x = 79, y = cy[1] + 4.6, label = "protein-coding",
           size = sz(6.1), colour = COL$ink2) +
  geom_segment(aes(x = cols$c1[2], y = cy[3], xend = cols$c2[1], yend = cy[4]),
               colour = COL$alert, linewidth = 0.5, linetype = "dashed",
               arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("text", x = 83, y = 36.5, label = "coverage", size = sz(6.1),
           colour = COL$alert, fontface = "bold") +
  annotate("text", x = 83, y = 31.5, label = "collapse", size = sz(6.1),
           colour = COL$alert, fontface = "bold") +
  geom_segment(aes(x = cols$c2[2], y = cy[3], xend = cols$card[1], yend = cy[3]),
               colour = COL$ink3, linewidth = 0.45,
               arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  # interpretation card
  annotate("rect", xmin = cols$card[1], xmax = cols$card[2], ymin = 32, ymax = 92,
           fill = "white", colour = COL$ink3, linewidth = 0.45) +
  annotate("rect", xmin = cols$card[1], xmax = cols$card[1] + 1.9, ymin = 32, ymax = 92,
           fill = COL$ink, colour = NA) +
  annotate("text", x = mean(cols$card) + 1, y = 86.5, label = "Limiting step",
           size = sz(8.0), fontface = "bold", colour = COL$ink) +
  annotate("text", x = mean(cols$card) + 1, y = 78.5, label = "protein-layer coverage,",
           size = sz(6.6), colour = COL$ink) +
  annotate("text", x = mean(cols$card) + 1, y = 73.5, label = "not genetic significance",
           size = sz(6.6), colour = COL$ink) +
  annotate("segment", x = cols$card[1] + 4, xend = cols$card[2] - 3, y = 66, yend = 66,
           colour = COL$grid, linewidth = 0.8) +
  annotate("text", x = mean(cols$card) + 1, y = 58, label = paste0("eQTL ", NE, " pQTL"),
           size = sz(7.0), fontface = "bold", colour = COL$ink) +
  annotate("text", x = mean(cols$card) + 1, y = 49, label = "5 of 6 core candidates",
           size = sz(6.6), colour = COL$ink) +
  annotate("text", x = mean(cols$card) + 1, y = 44, label = "not measured in",
           size = sz(6.6), colour = COL$ink) +
  annotate("text", x = mean(cols$card) + 1, y = 39, label = "either panel",
           size = sz(6.6), colour = COL$ink) +
  coord_fixed(ratio = 1, xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(2, 3, 2, 2))
save_panel(p1a, file.path(outdir, "Fig1a_conceptual_overview.pdf"),
           mm_snap(154.18), mm_snap(70.36))

# ============================================================
#  Figure 1b -- study design and analytical pipeline
# ============================================================
cat("\n>>> Figure 1b: study design\n")
N_CIS_PQTL <- nrow(fread("results/phase2_decode/phase2_decode_pqtl_mr_combined.csv")) +
              nrow(fread("results/phase2_pqtl/phase2_pqtl_mr_combined.csv"))
stopifnot(N_CIS_PQTL == 9)

BW <- 176.8; BH <- 100
rowsy <- c(82, 52, 22)
band  <- data.frame(ymin = rowsy - 13.5, ymax = rowsy + 13.5,
                    fill = c(COL$tint_blue, COL$tint_blue, COL$tint_grey))
bx    <- list(c(36, 78), c(84, 126), c(132, 174))
boxes_b <- data.frame(
  xmin = rep(c(36, 84, 132), 3), xmax = rep(c(78, 126, 174), 3),
  ycen = rep(rowsy, each = 3),
  bcol = rep(c(COL$genetic, COL$genetic, COL$mute), each = 3),
  l1 = c("SMR: eQTL to GWAS", "Colocalization + SuSiE", "pQTL MR + colocalization",
         "FinnGen R13 replication", "GTEx colon eQTL", "scRNA-seq + Visium spatial",
         "Drug repurposing", "IDG druggability", "PheWAS safety screen"),
  l2 = c("75 genes pass Bonferroni", "43 of 76 probe-gene pairs", "1 of 9 genes with cis-pQTL",
         "1 of 8 Tier 1 genes", "tissue-of-action", "6 microenvironment",
         "8 of 62 genes with", "19 Tier 1 genes", "6 core genes on two"),
  l3 = c("(62 protein-coding)", "colocalize", "instruments is convergent",
         "replicated", "sensitivity", "zones",
         "known drugs", "graded", "transparent dimensions"),
  stringsAsFactors = FALSE)
boxes_b$ymin <- boxes_b$ycen - 11.5; boxes_b$ymax <- boxes_b$ycen + 11.5

labs <- c("1. Discovery", "2. Replication &\ntissue validation", "3. Translation")
band$lab <- labs

chev <- data.frame(x = c(78, 126), y = rep(rowsy, each = 2),
                   yend = rep(rowsy, each = 2), xend = c(84, 132))

p1b <- ggplot() +
  geom_rect(data = band, aes(xmin = 0, xmax = BW, ymin = ymin, ymax = ymax, fill = fill),
            colour = NA) +
  scale_fill_identity() +
  geom_text(data = band, aes(x = 3, y = (ymin + ymax) / 2, label = lab),
            hjust = 0, size = sz(7.6), fontface = "bold", colour = COL$ink,
            lineheight = 1.05) +
  geom_rect(data = boxes_b, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", colour = COL$ink3, linewidth = 0.4) +
  geom_rect(data = boxes_b, aes(xmin = xmin, xmax = xmin + 1.4, ymin = ymin, ymax = ymax,
                                fill = bcol), colour = NA) +
  geom_text(data = boxes_b, aes(x = (xmin + xmax) / 2 + 0.7, y = ymax - 4.6, label = l1),
            size = sz(7.2), fontface = "bold", colour = COL$ink) +
  geom_text(data = boxes_b, aes(x = (xmin + xmax) / 2 + 0.7, y = ymax - 11.0, label = l2),
            size = sz(6.4), colour = COL$ink2) +
  geom_text(data = boxes_b, aes(x = (xmin + xmax) / 2 + 0.7, y = ymax - 16.8, label = l3),
            size = sz(6.4), colour = COL$ink2) +
  geom_segment(data = chev, aes(x = x, y = y, xend = xend, yend = yend),
               colour = COL$ink3, linewidth = 0.5,
               arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  annotate("text", x = 3, y = 1.2, hjust = 0, vjust = 0, size = sz(5.9), colour = COL$ink3,
           lineheight = 1.25,
           label = paste0("Data: 78,473-case CRC GWAS (GCST90255675); eQTLGen (n = 31,684); ",
                          "deCODE and UKB-PPP plasma proteomes;",
                          "\nFinnGen R13; TCGA-COAD and READ; Visium spatial transcriptomics.")) +
  coord_fixed(ratio = 1, xlim = c(0, BW), ylim = c(0, BH), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(3, 4, 2, 3))
save_panel(p1b, file.path(outdir, "Fig1b_study_design.pdf"), mm_snap(169.69), mm_snap(95.96))

# ============================================================
#  Figure 1c -- eQTL-GWAS colocalization screen
# ============================================================
cat("\n>>> Figure 1c: coloc screen\n")
coloc <- read_csv("results/phase3_coloc_all75.csv", show_col_types = FALSE)
if ("SYMBOL" %in% names(coloc)) coloc <- coloc %>% rename(gene = SYMBOL)
if ("pph4"   %in% names(coloc)) coloc <- coloc %>% rename(PPH4 = pph4)

biotype_map <- read_tsv("results/phase1_ensg_biotype.tsv", show_col_types = FALSE) %>%
  select(ensg = ENSG, gene_type, gene_name)

core6 <- c("UTP11", "BMP2", "PNKD", "LAMC1", "TMBIM1", "GPBAR1")

coloc <- coloc %>%
  left_join(biotype_map, by = "ensg") %>%
  mutate(
    gene_label = case_when(
      !grepl("^ENSG", gene) ~ gene,
      !is.na(gene_name) & gene_name != "" & !grepl("^ENSG", gene_name) ~ gene_name,
      !is.na(gene_type) & gene_type == "not_in_GENCODE_v47" ~ paste0(ensg, " (not in GENCODE)"),
      !is.na(gene_type) ~ paste0(ensg, " (", gene_type, ")"),
      TRUE ~ ensg),
    core = gene_label %in% core6)

n_total  <- nrow(coloc)
n_pass   <- sum(coloc$status == "PASS")
n_margin <- sum(coloc$status == "MARGINAL")
n_hidden <- sum(coloc$status %in% c("WEAK", "FAIL"))
n_show   <- n_pass + n_margin

lvl <- c(sprintf("PASS (PPH4 > 0.8; n = %d)", n_pass),
         sprintf("Marginal (0.5%s0.8; n = %d)", EN, n_margin))

cs <- coloc %>%
  filter(status %in% c("PASS", "MARGINAL")) %>%
  mutate(st = factor(status, levels = c("PASS", "MARGINAL"), labels = lvl),
         lab = fct_reorder(gene_label, PPH4))

XLO <- 0.155; XHI <- 1.04; XLAB <- 0.285
set.seed(11)
p1c <- ggplot(cs, aes(y = lab)) +
  geom_segment(aes(x = 0.30, xend = PPH4, yend = lab), colour = "#D8D8D8", linewidth = 0.55) +
  geom_point(aes(x = PPH4, fill = st, shape = st), size = 2.5, stroke = 0.5,
             colour = COL$genetic) +
  geom_text(aes(x = XLAB, label = lab, fontface = ifelse(core, "bold", "plain"),
                colour = ifelse(core, COL$alert, COL$ink2)),
            hjust = 1, size = sz(6.3), show.legend = FALSE) +
  scale_shape_manual(name = NULL, values = c(21, 21), drop = FALSE) +
  scale_fill_manual(name = NULL, values = setNames(c(COL$genetic, COL$genetic2), lvl), drop = FALSE) +
  scale_colour_identity() +
  geom_vline(xintercept = 0.8, linetype = "dashed", colour = COL$alert, linewidth = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = COL$ink3, linewidth = 0.5) +
  annotate("text", x = 0.79, y = n_show + 1.1, label = "PPH4 = 0.8", hjust = 1,
           size = sz(6.0), colour = COL$alert, fontface = "bold") +
  annotate("text", x = 0.495, y = n_show + 1.1, label = "0.5", hjust = 1,
           size = sz(6.0), colour = COL$ink3) +
  scale_x_log10(breaks = c(0.3, 0.5, 0.8, 1.0), labels = c("0.3", "0.5", "0.8", "1.0"),
                limits = c(XLO, XHI), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(1.6, 1.6))) +
  labs(x = "Posterior probability of colocalization (PPH4, log10 scale)", y = NULL,
       caption = sprintf("%d of %d probe-gene pairs shown; %d pairs with PPH4 %s 0.5 omitted.",
                         n_show, n_total, n_hidden, LE)) +
  th + theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key = element_blank(),
    panel.grid.major.y = element_blank(),
    # theme elements take POINTS, not the geom_text size unit: sz() here
    # rendered the caption at 2.14 pt (measured), i.e. unreadable in print.
    plot.caption = element_text(size = 6.1, colour = COL$ink3, hjust = 0,
                                margin = margin(t = 3)),
    plot.margin = margin(6, 8, 4, 4)) +
  guides(fill = guide_legend(override.aes = list(size = 2.6)))
save_panel(p1c, file.path(outdir, "Fig1c_coloc_pph4.pdf"), mm_snap(169.69), mm_snap(188.38))

cat("\n== rev8 Figure 1 done ==\n")