#!/usr/bin/env Rscript
# ============================================================
#  rev8_fig2.R -- Figure 2, redrawn to CNS art standard.
#  Panel sizes identical to the previous round.
#  Logic change (no number or claim altered):
#   the old Fig 2b put "pQTL tested and null" and "no cis-pQTL
#   instrument at all" in one class C. The panel now separates them,
#   which is the paper's own point: coverage, not significance, is
#   the bottleneck. The split is derived from the data, not typed in.
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2); library(dplyr); library(readr)
  library(tidyr); library(tibble); library(forcats); library(stringr)
  library(ggrepel); library(scales)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
source("scripts/theme_pub2.R")
outdir <- "results/figures_rev8_cns"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
th <- theme_pub2()

U <- "\u00d7"; LE <- "\u2264"; SUP6 <- "\u207b\u2076"; SUP8 <- "\u207b\u2078"; EN <- "\u2013"

# ============================================================
#  Figure 2a -- protein-layer coverage funnel
# ============================================================
cat("\n>>> Figure 2a: coverage funnel\n")
bio  <- fread("results/phase1_ensg_biotype.tsv")
n_coding <- sum(bio$gene_type == "protein_coding")
dec  <- unique(sub("^([0-9]+_[0-9]+_)([^_]+)_.*$", "\\2", list.files("data/decode", pattern = "\\.txt\\.gz$")))
ukb  <- unique(sub("_.*$", "", list.files("data/ukbppp_merged", pattern = "_merged\\.txt\\.gz$")))
coding <- bio$gene_name[bio$gene_type == "protein_coding"]
all_meas <- union(dec, ukb)
measured <- all_meas[all_meas %in% coding]
extra    <- setdiff(all_meas, coding)
with_inst <- c("CDKN1A","LTA","NCF2","TNF","CCM2","LIMA1","STAT6","TNFRSF1A","TNFRSF1B")
n_inst_meas <- sum(with_inst %in% measured)
stopifnot(n_coding == 62, length(measured) == 15, length(extra) == 2,
          length(all_meas) == 17, sum(with_inst %in% all_meas) == 9, n_inst_meas == 7)

wrap_lab <- function(x, w) paste(strwrap(x, width = w), collapse = "\n")
fun <- data.frame(
  y  = c(110, 80, 50, 20),
  n  = c(n_coding, length(measured), sum(with_inst %in% measured), 1),
  head = c("Protein-coding genes",
           "With plasma protein measurement",
           "With cis-pQTL instruments",
           "With pQTL-GWAS colocalization"),
  detail = c("of 75 genes passing Bonferroni-corrected SMR (p < 3.21e-6)",
             sprintf("deCODE %d + UKB-PPP %d (%d in both); plus 2 TNF-pathway genes outside the Bonferroni set = 17-gene analysis set",
                     length(intersect(dec, coding)), length(intersect(ukb, coding)), length(intersect(dec, ukb))),
             sprintf("%d of these 15; %d of the 17-gene analysis set",
                     n_inst_meas, sum(with_inst %in% all_meas)),
             "CCM2 only (PPH4 = 0.95)"),
  fill = c(COL$genetic, "#5C8FC0", COL$protein, COL$alert),
  stringsAsFactors = FALSE)
stopifnot(fun$n == c(62, 15, 7, 1))
fun$detail <- vapply(fun$detail, wrap_lab, character(1), w = 120)
fun$bar <- 84 * fun$n / 62
# the "five of six have no protein measurement" sentence is already the first
# line of the published panel-a legend, so it is not repeated inside the panel

p2a <- ggplot(fun) +
  geom_rect(aes(xmin = 0, xmax = bar, ymin = y - 4.6, ymax = y + 4.6, fill = fill), colour = NA) +
  scale_fill_identity() +
  geom_text(aes(x = bar + 2.2, y = y, label = n), hjust = 0, size = sz(9.6),
            fontface = "bold", colour = COL$ink) +
  geom_text(aes(x = 0, y = y + 10.4, label = head), hjust = 0, size = sz(6.9),
            fontface = "bold", colour = COL$ink) +
  geom_text(aes(x = 0, y = y - 6.4, label = detail), hjust = 0, size = sz(6.0),
            colour = COL$ink3, lineheight = 0.95, vjust = 1) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 126), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(3, 3, 2, 3))
save_panel(p2a, file.path(outdir, "Fig2a_coverage_funnel.pdf"), mm_snap(169.69), mm_snap(60))

# ============================================================
#  Figure 2b -- can the protein layer be tested at all?
# ============================================================
cat("\n>>> Figure 2b: protein-layer testability\n")
con <- read_csv("results/phase2_decode/eqtl_pqtl_contrast_matrix.csv", show_col_types = FALSE) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "")

con <- con %>%
  mutate(
    tested = !is.na(pqtl_wald_p),
    cls = case_when(
      !tested                                    ~ "D",
      pqtl_wald_p < 0.05 & direction_match == TRUE  ~ "A",
      pqtl_wald_p < 0.05 & direction_match == FALSE ~ "B",
      TRUE                                       ~ "C"),
    grp = case_when(cls == "A" ~ "A  pQTL concordant and significant",
                    cls == "B" ~ "B  pQTL discordant and significant",
                    cls == "C" ~ "C  cis-pQTL tested, not significant",
                    TRUE       ~ "D  no cis-pQTL instrument available"),
    fill = c(A = COL$genetic, B = COL$alert, C = "#C9C9C9", D = COL$protein)[cls])

tab <- con %>% count(cls, grp, fill) %>% arrange(cls)
stopifnot(sum(tab$n[tab$cls == "A"]) == 4, sum(tab$n[tab$cls == "B"]) == 1,
          sum(tab$n[tab$cls == "C"]) == 4, sum(tab$n[tab$cls == "D"]) == 8,
          sum(tab$n) == 17)
genes_of <- function(k) {
  g <- sort(con$SYMBOL[con$cls == k]); paste(strwrap(paste(g, collapse = ", "), width = 58),
                                             collapse = "\n") }
tab$genes <- vapply(tab$cls, genes_of, character(1))
tab$y  <- c(98, 74, 50, 20)
tab$bar <- 84 * tab$n / 17
tab$lbl <- sprintf("%d  (%.0f%%)", tab$n, 100 * tab$n / 17)
tab$short <- c("A  pQTL concordant and significant",
               "B  pQTL discordant and significant",
               "C  cis-pQTL tested, not significant",
               "D  no cis-pQTL instrument available")

p2b <- ggplot(tab) +
  annotate("rect", xmin = 0, xmax = 100, ymin = 40, ymax = 112,
           fill = COL$tint_blue, colour = NA) +
  annotate("rect", xmin = 0, xmax = 100, ymin = 0, ymax = 31,
           fill = COL$tint_amb, colour = NA) +
  annotate("text", x = 0, y = 114, hjust = 0, size = sz(6.7), fontface = "bold",
           colour = COL$ink, label = "cis-pQTL instrument available: 9 of 17 genes") +
  annotate("text", x = 0, y = 33.5, hjust = 0, size = sz(6.7), fontface = "bold",
           colour = COL$ink, label = "not testable at the protein layer: 8 of 17 genes") +
  geom_rect(aes(xmin = 0, xmax = bar, ymin = y - 3.6, ymax = y + 3.6, fill = fill),
            colour = NA) +
  scale_fill_identity() +
  geom_text(aes(x = bar + 2.2, y = y, label = lbl), hjust = 0, size = sz(7.4),
            fontface = "bold", colour = COL$ink) +
  geom_text(aes(x = 0, y = y + 8.4, label = short), hjust = 0, size = sz(6.3),
            colour = COL$ink2) +
  geom_text(aes(x = 0, y = y - 5.2, label = genes), hjust = 0, vjust = 1, size = sz(6.0),
            colour = COL$ink3, lineheight = 0.95) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 122), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(3, 3, 2, 3))
save_panel(p2b, file.path(outdir, "Fig2b_eqtl_pqtl_classification.pdf"), mm_snap(83.05), mm_snap(95))

# ============================================================
#  Figure 2c -- effect sizes in the two layers
# ============================================================
cat("\n>>> Figure 2c: effect-size comparison\n")
both <- con %>% filter(!is.na(b_SMR) & !is.na(pqtl_wald_b) & !is.na(pqtl_wald_p))
stopifnot(nrow(both) == 7)
ct <- cor.test(both$b_SMR, both$pqtl_wald_b, method = "spearman", use = "complete.obs")
stat <- sprintf("Spearman rho = %.3f, p = %.2f (n = %d)", ct$estimate, ct$p.value, nrow(both))

both <- both %>%
  mutate(sig = pqtl_wald_p < 0.05,
         dir = ifelse(direction_match, "concordant", "discordant"),
         rc  = ifelse(SYMBOL == "CCM2", "CCM2", "other"))

XR <- range(both$b_SMR); YR <- range(both$pqtl_wald_b)
xp <- diff(XR) * 0.16; yp <- diff(YR) * 0.10
XL <- c(XR[1] - xp, XR[2] + xp); YL <- c(YR[1] - yp, YR[2] + yp)

p2c <- ggplot(both) +
  # direction-concordant quadrants (I and III) carry a light tint
  annotate("rect", xmin = 0, xmax = XL[2], ymin = 0, ymax = YL[2],
           fill = COL$tint_blue, colour = NA) +
  annotate("rect", xmin = XL[1], xmax = 0, ymin = YL[1], ymax = 0,
           fill = COL$tint_blue, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = COL$ink3, linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL$ink3, linewidth = 0.45) +
  geom_point(aes(x = b_SMR, y = pqtl_wald_b, colour = dir, shape = sig,
                 size = rc), stroke = 0.55) +
  geom_text_repel(aes(x = b_SMR, y = pqtl_wald_b, label = SYMBOL,
                      fontface = ifelse(SYMBOL == "CCM2", "bold", "plain")),
                  size = sz(6.6), colour = COL$ink, max.overlaps = Inf, force = 12,
                  box.padding = 0.6, point.padding = 0.4, min.segment.length = 0,
                  seed = 5, segment.size = 0.25, segment.colour = COL$ink3) +
  scale_colour_manual(name = "Direction",
                      values = c(concordant = COL$genetic, discordant = COL$alert)) +
  scale_shape_manual(name = "pQTL Wald p",
                     values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "< 0.05", `FALSE` = ">= 0.05")) +
  scale_size_manual(name = NULL, values = c(CCM2 = 3.3, other = 2.5), guide = "none") +
  annotate("text", x = XL[1] + 0.05 * diff(XR), y = YL[2] - 0.03 * diff(YR), hjust = 0, vjust = 1, size = sz(6.4),
           colour = COL$ink2, lineheight = 1.15, label = paste0(
             stat, "\n", sum(both$direction_match), " of ", nrow(both),
             " direction-concordant\nshaded quadrants = direction-concordant")) +
  scale_x_continuous(limits = XL, expand = c(0, 0)) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  labs(x = "eQTL SMR effect size b (blood, per SD of expression)",
       y = "pQTL Wald ratio b (plasma, per SD of protein)") +
  th + theme(legend.position = "bottom", legend.box = "vertical",
             legend.box.just = "left", legend.spacing.y = unit(1.5, "pt"),
             plot.margin = margin(4, 8, 3, 4))
save_panel(p2c, file.path(outdir, "Fig2c_eqtl_pqtl_scatter.pdf"), mm_snap(83.05), mm_snap(95))

cat("\n== rev8 Figure 2 done ==\n")
