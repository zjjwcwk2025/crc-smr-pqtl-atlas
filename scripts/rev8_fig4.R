#!/usr/bin/env Rscript
# ============================================================
#  rev8_fig4.R -- Figure 4, redrawn to CNS art standard (rev8).
#  Panel order changed by decision: a NEW panel 4a (IDG tractability
#  gradient) is added; the former 4a (PheWAS) becomes 4b and the former
#  4b (drug repurposing) becomes 4c. Sizes of the two moved panels are
#  unchanged; only the new panel gets a new size.
#  Every number is read from product files and asserted with stopifnot().
#  Bars are drawn with geom_rect, never geom_col: geom_col guesses its
#  orientation from the mapped aesthetics and silently collapses when both
#  axes are continuous.
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2); library(dplyr); library(readr)
  library(tidyr); library(tibble); library(forcats); library(stringr)
  library(patchwork)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
source("scripts/theme_pub2.R")
outdir <- "results/figures_rev8_cns"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
th <- theme_pub2()

TIER_COL <- c(A = "#084594", B = "#2171B5", C = "#6BAED6", D = "#C6DBEF")

# ============================================================
#  Figure 4a -- IDG tractability of the six core candidates (NEW)
# ============================================================
cat("\n>>> Figure 4a: IDG tractability gradient\n")
idg <- fread("results/phase_tier1_idg_druggability.csv")
stopifnot(nrow(idg) == 6,
          setequal(idg$gene, c("UTP11","BMP2","PNKD","LAMC1","TMBIM1","GPBAR1")))
idg[, base := sub("\\*$", "", tdl)]
idg[, colr := c(Tclin = "#084594", Tchem = "#2171B5", Tbio = "#A6A6A6")[base]]
stopifnot(!any(is.na(idg$colr)))
idg[, yp := match(gene, idg[order(tdl_score, gene)]$gene)]     # 1 = bottom
idg[, lab := sprintf("%s   |   antibody entries: %d   |   %s", tdl, antibody, sm)]
XL4A <- c(0, 4.35)
YL4A <- c(0.45, max(idg$yp) + 1.20)

p4a <- ggplot(idg) +
  annotate("rect", xmin = XL4A[1], xmax = XL4A[2],
           ymin = idg[gene == "BMP2", yp] - 0.45, ymax = idg[gene == "BMP2", yp] + 0.45,
           fill = COL$tint_blue, colour = NA) +
  geom_rect(aes(xmin = 0, xmax = tdl_score, ymin = yp - 0.26, ymax = yp + 0.26,
                fill = colr), colour = NA) +
  scale_fill_identity() +
  geom_text(aes(x = tdl_score + 0.09, y = yp, label = sprintf("%.1f", tdl_score)),
            hjust = 0, size = sz(6.6), fontface = "bold", colour = COL$ink) +
  geom_text(aes(x = XL4A[2] + 0.30, y = yp, label = lab), hjust = 0,
            size = sz(6.2), colour = COL$ink2) +
  geom_text(aes(x = -0.12, y = yp, label = gene,
                fontface = ifelse(gene == "BMP2", "bold", "plain")),
            hjust = 1, size = sz(7.2), colour = COL$ink) +
  annotate("text", x = XL4A[2] + 0.30, y = max(idg$yp) + 0.80, hjust = 0,
           label = "IDG class   |   antibody tractability   |   structure / ligand evidence",
           size = sz(6.2), fontface = "bold", colour = COL$ink) +
  scale_x_continuous(expand = c(0, 0), breaks = 0:4) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = XL4A, ylim = YL4A, clip = "off") +
  labs(x = "IDG target development level score (Pharos)", y = NULL) +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.line.y = element_blank(), panel.grid.major.y = element_blank(),
             plot.margin = margin(4, 92, 3, 20, unit = "mm"))
save_panel(p4a, file.path(outdir, "Fig4a_idg_druggability.pdf"),
           mm_snap(169.69), mm_snap(58))

# ============================================================
#  Figure 4b -- PheWAS safety (was 4a, same size)
#  Cells are restricted to the organ-system risk categories, i.e. exactly the
#  categories whose counts add up to n_non_crc_high_risk (asserted below).
#  The previous version also carried "Other", "Cancer (non-CRC)", "CRC-related"
#  and "Infectious"; with those in the panel BMP2 no longer looked like the
#  narrowest row, although the caption says it is.
# ============================================================
cat("\n>>> Figure 4b: PheWAS safety\n")
ph <- fread("results/phase5c_phewas_safety.csv")
stopifnot(nrow(ph) == 6)
CATS <- c("Autoimmune/Inflammatory", "Cardiovascular", "Metabolic/Endocrine",
          "Neurological/Psychiatric", "Renal/Hepatic")
cells <- ph[, {
  kv <- strsplit(category_breakdown, "\\s*\\|\\s*")[[1]]
  data.table(cat = trimws(sub("=.*$", "", kv)),
             n   = as.numeric(sub("^.*=", "", kv)))
}, by = .(gene, n_non_crc_high_risk)]
cells <- cells[cat %in% CATS]
chk <- cells[, .(total = sum(n)), by = gene]
stopifnot(all(chk$total == ph$n_non_crc_high_risk[match(chk$gene, ph$gene)]))
cat("   row sums == n_non_crc_high_risk for all 6 genes\n")
grid <- merge(CJ(gene = ph$gene, cat = CATS, unique = TRUE), cells,
              by = c("gene", "cat"), all.x = TRUE)
grid[is.na(n), n := 0]
tot <- ph[, .(gene, total = n_non_crc_high_risk)]
tot[, yp := match(gene, ph[order(-n_non_crc_high_risk, gene)]$gene)]
grid[, yp := tot$yp[match(gene, tot$gene)]]
grid[, cat := factor(cat, levels = CATS)]
grid[, xp := as.numeric(cat)]
stopifnot(max(grid$n) == 6, min(tot$total) == 2, tot[gene == "BMP2", total] == 2)
YMAX <- max(tot$yp) + 1.05

p4b <- ggplot(grid) +
  annotate("rect", xmin = 0.5, xmax = length(CATS) + 0.5,
           ymin = tot[gene == "BMP2", yp] - 0.45, ymax = tot[gene == "BMP2", yp] + 0.45,
           fill = COL$tint_blue, colour = NA) +
  geom_tile(aes(x = xp, y = yp, fill = n), colour = "white", linewidth = 0.7) +
  geom_text(aes(x = xp, y = yp, label = ifelse(n > 0, n, ""), colour = n >= 3),
            size = sz(6.8), show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = COL$ink, `TRUE` = "white")) +
  geom_text(data = tot, aes(x = length(CATS) + 0.62, y = yp, label = total),
            hjust = 0, size = sz(6.8), colour = COL$ink,
            fontface = ifelse(tot$gene == "BMP2", "bold", "plain")) +
  geom_text(data = tot, aes(x = 0.62, y = yp, label = gene), hjust = 1,
            size = sz(7.0), colour = COL$ink,
            fontface = ifelse(tot$gene == "BMP2", "bold", "plain")) +
  annotate("text", x = 0.62, y = max(tot$yp) + 0.82, hjust = 1, label = "Gene",
           size = sz(6.4), fontface = "bold", colour = COL$ink) +
  annotate("text", x = length(CATS) + 0.62, y = max(tot$yp) + 0.82, hjust = 0,
           label = "High-risk\nphenotypes", size = sz(6.4), fontface = "bold",
           colour = COL$ink, lineheight = 0.95) +
  scale_fill_gradient(low = "#F2F7FC", high = "#08306B", breaks = c(0, 2, 4, 6),
                      name = "Genetic associations\n(blank = none)",
                      guide = guide_colourbar(direction = "horizontal",
                                              display = "rectangles",
                                              barwidth = unit(26, "mm"),
                                              barheight = unit(2.8, "mm"),
                                              title.position = "top")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq_along(CATS),
                     labels = str_replace(CATS, "/", "/\n")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, length(CATS) + 0.5), ylim = c(0.45, YMAX), clip = "off") +
  labs(x = NULL, y = NULL) +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.line.y = element_blank(), panel.grid = element_blank(),
             axis.ticks.x = element_blank(), axis.line.x = element_blank(),
             legend.position = "bottom", legend.title = element_text(size = 6.4),
             plot.margin = margin(6, 34, 3, 16, unit = "mm"))
save_panel(p4b, file.path(outdir, "Fig4b_phewas_safety.pdf"),
           mm_snap(169.69), mm_snap(86.08))

# ============================================================
#  Figure 4c -- drug repurposing (was 4b, same size)
#  Top strip: tier composition of the 62 protein-coding genes.
#  Bottom    : the same genes ranked by repurposing score, coloured by tier.
# ============================================================
cat("\n>>> Figure 4c: drug repurposing\n")
drug <- fread("results/phase5b_deep_drug_annotation.csv")
bio  <- fread("results/phase1_ensg_biotype.tsv")
coding <- bio[gene_type == "protein_coding", gene_name]
dpc <- drug[gene %in% coding]
stopifnot(nrow(dpc) == 62)
dpc[, letter := substr(drug_tier, 1, 1)]
dpc[, tname  := sub("^[A-D]\\s*\\(\\s*([A-Za-z]+).*$", "\\1", drug_tier)]
stopifnot(setequal(dpc$letter, c("A", "B", "C", "D")))
tnum <- dpc[, .N, by = .(letter, tname)][order(letter)]
stopifnot(sum(tnum$N) == 62, tnum[N == max(N), N] == 51)
LVL <- sprintf("%s %s (n = %d)", tnum$letter, tnum$tname, tnum$N)
tnum[, lvl := factor(LVL, levels = LVL)]
dpc  <- merge(dpc, tnum[, .(letter, lvl)], by = "letter", all.x = TRUE)
stopifnot(sum(dpc$n_drugs > 0) == 8, sum(dpc$n_approved > 0) == 4)
TCOL <- setNames(unname(TIER_COL[tnum$letter]), LVL)
cat(sprintf("   tiers: %s  (total %d); 8 (13%%) with agents, 4 approved\n",
            paste(LVL, collapse = " | "), sum(tnum$N)))

tnum[, xmax := cumsum(N)]
tnum[, xmin := c(0, head(xmax, -1))]
strip <- ggplot(tnum) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = 0.5, fill = lvl),
            colour = "white", linewidth = 0.6) +
  annotate("text", x = 13, y = 0, hjust = 0, size = sz(6.2), colour = COL$ink,
           label = "8 of 62 (13%) have known pharmacological agents; 4 with approved drugs") +
  annotate("text", x = -1.2, y = 0, hjust = 1, size = sz(6.4), colour = COL$ink,
           fontface = "bold", lineheight = 0.95,
           label = "62 protein-coding\nSMR-significant genes") +
  scale_fill_manual(name = NULL, values = TCOL) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 62), ylim = c(-1.05, 1.05), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "Liberation Sans") +
  theme(legend.position = "none",
        plot.margin = margin(2, 32, 2, 38, unit = "mm"))

CORE <- c("UTP11", "BMP2", "PNKD", "LAMC1", "TMBIM1", "GPBAR1")
top <- dpc[order(-repurposing_score)][1:25]
top[, yp := rev(seq_len(.N))]
top[, drug_lab := ifelse(n_drugs > 0, sprintf("%d (%d appr.)", n_drugs, n_approved), "-")]
stopifnot(all(CORE %in% top$gene))
Y25 <- max(top$yp)

rank <- ggplot(top) +
  geom_rect(aes(xmin = 0, xmax = repurposing_score, ymin = yp - 0.36, ymax = yp + 0.36,
                fill = lvl), colour = NA) +
  scale_fill_manual(name = NULL, values = TCOL) +
  geom_text(aes(x = repurposing_score + 0.010, y = yp,
                label = sprintf("%.3f", repurposing_score)),
            hjust = 0, size = sz(6.0), colour = COL$ink2) +
  geom_text(aes(x = -0.02, y = yp, label = gene,
                fontface = ifelse(gene %in% CORE, "bold", "plain")),
            hjust = 1, size = sz(6.4), colour = COL$ink) +
  geom_text(aes(x = 0.855, y = yp, label = drug_lab), hjust = 0, size = sz(6.0),
            colour = COL$ink2) +
  annotate("text", x = -0.02, y = Y25 + 0.95, hjust = 1, label = "Gene",
           size = sz(6.2), fontface = "bold", colour = COL$ink) +
  annotate("text", x = 0.855, y = Y25 + 0.95, hjust = 0,
           label = "agents (approved)", size = sz(6.2), fontface = "bold", colour = COL$ink) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 0.80), ylim = c(0.3, Y25 + 1.35), clip = "off") +
  labs(x = "Integrated repurposing score (components and weights in Methods)", y = NULL) +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.line.y = element_blank(), panel.grid.major.y = element_blank(),
             legend.position = "bottom",
             plot.margin = margin(4, 32, 3, 38, unit = "mm"))

p4c <- strip / rank + plot_layout(heights = c(1, 6.6))
save_panel(p4c, file.path(outdir, "Fig4c_drug_repurposing.pdf"),
           mm_snap(169.69), mm_snap(131.94))

cat("\n== rev8 Figure 4 done ==\n")
