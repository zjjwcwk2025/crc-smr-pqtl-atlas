#!/usr/bin/env Rscript
# ============================================================
#  rev8_fig3.R -- Figure 3, redrawn to CNS art standard (rev8).
#  Panel sizes are identical to the previous approved round.
#  No number, threshold or claim is changed: every value is read
#  from the same product files the previous round used, and every
#  headline number is asserted with stopifnot().
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2); library(dplyr); library(readr)
  library(tidyr); library(tibble); library(forcats); library(stringr)
  library(ggrepel); library(scales); library(patchwork)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
source("scripts/theme_pub2.R")
outdir <- "results/figures_rev8_cns"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
th <- theme_pub2()

# ============================================================
#  Figure 3a -- convergence at CCM2
#  left  : pQTL effect on CRC risk (Wald ratio b, 95% CI) + numeric CI column
#  right : pQTL-GWAS colocalization (PPH4) on a log10 axis, so that the
#          four sub-threshold genes stay visible instead of collapsing
#          onto the axis. Same data, same threshold (0.8), same order.
# ============================================================
cat("\n>>> Figure 3a: convergence at CCM2\n")
mr <- fread("results/phase2_decode/phase2_decode_pqtl_mr_combined.csv")
co <- fread("results/phase2_decode/phase2_decode_pqtl_coloc.csv")
dd <- merge(mr[, .(gene, wald_b, wald_se, wald_pval, n_instruments)],
            co[, .(gene, PPH4)], by = "gene")
dd[, lo := wald_b - 1.96 * wald_se]
dd[, hi := wald_b + 1.96 * wald_se]
LV <- c("TNFRSF1B", "TNFRSF1A", "STAT6", "LIMA1", "CCM2")   # bottom -> top
dd[, yp := match(gene, LV)]
dd <- dd[order(yp)]
stopifnot(nrow(dd) == 5, !any(is.na(dd$yp)),
          sum(dd$PPH4 > 0.8) == 1,
          dd[gene == "CCM2", PPH4] > 0.8,
          dd[gene == "CCM2", wald_pval] < 0.05)
dd[, ci_lab := sprintf("%.3f\n(%.3f, %.3f)", wald_b, lo, hi)]

XR <- range(c(dd$lo, dd$hi)); XP <- diff(XR) * 0.07
XL <- c(XR[1] - XP, XR[2] + XP)
ycc <- dd[gene == "CCM2", yp]
YL <- c(0.50, 5.72)

pF <- ggplot(dd) +
  annotate("rect", xmin = XL[1], xmax = XL[2], ymin = ycc - 0.46, ymax = ycc + 0.46,
           fill = COL$tint_blue, colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL$ink3, linewidth = 0.45) +
  geom_errorbar(aes(xmin = lo, xmax = hi, y = yp), orientation = "y",
                width = 0.18, linewidth = 0.62, colour = COL$ink2) +
  geom_point(aes(x = wald_b, y = yp), colour = COL$genetic, size = 2.5) +
  geom_text(aes(x = XL[1] - 0.035 * diff(XR), y = yp, label = gene,
                fontface = ifelse(gene == "CCM2", "bold", "plain")),
            hjust = 1, size = sz(7.2), colour = COL$ink) +
  annotate("text", x = XL[2] + 0.055 * diff(XR), y = 5.62, hjust = 0, vjust = 1,
           label = "b (95% CI)", size = sz(6.2), fontface = "bold", colour = COL$ink) +
  geom_text(aes(x = XL[2] + 0.055 * diff(XR), y = yp, label = ci_lab),
            hjust = 0, size = sz(6.2), colour = COL$ink2, lineheight = 0.90) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  coord_cartesian(xlim = XL, clip = "off") +
  labs(x = "pQTL effect on CRC risk (Wald ratio b, 95% CI)", y = NULL) +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.line.y = element_blank(),
             plot.margin = margin(4, 36, 3, 15, unit = "mm"))

pR <- ggplot(dd) +
  annotate("rect", xmin = 1e-11, xmax = 2, ymin = ycc - 0.46, ymax = ycc + 0.46,
           fill = COL$tint_blue, colour = NA) +
  geom_segment(aes(x = 1e-11, xend = PPH4, y = yp, yend = yp),
               colour = COL$grid, linewidth = 1.6) +
  geom_vline(xintercept = 0.8, linetype = "dashed", colour = COL$alert, linewidth = 0.55) +
  geom_point(aes(x = PPH4, y = yp), colour = COL$genetic, size = 2.5) +
  scale_x_log10(limits = c(1e-11, 2),
                breaks = c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1),
                labels = c("1e-10", "1e-8", "1e-6", "1e-4", "0.01", "1")) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  labs(x = "pQTL-GWAS colocalization (PPH4, log10 scale)", y = NULL) +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.line.y = element_blank(), plot.margin = margin(4, 8, 3, 2, unit = "mm"))

p3a <- pF + pR + plot_layout(widths = c(1.45, 1))
save_panel(p3a, file.path(outdir, "Fig3a_ccm2_convergence.pdf"), mm_snap(169.69), mm_snap(64.56))

# ============================================================
#  Figure 3b -- blood vs colon SMR effect sizes
#  Colour marks the quantity the caption states (direction concordance,
#  21 of 25). Filled = concordant, open triangle = discordant, so the
#  encoding survives greyscale printing. y = x reference added.
# ============================================================
cat("\n>>> Figure 3b: blood vs colon SMR\n")
gtex <- fread("results/gtex_colon_transverse_smr.smr")
e1   <- fread("results/phase1_eqtl_smr.smr")
stopifnot(nrow(gtex) == 25, max(gtex$p_eQTL) < 5e-8)
m <- merge(gtex[, .(probeID, Gene, b_g = b_SMR, se_g = se_SMR, p_g = p_SMR)],
           e1[, .(probeID, b_b = b_SMR, se_b = se_SMR, p_b = p_SMR)], by = "probeID")
stopifnot(nrow(m) == 25)
ct <- cor.test(m$b_g, m$b_b, method = "spearman")
m[, conc := ifelse(sign(b_g) == sign(b_b), "Concordant", "Discordant")]
n_conc <- sum(m$conc == "Concordant")
stopifnot(abs(ct$estimate - 0.58) < 0.01,
          abs(ct$p.value - 0.0030) < 0.0006,
          n_conc == 21)
cat(sprintf("   rho = %.4f  p = %.4f  concordant = %d/%d\n",
            ct$estimate, ct$p.value, n_conc, nrow(m)))

XR <- range(m$b_b); YR <- range(m$b_g)
xp <- diff(XR) * 0.12; yp <- diff(YR) * 0.12
XL <- c(XR[1] - xp, XR[2] + xp); YL <- c(YR[1] - yp, YR[2] + yp)

stat_lab <- sprintf("Spearman rho = %.2f, p = %.3f\n%d of %d direction-concordant",
                    ct$estimate, ct$p.value, n_conc, nrow(m))

# ------------------------------------------------------------
#  Label seeding (2026-09-29).  In the approved layout five gene names
#  sat underneath the regression line: measured on the rendered PDF, up
#  to 18.7 % of a label's glyph ink overlapped the trend line.  ggrepel
#  cannot see a smooth as an obstacle, so those five are seeded away
#  from it.  force / box.padding / point.padding / max.overlaps / seed
#  are untouched, so the other twenty labels keep the approved layout.
#  After the change no label is cut by the trend line except B3GNTL1,
#  whose box still grazes it (4 % of its ink).
# ------------------------------------------------------------
NUDGE_X <- c("HLA-K" = -0.0289, "B3GNTL1" = -0.0340, "RP11-378A13.1" = -0.0554,
             "ZFP57" = 0.0339, "ARPC2" = -0.0161)
NUDGE_Y <- c("HLA-K" = 0.0931, "B3GNTL1" = -0.0127, "RP11-378A13.1" = 0.0284,
             "ZFP57" = 0.0172, "ARPC2" = 0.0281)
m[, `:=`(nx = 0, ny = 0)]
m[Gene %in% names(NUDGE_X), nx := NUDGE_X[Gene]]
m[Gene %in% names(NUDGE_Y), ny := NUDGE_Y[Gene]]
stopifnot(nrow(m) == 25, all(c("nx", "ny") %in% names(m)))

p3b <- ggplot(m) +
  annotate("rect", xmin = 0, xmax = XL[2], ymin = 0, ymax = YL[2],
           fill = COL$tint_blue, colour = NA) +
  annotate("rect", xmin = XL[1], xmax = 0, ymin = YL[1], ymax = 0,
           fill = COL$tint_blue, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = COL$ink3, linewidth = 0.40) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL$ink3, linewidth = 0.40) +
  geom_abline(slope = 1, intercept = 0, colour = COL$ink3, linewidth = 0.45, linetype = "dotted") +
  geom_smooth(aes(x = b_b, y = b_g), method = "lm", formula = y ~ x,
              colour = COL$ink2, fill = COL$mute, alpha = 0.25, linewidth = 0.60) +
  geom_point(aes(x = b_b, y = b_g, colour = conc, shape = conc), size = 2.6, stroke = 0.65) +
  geom_text_repel(aes(x = b_b, y = b_g, label = Gene), size = sz(6.4), colour = COL$ink,
                  max.overlaps = Inf, force = 14, box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0, seed = as.integer(Sys.getenv("FIG3B_SEED", "11")), segment.size = 0.22,
                  nudge_x = m$nx, nudge_y = m$ny,
                  segment.colour = COL$ink3) +
  annotate("label", x = XL[1] + 0.025 * diff(XR), y = YL[2] - 0.02 * diff(YR),
           hjust = 0, vjust = 1, label = stat_lab, size = sz(6.4), colour = COL$ink,
           fill = "white", linewidth = 0.20, label.r = unit(1.2, "pt"), lineheight = 1.10) +
  scale_colour_manual(name = NULL,
                      values = c(Concordant = COL$genetic, Discordant = COL$alert)) +
  scale_shape_manual(name = NULL,
                     values = c(Concordant = 16, Discordant = 2)) +
  scale_x_continuous(limits = XL, expand = c(0, 0)) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  labs(x = "SMR effect size b (eQTLGen blood, per SD of expression)",
       y = "SMR effect size b (GTEx colon transverse)") +
  th + theme(legend.position = "bottom", plot.margin = margin(4, 8, 3, 4))
save_panel(p3b, file.path(outdir, "Fig3b_gtex_colon_scatter.pdf"), mm_snap(169.69), mm_snap(134.76))

# ============================================================
#  Figure 3c -- forest of the 10 genes significant in both tissues
#  colour + shape both encode the tissue, so the panel reads in
#  greyscale; legend moved out of the data area.
# ============================================================
cat("\n>>> Figure 3c: dual-significant forest\n")
dual <- m[p_g < 3.21e-6 & p_b < 3.21e-6]
stopifnot(nrow(dual) == 10)
long <- rbind(
  data.table(gene = dual$Gene, tissue = "eQTLGen blood", b = dual$b_b, se = dual$se_b),
  data.table(gene = dual$Gene, tissue = "GTEx colon",   b = dual$b_g, se = dual$se_g))
long[, lo := b - 1.96 * se][, hi := b + 1.96 * se]
ordg <- long[, .(mb = mean(b)), by = gene][order(mb), gene]
long[, gene := factor(gene, levels = ordg)]

p3c <- ggplot(long) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL$ink3, linewidth = 0.45) +
  geom_errorbar(aes(xmin = lo, xmax = hi, y = gene, colour = tissue),
                orientation = "y", width = 0.30, linewidth = 0.58,
                position = position_dodge(width = 0.62)) +
  geom_point(aes(x = b, y = gene, colour = tissue, shape = tissue),
             size = 2.5, stroke = 0.65, position = position_dodge(width = 0.62)) +
  scale_colour_manual(name = NULL,
                      values = c(`eQTLGen blood` = COL$genetic, `GTEx colon` = "#6BAED6")) +
  scale_shape_manual(name = NULL,
                     values = c(`eQTLGen blood` = 16, `GTEx colon` = 17)) +
  labs(x = "SMR effect size b (95% CI)", y = NULL) +
  th + theme(legend.position = "bottom", legend.box = "horizontal",
             plot.margin = margin(4, 8, 3, 4))
save_panel(p3c, file.path(outdir, "Fig3c_gtex_colon_forest.pdf"), mm_snap(169.69), mm_snap(113.24))

cat("\n== rev8 Figure 3 done ==\n")
