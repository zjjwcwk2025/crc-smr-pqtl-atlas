#!/usr/bin/env Rscript
# figS26S27_regen.R
# Regenerate FigS26 and FigS27 only, carrying three fixes:
#   1. FigS26 x axis: drop tick labels that would collide (guide_axis check.overlap)
#   2. FigS27 MR-PRESSO annotation: read the NbDistribution = 5000 values from
#      results/phase5d_mr_presso.csv. The inline MR-PRESSO call in
#      48_mr_methods_pleiotropy.R runs at the package default (1000), which
#      cannot resolve n/NbDistribution <= 0.046 at these instrument counts and
#      returns NULL for every gene, so the figure previously read
#      "MR-PRESSO p = NA (0 outliers)" while its caption reported
#      5/10/15/22 outliers.
#   3. FigS27 MR-Egger intercept p value: print "< 0.001" instead of "0.000".
# 48_mr_methods_pleiotropy.R is NOT re-run because it would recompute MR-PRESSO
# for TNF (2,288 instruments), which is not computationally feasible; the plots
# are rebuilt from the cached results/phase5d_mr_methods_summary.csv.
suppressMessages({ library(data.table); library(ggplot2); library(dplyr); library(tidyr) })
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "results/figures"
theme_nc <- theme_classic(base_size = 9, base_family = "Liberation Sans") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8.5, color = "black"),
        plot.title = element_blank(), plot.subtitle = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "none")
decode_sig <- c("CCM2", "LIMA1", "STAT6")
ukb_sig <- c("CDKN1A", "TNF")
all_genes <- c(decode_sig, ukb_sig)

res_tbl <- as.data.table(fread("results/phase5d_mr_methods_summary.csv"))
pr <- as.data.table(fread("results/phase5d_mr_presso.csv"))
drop_cols <- intersect(c("presso_global_p", "presso_outlier_p", "presso_n_outliers"), names(res_tbl))
if (length(drop_cols)) res_tbl[, (drop_cols) := NULL]
prx <- pr[, c("gene", "presso_global_p", "presso_n_outliers", "presso_global_p_is_bound", "n_inst"), with = FALSE]
    setnames(prx, "n_inst", "presso_n_inst")
    res_tbl <- merge(res_tbl, prx, by = "gene", all.x = TRUE)

## --- FigS26: multi-method comparison forest -----------------------
forest <- rbindlist(lapply(all_genes, function(g) {
  r <- res_tbl[gene == g]
  rbindlist(list(
    data.table(gene = g, method = "IVW (fixed)", b = r$ivw_b, se = r$ivw_se),
    data.table(gene = g, method = "Weighted median", b = r$median_b, se = r$median_se),
    data.table(gene = g, method = "Weighted mode", b = r$modal_b, se = r$modal_se),
    data.table(gene = g, method = "MR-Egger (slope)", b = r$egger_b, se = r$egger_se)))
}))
forest[, lo := b - 1.96 * se][, hi := b + 1.96 * se]
forest[, method := factor(method, levels = rev(c("IVW (fixed)", "Weighted median",
                                                 "Weighted mode", "MR-Egger (slope)")))]
forest[, gene := factor(gene, levels = all_genes)]
p26 <- ggplot(forest, aes(x = b, y = method, color = method)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4,
                  linewidth = 0.6, position = position_dodge(width = 0.4)) +
  facet_wrap(~ gene, scales = "free", ncol = 3) +
  scale_color_manual(values = c("IVW (fixed)" = "#2166AC", "Weighted median" = "#1B7837",
                                "Weighted mode" = "#D6604D", "MR-Egger (slope)" = "#B2182B")) +
  scale_x_continuous(n.breaks = 4, guide = guide_axis(check.overlap = TRUE)) +
  labs(x = "Causal estimate on CRC risk (per SD), 95% CI", y = NULL) +
  theme_nc
ggsave(file.path(outdir, "FigS26_pqtl_mr_methods_forest.pdf"), p26,
       device = cairo_pdf, width = 6.30, height = 4.01)
cat("-> FigS26_pqtl_mr_methods_forest.pdf\n")

## --- FigS27: funnel + pleiotropy annotations ----------------------
load_instruments <- function(gene, source) {
  f <- if (source == "deCODE") sprintf("results/phase2_decode/%s_instruments.csv", gene)
       else sprintf("results/phase2_pqtl/%s_instruments.csv", gene)
  d <- fread(f, header = TRUE)
  d[, .(gene = gene, source = source, e_beta = pqtl.beta, e_se = pqtl.se,
        o_beta = gwas.beta, o_se = gwas.se)]
}
inst <- rbindlist(lapply(all_genes, function(g)
  load_instruments(g, if (g %in% decode_sig) "deCODE" else "UKB-PPP")))
inst <- inst[is.finite(e_beta) & is.finite(e_se) & e_beta != 0 &
               is.finite(o_beta) & is.finite(o_se)]
inst[, ratio := o_beta / e_beta]
inst[, ratio_se := sqrt(o_se^2 / e_beta^2 + o_beta^2 * e_se^2 / e_beta^4)]
inst[, precision := 1 / ratio_se]
funnel <- inst[, .(gene, ratio, precision)]
funnel[, gene := factor(gene, levels = all_genes)]

annot <- res_tbl[, .(gene, ivw_b, egger_b,
                     egg_int = egger_int_b, egg_int_p = egger_int_p,
                     cochran_q,
                     presso_p = presso_global_p,
                     p_bound = presso_global_p_is_bound,
                     n_out = presso_n_outliers,
                     n_inst = presso_n_inst)]
annot[, gene := factor(gene, levels = all_genes)]
# Two lines for MR-PRESSO: the widest single line must stay inside the facet
# (a 39-character line was clipped at the panel edge).
annot[, press_line := ifelse(is.na(n_inst), "MR-PRESSO p = NA",
                      ifelse(n_inst > 300 & is.na(presso_p),
                             sprintf("MR-PRESSO p = not feasible\n(n = %d instruments)", n_inst),
                             sprintf("MR-PRESSO p = %s\n(%d outlier instruments)",
                                     ifelse(!is.na(p_bound) & p_bound == 1, "<0.001",
                                            sprintf("%.3f", presso_p)), n_out)))]
annot[, lab := sprintf("IVW = %.3f\nEgger slope = %.3f\nEgger int. = %.4f (p %s)\nCochran Q = %.1f\n%s",
                       ivw_b, egger_b, egg_int,
                       ifelse(egg_int_p < 0.001, "< 0.001", sprintf("= %.3f", egg_int_p)),
                       cochran_q, press_line)]
print(annot[, .(gene, ivw_b, egger_b, egg_int_p, cochran_q, press_line)])

p27 <- ggplot(funnel, aes(x = ratio, y = precision)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey70") +
  geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
  geom_point(alpha = 0.5, size = 1.1, color = "#377EB8") +
  geom_vline(aes(xintercept = ivw_b, color = "IVW"), data = annot,
             linetype = 2, linewidth = 0.8) +
  geom_vline(aes(xintercept = egger_b, color = "MR-Egger slope"), data = annot,
             linetype = 2, linewidth = 0.8) +
  facet_wrap(~ gene, scales = "free", ncol = 3) +
  scale_color_manual(values = c("IVW" = "#2166AC", "MR-Egger slope" = "#B2182B"),
                     name = NULL, breaks = c("IVW", "MR-Egger slope")) +
  geom_label(data = annot, aes(x = Inf, y = Inf, label = lab),
             inherit.aes = FALSE, hjust = 1.04, vjust = 1.06, size = 2.5,
             colour = "grey20", fill = "white", alpha = 0.72,
             label.size = 0, label.padding = unit(0.08, "lines")) +
  scale_x_continuous(n.breaks = 4) +
  labs(x = "Per-SNP causal estimate (Wald ratio)", y = "Instrument precision (1/SE)") +
  theme_nc + theme(legend.position = "bottom",
                   legend.key.width = unit(0.75, "cm"),
                   legend.key.height = unit(0.22, "cm"))
ggsave(file.path(outdir, "FigS27_pqtl_mr_funnel_pleiotropy.pdf"), p27,
       device = cairo_pdf, width = 6.30, height = 4.01)
cat("-> FigS27_pqtl_mr_funnel_pleiotropy.pdf\n")
