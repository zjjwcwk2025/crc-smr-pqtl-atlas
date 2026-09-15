#!/usr/bin/env Rscript
# 48_mr_methods_pleiotropy.R
# Group B — strengthen the pQTL MR with multi-method estimates and formal
# pleiotropy/heterogeneity diagnostics (all computed from instrument-level data):
#   FigS26_pqtl_mr_methods_forest.pdf   — per-gene multi-method comparison forest
#                                         (Wald, IVW, weighted median, weighted mode, MR-Egger)
#   FigS27_pqtl_mr_funnel_pleiotropy.pdf— per-gene funnel plots + MR-Egger intercept /
#                                         MR-PRESSO global test / Cochran Q annotations
#
# Uses the MendelianRandomization + MRPRESSO R packages (installed) on the same
# harmonised instrument CSVs used for the discovery pQTL MR. No hard-coded effects.

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(MendelianRandomization)
  library(MRPRESSO)
})

setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "results/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

theme_nc <- theme_classic(base_size = 9, base_family = "Liberation Sans") +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.text  = element_text(size = 8.5, color = "black"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

decode_sig <- c("CCM2", "LIMA1", "STAT6")
ukb_sig    <- c("CDKN1A", "TNF")
all_genes  <- c(decode_sig, ukb_sig)

load_instruments <- function(gene, source) {
  f <- if (source == "deCODE") sprintf("results/phase2_decode/%s_instruments.csv", gene)
       else sprintf("results/phase2_pqtl/%s_instruments.csv", gene)
  d <- fread(f, header = TRUE)
  d[, .(gene = gene, source = source,
        e_beta = pqtl.beta, e_se = pqtl.se,
        o_beta = gwas.beta, o_se = gwas.se)]
}

inst <- rbindlist(lapply(c(decode_sig, ukb_sig), function(g)
  load_instruments(g, if (g %in% decode_sig) "deCODE" else "UKB-PPP")))
inst <- inst[is.finite(e_beta) & is.finite(e_se) & e_beta != 0 &
               is.finite(o_beta) & is.finite(o_se)]

# ── compute per-SNP Wald ratio + delta-method SE (for funnel) ──────
inst[, ratio := o_beta / e_beta]
inst[, ratio_se := sqrt(o_se^2 / e_beta^2 + o_beta^2 * e_se^2 / e_beta^4)]
inst[, precision := 1 / ratio_se]

# ── run methods per gene using MendelianRandomization + MRPRESSO ────
run_mr <- function(dt) {
  bx <- dt$e_beta; bxse <- dt$e_se; by <- dt$o_beta; byse <- dt$o_se
  obj <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse)
  out <- list()
  ivw  <- mr_ivw(obj)
  out$ivw <- c(b = as.numeric(ivw@Estimate), se = as.numeric(ivw@StdError), p = as.numeric(ivw@Pvalue))
  egg  <- tryCatch(mr_egger(obj), error = function(e) NULL)
  if (!is.null(egg)) {
    out$egger <- c(b = as.numeric(egg@Estimate), se = as.numeric(egg@StdError.Est),
                   p = as.numeric(egg@Pvalue.Est))
    out$egger_int <- c(b = as.numeric(egg@Intercept), se = as.numeric(egg@StdError.Int),
                       p = as.numeric(egg@Pleio.pval))
  } else {
    out$egger <- c(b = NA, se = NA, p = NA)
    out$egger_int <- c(b = NA, se = NA, p = NA)
  }
  med  <- mr_median(obj)
  out$median <- c(b = as.numeric(med@Estimate), se = as.numeric(med@StdError), p = as.numeric(med@Pvalue))
  modal <- tryCatch(mr_mbe(obj), error = function(e) NULL)
  if (!is.null(modal)) out$modal <- c(b = as.numeric(modal@Estimate), se = as.numeric(modal@StdError),
                                      p = as.numeric(modal@Pvalue))
  else out$modal <- c(b = NA, se = NA, p = NA)
  # Cochran's Q from fixed-effect IVW
  r <- dt$ratio; se_r <- dt$ratio_se; bivw <- out$ivw["b"]
  out$cochran_q <- sum((r - bivw)^2 / se_r^2)
  out$n_inst <- nrow(dt)
  # MR-PRESSO (global + distortion / outlier)
  presso <- tryCatch(
    mr_presso(BetaOutcome = by, BetaExposure = bx, SdOutcome = byse,
              SdExposure = bxse, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
              NbDistribution = 1000, SignifThreshold = 0.05),
    error = function(e) NULL)
  if (!is.null(presso)) {
    gt <- presso$`MR-PRESSO results`$`Global Test`
    ot <- presso$`MR-PRESSO results`$`Outlier Test`
    pv <- suppressWarnings(as.numeric(ot$Pvalue))
    out$presso_global_p <- as.numeric(gt$Pvalue)
    out$presso_outlier_p <- if (length(pv) > 0) max(pv, na.rm = TRUE) else NA
    out$presso_n_outliers <- sum(pv < 0.05, na.rm = TRUE)
  } else {
    out$presso_global_p <- NA; out$presso_outlier_p <- NA; out$presso_n_outliers <- NA
  }
  out
}

res <- lapply(all_genes, function(g) {
  cat(sprintf("Running methods for %s ...\n", g))
  r <- run_mr(inst[gene == g])
  r$gene <- g
  r
})
names(res) <- all_genes

# Flatten each per-gene named vector (r$ivw, r$egger, ...) into a single
# 1-row-per-gene summary table. as.data.table on a list of named vectors would
# produce nested `ivw`/`egger` columns rather than `ivw.b` etc.; build explicitly.
summ_row <- function(g, r) {
  data.table(
    gene = g,
    ivw_b = as.numeric(r$ivw["b"]), ivw_se = as.numeric(r$ivw["se"]), ivw_p = as.numeric(r$ivw["p"]),
    egger_b = as.numeric(r$egger["b"]), egger_se = as.numeric(r$egger["se"]), egger_p = as.numeric(r$egger["p"]),
    egger_int_b = as.numeric(r$egger_int["b"]), egger_int_se = as.numeric(r$egger_int["se"]),
    egger_int_p = as.numeric(r$egger_int["p"]),
    median_b = as.numeric(r$median["b"]), median_se = as.numeric(r$median["se"]), median_p = as.numeric(r$median["p"]),
    modal_b = as.numeric(r$modal["b"]), modal_se = as.numeric(r$modal["se"]), modal_p = as.numeric(r$modal["p"]),
    cochran_q = r$cochran_q, n_inst = r$n_inst,
    presso_global_p = r$presso_global_p, presso_outlier_p = r$presso_outlier_p,
    presso_n_outliers = r$presso_n_outliers
  )
}
res_tbl <- rbindlist(lapply(all_genes, function(g) summ_row(g, res[[g]])))

cat("\n=== method estimates ===\n")
print(res_tbl[, .(gene, ivw_b, egger_b, egger_int = egger_int_b, egger_int_p,
                  median_b, modal_b, cochran_q, presso_global_p, presso_n_outliers)])

fwrite(res_tbl, "results/phase5d_mr_methods_summary.csv")
cat("\nSaved results/phase5d_mr_methods_summary.csv\n")

# ══════════════════════════════════════════════════════════════════
# FigS26 — multi-method comparison forest (one facet per gene)
# ══════════════════════════════════════════════════════════════════
melt_methods <- function(gene, r) {
  rbindlist(list(
    data.table(gene = gene, method = "IVW (fixed)",       b = r$ivw["b"],     se = r$ivw["se"]),
    data.table(gene = gene, method = "Weighted median",   b = r$median["b"],  se = r$median["se"]),
    data.table(gene = gene, method = "Weighted mode",     b = r$modal["b"],   se = r$modal["se"]),
    data.table(gene = gene, method = "MR-Egger (slope)",  b = r$egger["b"],   se = r$egger["se"])
  ))
}
forest <- rbindlist(lapply(all_genes, function(g) melt_methods(g, res[[g]])))
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
  labs(x = "Causal estimate on CRC risk (per SD), 95% CI", y = NULL,
       title = "pQTL MR — multi-method comparison",
       subtitle = "IVW / weighted median / weighted mode / MR-Egger (all instruments)") +
  theme_nc

ggsave(file.path(outdir, "FigS26_pqtl_mr_methods_forest.pdf"), p26,
       device = cairo_pdf, width = 6.30, height = 4.01)
cat("-> FigS26_pqtl_mr_methods_forest.pdf\n")

# ══════════════════════════════════════════════════════════════════
# FigS27 — funnel plots + pleiotropy annotations (one facet per gene)
# ══════════════════════════════════════════════════════════════════
funnel <- inst[, .(gene, ratio, precision)]
funnel[, gene := factor(gene, levels = all_genes)]

# annotation table per gene
annot <- res_tbl[, .(gene,
                     ivw_b, egger_b,
                     egg_int = egger_int_b, egg_int_p = egger_int_p,
                     cochran_q, presso_p = presso_global_p,
                     n_out = presso_n_outliers)]
annot[, gene := factor(gene, levels = all_genes)]
annot[, lab := sprintf("IVW = %.3f | Egger slope = %.3f\nEgger intercept = %.4f (p=%.3f)\nCochran Q = %.1f | MR-PRESSO global p = %.3f (outliers: %d)",
                       ivw_b, egger_b, egg_int, egg_int_p, cochran_q, presso_p,
                       ifelse(is.na(n_out), 0, n_out))]

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
  geom_text(data = annot, aes(x = Inf, y = Inf, label = lab),
            inherit.aes = FALSE, hjust = 1.02, vjust = 1.1, size = 2.5, colour = "grey25") +
  labs(x = "Per-SNP causal estimate (Wald ratio)", y = "Instrument precision (1/SE)",
       title = "pQTL MR — funnel plot & pleiotropy diagnostics",
       subtitle = "Symmetry around IVW line ~ no directional pleiotropy; Egger intercept & MR-PRESSO global p test pleiotropy") +
  theme_nc

ggsave(file.path(outdir, "FigS27_pqtl_mr_funnel_pleiotropy.pdf"), p27,
       device = cairo_pdf, width = 6.30, height = 4.01)
cat("-> FigS27_pqtl_mr_funnel_pleiotropy.pdf\n")

# ── postcondition checks ──────────────────────────────────────────
f26 <- file.path(outdir, "FigS26_pqtl_mr_methods_forest.pdf")
f27 <- file.path(outdir, "FigS27_pqtl_mr_funnel_pleiotropy.pdf")
stopifnot(file.exists(f26), file.exists(f27))
cat("Postcondition: FigS26 + FigS27 exist.\n")
