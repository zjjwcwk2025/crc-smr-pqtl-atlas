#!/usr/bin/env Rscript
suppressMessages({ library(data.table); library(dplyr) })
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
tdir <- "manuscript/tables"

pfmt <- function(p, d = 2) {
  if (length(p) == 0 || is.na(p)) return("---")
  if (!is.finite(p) || p == 0) return("$<10^{-300}$")
  if (p >= 1e-3) return(sprintf("%.3f", p))
  s <- formatC(p, format = "e", digits = d)
  m <- regmatches(s, regexec("^([0-9.]+)e([+-]?)([0-9]+)$", s))[[1]]
  if (length(m) == 0) return(paste0("$", s, "$"))
  paste0("$", m[2], "\\times10^{", if (m[3] == "-") "-" else "", as.integer(m[4]), "}$")
}
pph4f <- function(x) {
  if (length(x) == 0 || is.na(x)) return("---")
  if (x >= 1e-3) return(sprintf("%.4f", x))
  pfmt(x)
}
bfmt <- function(b, d = 2) if (is.na(b)) "---" else sprintf(paste0("%.", d, "f"), b)

## ---------------- Table 1: six core non-MHC Tier 1 genes ----------------
sus  <- fread("results/phase3_susie_combined.csv")
top  <- fread("results/phase3_coloc_top10.csv")
saf  <- fread("results/phase5c_phewas_safety.csv")
drug <- fread("results/phase5b_deep_drug_annotation.csv")
six  <- c("UTP11","BMP2","PNKD","LAMC1","TMBIM1","GPBAR1")
t1 <- data.table(gene = six)
t1 <- merge(t1, sus[gene %in% six, .(gene, chr, PPH4, n_cs, top_pip)], by = "gene")
t1 <- merge(t1, top[gene %in% six, .(gene, p_SMR)], by = "gene")
t1 <- merge(t1, saf[gene %in% six, .(gene, unique_risk_categories, n_non_crc_high_risk, n_genetic_diseases)], by = "gene")
t1 <- merge(t1, drug[gene %in% six, .(gene, top_drug, top_stage, n_drugs_unique)], by = "gene")
setorder(t1, -PPH4)

rows <- vapply(seq_len(nrow(t1)), function(i) {
  d <- t1[i]
  dr <- if (is.na(d$n_drugs_unique) || d$n_drugs_unique == 0) "---" else sprintf("%s (%s)", d$top_drug, tolower(d$top_stage))
  sprintf("  %-6s & %s & %s & %d & %.2f & %s & %s & %d & %d \\\\",
          d$gene, d$chr, sprintf("%.4f", d$PPH4), d$n_cs, d$top_pip,
          pfmt(d$p_SMR), dr, d$unique_risk_categories, d$n_non_crc_high_risk)
}, character(1))

cat("\\begin{table}[htbp]\n\\centering\n", file = file.path(tdir, "Table1_tier1_genes.tex"))
cat("\\caption{The six core non-MHC Tier 1 genes with Bayesian colocalization (PPH4 $>$ 0.8) and SuSiE fine-mapping validation. Risk categories and high-risk phenotypes are from Phenome-Wide Association Study (PheWAS) screening in the Open Targets Platform (see Methods and Supplementary Table S6): risk categories = number of distinct organ-system categories carrying at least one genetic association for the gene; high-risk phenotypes = number of non-CRC phenotypes with genetic evidence in predefined high-risk categories.}\n", append = TRUE, file = file.path(tdir, "Table1_tier1_genes.tex"))
cat("\\label{tab:tier1}\n\\small\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{lllllllll}\n\\toprule\n", append = TRUE, file = file.path(tdir, "Table1_tier1_genes.tex"))
cat("  Gene & Chr & PPH4 & $n_{\\text{CS}}$ & PIP & $p_{\\text{SMR}}$ & Known drug & Risk categories & High-risk phenotypes \\\\\n\\midrule\n", append = TRUE, file = file.path(tdir, "Table1_tier1_genes.tex"))
cat(paste0(rows, collapse = "\n"), "\n", append = TRUE, file = file.path(tdir, "Table1_tier1_genes.tex"))
cat("\\bottomrule\n\\end{tabular}\n}\n\\end{table}\n", append = TRUE, file = file.path(tdir, "Table1_tier1_genes.tex"))
cat(">>> Table1 written\n"); print(t1[, .(gene, PPH4, unique_risk_categories, n_non_crc_high_risk, n_genetic_diseases)])

## ---------------- Table 2: cis-pQTL coverage ----------------
dec_sym <- unique(sub("^([0-9]+_[0-9]+_)([^_]+)_.*$", "\\2", list.files("data/decode", pattern = "\\.txt\\.gz$")))
ukb_sym <- unique(sub("_.*$", "", list.files("data/ukbppp_merged", pattern = "_merged\\.txt\\.gz$")))
mr_d <- fread("results/phase2_decode/phase2_decode_pqtl_mr_combined.csv")
mr_u <- fread("results/phase2_pqtl/phase2_pqtl_mr_combined.csv")
co_d <- fread("results/phase2_decode/phase2_decode_pqtl_coloc.csv")
co_u <- fread("results/phase3c_pqtl_coloc.csv")
fs   <- fread("results/fstatistics_summary.csv")

d2 <- data.table(gene = union(dec_sym, ukb_sym))
d2[, panel := ifelse(gene %in% dec_sym & gene %in% ukb_sym, "both",
              ifelse(gene %in% dec_sym, "deCODE", "UKB-PPP"))]
d2 <- merge(d2, mr_d[, .(gene, n_inst_d = n_instruments, wb_d = wald_b, wp_d = wald_pval)], by = "gene", all.x = TRUE)
d2 <- merge(d2, mr_u[, .(gene, n_inst_u = n_instruments, wb_u = b, wp_u = pval)], by = "gene", all.x = TRUE)
d2 <- merge(d2, co_d[, .(gene, pph4_d = PPH4)], by = "gene", all.x = TRUE)
d2 <- merge(d2, co_u[, .(gene, pph4_u = PPH4)], by = "gene", all.x = TRUE)
d2[, n_inst := ifelse(!is.na(n_inst_d), n_inst_d, n_inst_u)]
d2[, wb := ifelse(!is.na(wb_d), wb_d, wb_u)]
d2[, wp := ifelse(!is.na(wp_d), wp_d, wp_u)]
d2[, pph4 := ifelse(!is.na(pph4_d), pph4_d, pph4_u)]
d2[, fmed := fs$f_median[match(paste(gene, ifelse(!is.na(n_inst_d), "deCODE", "UKB-PPP")), paste(fs$gene, fs$source))]]
d2[, verdict := ifelse(!is.na(pph4) & pph4 > 0.8, "colocalizes",
                 ifelse(!is.na(n_inst) & n_inst > 0, "no coloc",
                 "no cis-pQTL instruments"))]
d2[, mr_panel := ifelse(!is.na(n_inst_d), "deCODE", ifelse(!is.na(n_inst_u), "UKB-PPP", "---"))]
d2[, orank := ifelse(verdict == "colocalizes", 1, ifelse(!is.na(n_inst), 2, 3))]
setorder(d2, orank, -pph4, gene)

rows2 <- vapply(seq_len(nrow(d2)), function(i) {
  d <- d2[i]
  ncol <- if (!is.na(d$n_inst)) as.character(d$n_inst) else "none"
  mrcol <- if (!is.na(d$wp)) sprintf("%+.3f (%s)", d$wb, pfmt(d$wp)) else "---"
  pcol  <- pph4f(d$pph4)
  fcol  <- if (!is.na(d$fmed)) sprintf("%.1f", d$fmed) else "---"
  sprintf("  %-9s & %-14s & %-8s & %-6s & %-6s & %-20s & %-10s & %s \\\\",
          d$gene, d$panel, d$mr_panel, ncol, fcol, mrcol, pcol, d$verdict)
}, character(1))
f2 <- file.path(tdir, "Table2_pqtl_coverage.tex")
cat("\\begin{table}[htbp]\n\\centering\n", file = f2)
cat("\\caption{Plasma protein coverage and pQTL--GWAS colocalization for the 17 genes with matching protein measurements. Cis-pQTL instruments were defined at $p < 5\\times10^{-8}$ within $\\pm$1 Mb of the gene; $n$ is the number of instruments used in the Mendelian randomization analysis. Wald-ratio estimates are the ratio of the CRC GWAS effect to the cis-pQTL effect on protein level. PPH4, posterior probability of a shared causal variant between the cis-pQTL and CRC GWAS signals (coloc.abf); PPH4 $>$ 0.8 indicates colocalization. Median $F$ is the median instrument $F$ statistic across the cis-pQTL instruments. No genome-wide significant cis-pQTL was detected for the seven deCODE-only proteins or for TET2; TET2 was nonetheless tested by colocalization and showed no evidence of a shared causal variant (PPH4 $= 0.006$). Panel names refer to the proteomic panel in which the protein was measured; MR panel is the panel supplying the cis-pQTL instruments.}\n", append = TRUE, file = f2)
cat("\\label{tab:pqtl_coverage}\n\\scriptsize\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{llllllll}\n\\toprule\n", append = TRUE, file = f2)
cat("  Gene & Panel & $n$ cis-pQTL & Median $F$ & Wald $b$ ($p$) & PPH4 & Protein-level support \\\\\n\\midrule\n", append = TRUE, file = f2)
cat(paste0(rows2, collapse = "\n"), "\n", append = TRUE, file = f2)
cat("\\bottomrule\n\\end{tabular}\n}\n\\end{table}\n", append = TRUE, file = f2)
cat(">>> Table2 written\n"); print(d2[, .(gene, panel, n_inst, wb, wp, pph4, verdict)])

## ---------------- Table 3: MR sensitivity ----------------
ms <- fread("results/phase5d_mr_methods_summary.csv")
ms[, cochran_p := pchisq(cochran_q, df = n_inst - 1, lower.tail = FALSE)]
pr <- if (file.exists("results/phase5d_mr_presso.csv")) fread("results/phase5d_mr_presso.csv") else data.table(gene = character(0))
ms[, c("presso_global_p", "presso_outlier_p", "presso_n_outliers") := NULL]
if (nrow(pr) > 0 && !all(c("nb_distribution", "presso_global_p_is_bound") %in% names(pr))) {
  cat(">>> NOTE: results/phase5d_mr_presso.csv predates the NbDistribution=5000 rerun; MR-PRESSO columns set to pending\n")
  pr <- data.table(gene = character(0))
}
if (nrow(pr) > 0) {
  if (!"presso_distortion_p" %in% names(pr)) pr[, presso_distortion_p := NA_real_]
  if (!"presso_distortion_p_is_bound" %in% names(pr)) pr[, presso_distortion_p_is_bound := NA_integer_]
  if (!"presso_corrected_b" %in% names(pr)) pr[, presso_corrected_b := NA_real_]
  if (!"presso_corrected_p" %in% names(pr)) pr[, presso_corrected_p := NA_real_]
  ms <- merge(ms, pr[, .(gene, presso_global_p, presso_global_p_is_bound, presso_n_outliers,
                         presso_distortion_p, presso_distortion_p_is_bound,
                         presso_corrected_b, presso_corrected_p,
                         presso_nb = nb_distribution, presso_note = note)], by = "gene", all.x = TRUE)
}
setorder(ms, n_inst)
getc <- function(nm, default) if (nm %in% names(ms)) ms[[nm]] else rep(default, nrow(ms))
gp  <- getc("presso_global_p", NA_real_)
gib <- getc("presso_global_p_is_bound", NA_integer_)
no  <- getc("presso_n_outliers", NA_integer_)
nt  <- as.character(getc("presso_note", ""))
dp  <- getc("presso_distortion_p", NA_real_)
dib <- getc("presso_distortion_p_is_bound", NA_integer_)
nbv <- getc("presso_nb", NA_integer_)
cb  <- getc("presso_corrected_b", NA_real_)
cp  <- getc("presso_corrected_p", NA_real_)

# MR-PRESSO folds the outlier p-value back through a Bonferroni-style factor nrow(data),
# so the achievable resolution is n/NbDistribution. Bounded p-values arrive as strings "<x".
pval_cell <- function(p, isb, nb) {
  if (is.na(p)) return("---")
  if (!is.na(isb) && isb == 1 && !is.na(nb) && nb > 0) {
    v <- 1 / nb
    if (v >= 1e-3) return(sprintf("$<%.3f$", v))
    s <- formatC(v, format = "e", digits = 1)
    m <- regmatches(s, regexec("^([0-9.]+)e([+-]?)([0-9]+)$", s))[[1]]
    if (length(m) == 0) return(sprintf("$<%s$", s))
    return(paste0("$<", m[2], "\\times10^{", if (m[3] == "-") "-" else "", as.integer(m[4]), "}$"))
  }
  pfmt(p)
}

rows3 <- vapply(seq_len(nrow(ms)), function(i) {
  d <- ms[i]
  gcell <- pval_cell(gp[i], gib[i], nbv[i])
  dcell <- pval_cell(dp[i], dib[i], nbv[i])
  ccell <- if (is.na(cb[i])) "---" else sprintf("$%+.3f$ (%s)", cb[i], pfmt(cp[i]))
  if (!is.na(nt[i]) && grepl("^skipped", nt[i])) {
    gcell <- "skipped$^{\\ddagger}$"; dcell <- "skipped$^{\\ddagger}$"; ccell <- "skipped$^{\\ddagger}$"
  }
  sprintf("  %-8s & %d & %+.3f (%s) & %+.3f (%s) & %+.3f (%s) & %+.4f (%s) & %s & %s & %s & %s & %s \\\\",
          d$gene, d$n_inst, d$ivw_b, pfmt(d$ivw_p), d$median_b, pfmt(d$median_p),
          d$modal_b, pfmt(d$modal_p), d$egger_int_b, pfmt(d$egger_int_p), pfmt(d$cochran_p),
          gcell, dcell, if (is.na(no[i])) "---" else as.character(no[i]), ccell)
}, character(1))

f3 <- file.path(tdir, "Table3_mr_sensitivity.tex")
cat("\\begin{table}[htbp]\n\\centering\n", file = f3)
cat("\\caption{Sensitivity analyses for the five genes with multi-instrument pQTL Mendelian randomization estimates. All estimates are on the log-odds scale per 1 SD of genetically predicted protein level. Cochran's $Q$ $p$-values test for heterogeneity across instruments; a non-significant MR-Egger intercept indicates no detectable directional pleiotropy. MR-PRESSO was run with \\texttt{NbDistribution = 5000} so that the outlier test resolves $n/\\mathrm{NbDistribution} \\le 0.046$ for every gene tested; at the package default of 1000 the same quantity is $0.061$--$0.228$, above the 0.05 significance threshold, so the default cannot identify outliers at this instrument count. The outlier-corrected column reports the IVW estimate after removing the instruments flagged by the MR-PRESSO outlier test; the Outliers column gives their number. Bounded $p$-values are reported as inequalities. $^{\\ddagger}$Not computationally feasible ($n = 2{,}288$ instruments).}\n", append = TRUE, file = f3)
cat("\\label{tab:mr_sensitivity}\n\\scriptsize\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{lllllllllll}\n\\toprule\n", append = TRUE, file = f3)
cat("  Gene & $n$ & IVW $b$ ($p$) & Weighted median $b$ ($p$) & Mode $b$ ($p$) & Egger intercept ($p$) & Cochran $Q$ $p$ & MR-PRESSO $p$ & MR-PRESSO distortion $p$ & Outliers & Outlier-corrected $b$ ($p$) \\\\\n\\midrule\n", append = TRUE, file = f3)
cat(paste0(rows3, collapse = "\n"), "\n", append = TRUE, file = f3)
cat("\\bottomrule\n\\end{tabular}\n}\n\\end{table}\n", append = TRUE, file = f3)
cat(">>> Table3 written\n")
print(ms[, .(gene, n_inst, ivw_b, ivw_p, egger_int_p, cochran_p)])
if ("presso_global_p" %in% names(ms)) print(ms[, .(gene, presso_nb, presso_global_p, presso_global_p_is_bound, presso_n_outliers, presso_distortion_p, presso_distortion_p_is_bound)])
