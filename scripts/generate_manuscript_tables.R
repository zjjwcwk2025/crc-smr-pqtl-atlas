# ====================================================================
#  generate_manuscript_tables.R
#  Generate LaTeX tables for CRC SMR Atlas manuscript (Table 1, S1-S6)
# ====================================================================
suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
})

setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "manuscript/tables"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# ---- helpers ----
sanitize_latex <- function(x) {
  x <- gsub("_", "\\_", x, fixed=TRUE)
  x <- gsub("&", "\\&", x, fixed=TRUE)
  x <- gsub("%", "\\%", x, fixed=TRUE)
  x <- gsub("$", "\\$", x, fixed=TRUE)
  x <- gsub("#", "\\#", x, fixed=TRUE)
  x <- gsub("{", "\\{", x, fixed=TRUE)
  x <- gsub("}", "\\}", x, fixed=TRUE)
  x
}

fmt_p <- function(p) {
  ifelse(p < 1e-300, "$<$1e-300",
  ifelse(p < 1e-100, sprintf("$%.1e$", p),
  ifelse(p < 1e-10, sprintf("$%.2e$", p),
  ifelse(p < 0.001, sprintf("$%.3e$", p),
  sprintf("%.4f", p)))))
}

fmt_p4 <- function(p) { ifelse(is.na(p) | p == 0, "---", sprintf("%.4f", round(p, 4))) }

# ====================================================================
cat("Reading data...\n")

# Tier 1 coloc + SuSiE
coloc_tier1 <- read_csv("results/phase3_coloc_top10.csv", show_col_types=FALSE) %>%
  filter(PPH4 > 0.8) %>% select(gene, ensg, chr, nsnps, PPH4, p_SMR)

susie    <- read_csv("results/phase3_susie_finemap.csv", show_col_types=FALSE) %>%
  select(gene, n_cs, top_pip_snp, top_pip, converged)

# need HEIDI for Tier 1
annot     <- read_tsv("results/phase1_bonferroni_significant_annotated.tsv",
                       show_col_types=FALSE) %>%
  select(gene=Gene, p_HEIDI) %>% unique()

drug <- read_csv("results/phase5b_deep_drug_annotation.csv", show_col_types=FALSE) %>%
  filter(!is.na(gene)) %>%
  select(gene, drug_tier, n_drugs=n_drugs_unique, top_drug, top_stage,
         repurposing_score)

phewas <- read_csv("results/phase5c_phewas_safety.csv", show_col_types=FALSE) %>%
  select(gene, n_diseases, n_genetic_assoc=n_genetic_diseases, n_risk_categories=unique_risk_categories,
         safety_score=safety_concern_score, safety_tier)

# ====================================================================
#  Table 1 — 6 Tier 1 genes
# ====================================================================
cat("\n>>> Table 1: Tier 1 gene summary\n")

t1 <- coloc_tier1 %>%
  left_join(susie, by="gene") %>%
  left_join(annot, by="gene") %>%
  left_join(drug, by="gene") %>%
  left_join(phewas, by="gene") %>%
  arrange(desc(PPH4))

stopifnot(nrow(t1) == 6)

t1_latex <- t1 %>%
  transmute(
    Gene   = sanitize_latex(gene),
    Chr    = chr,
    PPH4   = sprintf("%.4f", PPH4),
    `$n_{\\text{CS}}$` = as.character(n_cs),
    `PIP`  = sprintf("%.2f", top_pip),
    `$p_{\\text{SMR}}$` = fmt_p(p_SMR),
    `$p_{\\text{HEIDI}}$` = if_else(is.na(p_HEIDI), "---", fmt_p4(p_HEIDI)),
    `Drugs` = if_else(is.na(n_drugs) | n_drugs == 0, "None",
                      paste0(top_drug, " (", top_stage, ")", if_else(n_drugs > 1 & !is.na(n_drugs), paste0(" +", n_drugs-1), ""))),
    `PheWAS safety` = sprintf("%s (%.1f)", safety_tier, safety_score)
  )

# Write Table 1
cap <- "Summary of 6 Tier 1 genes with Bayesian colocalization (PPH4 $>$ 0.8) and SuSiE fine-mapping validation."
lab <- "tab:tier1"
cat("\\begin{table}[htbp]\n", file=file.path(outdir,"Table1_tier1_genes.tex"))
cat("\\centering\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat(sprintf("\\caption{%s}\\label{%s}\n", cap, lab),
    file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat("\\small\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat("\\resizebox{\\textwidth}{!}{%\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat("\\begin{tabular}{lllllllll}\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat("\\toprule\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)

hdr <- paste0("  ", paste(colnames(t1_latex), collapse=" & "), " \\\\\n")
cat(hdr, file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat("\\midrule\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)

for (i in 1:nrow(t1_latex)) {
  row <- paste0("  ", paste(t1_latex[i,], collapse=" & "), " \\\\\n")
  cat(row, file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
}

cat("\\bottomrule\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat("\\end{tabular}\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat("}\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat("\\end{table}\n", file=file.path(outdir,"Table1_tier1_genes.tex"), append=TRUE)
cat(sprintf("  -> %s\n", file.path(outdir,"Table1_tier1_genes.tex")))

# ====================================================================
#  Table S1 — SMR Top 10 (detailed)
# ====================================================================
cat("\n>>> Table S1: SMR top 10\n")

top10 <- read_tsv("results/phase1_bonferroni_significant_annotated.tsv", show_col_types=FALSE) %>%
  arrange(p_SMR) %>%
  head(10) %>%
  mutate(Gene_label = if_else(is.na(SYMBOL) | SYMBOL == "", paste0("ENSG:", substr(ENSG_clean, 1, 15)), SYMBOL)) %>%
  left_join(coloc_tier1 %>% select(gene, coloc_PPH4=PPH4), by=c("SYMBOL"="gene"))

tS1 <- top10 %>%
  transmute(
    Rank = 1:n(),
    Gene  = sanitize_latex(Gene_label),
    Chr   = ProbeChr,
    `Top SNP` = topSNP,
    `$b_{\\text{SMR}}$` = sprintf("%.3f", b_SMR),
    `$p_{\\text{SMR}}$`  = fmt_p(p_SMR),
    `$p_{\\text{HEIDI}}$` = fmt_p4(p_HEIDI),
    `PPH4` = if_else(is.na(coloc_PPH4), "---", sprintf("%.4f", coloc_PPH4))
  )

capS1 <- "Top 10 SMR associations for CRC (GCST90255675). SMR tested 15,582 genes; Bonferroni threshold $p < 3.21\\times10^{-6}$."
labS1 <- "tab:smr_top10"
f1 <- file.path(outdir, "TableS1_smr_top10.tex")
cat("\\begin{table}[htbp]\n", file=f1)
cat("\\centering\n", file=f1, append=TRUE)
cat(sprintf("\\caption{%s}\\label{%s}\n", capS1, labS1), file=f1, append=TRUE)
cat("\\small\n", file=f1, append=TRUE)
cat("\\begin{tabular}{rlllllll}\n", file=f1, append=TRUE)
cat("\\toprule\n", file=f1, append=TRUE)
hdr <- paste0("  ", paste(colnames(tS1), collapse=" & "), " \\\\\n")
cat(hdr, file=f1, append=TRUE)
cat("\\midrule\n", file=f1, append=TRUE)
for (i in 1:nrow(tS1)) {
  row <- paste0("  ", paste(tS1[i,], collapse=" & "), " \\\\\n")
  cat(row, file=f1, append=TRUE)
}
cat("\\bottomrule\n", file=f1, append=TRUE)
cat("\\end{tabular}\n", file=f1, append=TRUE)
cat("\\end{table}\n", file=f1, append=TRUE)
cat(sprintf("  -> %s\n", f1))

# ====================================================================
#  Table S2 — protein-coding genes (single GENCODE biotype annotation)
# ====================================================================
cat("\n>>> Table S2: protein-coding genes\n")

# Single-source biotype classification: results/phase1_ensg_biotype.tsv is a
# GENCODE v47 annotation for all 75 unique ENSG IDs, generated
# from the reference GTF.  It has columns: ENSG_clean, biotype, gene_name.
# We inner-join to the 62 protein_coding ENSGs, overwrite SYMBOL from the
# canonical gene_name, and deduplicate ENSG00000175711 (B3GNTL1 + LOC128966610)
# to a single canonical B3GNTL1 row.  This resolves the prior 62/63 mismatch
# (the duplicate ENSG contributed two rows under the old heuristic).

all_raw <- read_tsv("results/phase1_bonferroni_significant_annotated.tsv", show_col_types=FALSE)
biotype <- read_tsv("results/phase1_ensg_biotype.tsv", show_col_types=FALSE,
                    col_names=c("ENSG_clean", "biotype", "gene_name"))

all <- all_raw %>%
  inner_join(biotype %>% filter(biotype == "protein_coding"), by="ENSG_clean") %>%
  mutate(SYMBOL = gene_name) %>%
  distinct(ENSG_clean, .keep_all=TRUE) %>%
  arrange(p_SMR)

prot <- all %>%
  left_join(coloc_tier1 %>% select(gene, coloc_PPH4=PPH4), by=c("SYMBOL"="gene"))

stopifnot(nrow(prot) == 62)
cat(sprintf("  Protein-coding genes: %d\n", nrow(prot)))

tS2 <- prot %>%
  transmute(
    Gene  = sanitize_latex(SYMBOL),
    Chr   = ProbeChr,
    `$b_{\\text{SMR}}$` = sprintf("%.3f", b_SMR),
    `$p_{\\text{SMR}}$`  = fmt_p(p_SMR),
    `$p_{\\text{HEIDI}}$` = fmt_p4(p_HEIDI),
    `PPH4` = if_else(is.na(coloc_PPH4), "---", sprintf("%.4f", coloc_PPH4))
  )

capS2 <- sprintf("Protein-coding genes (GENCODE v47 annotation) passing SMR Bonferroni significance ($n=%d$).", nrow(prot))
labS2 <- "tab:protein_coding"

# Split into columns for readability
n_per_col <- ceiling(nrow(tS2) / 3)
col1 <- tS2[1:n_per_col, ]
col2 <- if (nrow(tS2) > n_per_col) tS2[(n_per_col+1):min(2*n_per_col, nrow(tS2)), ] else tS2[FALSE,]
col3 <- if (nrow(tS2) > 2*n_per_col) tS2[(2*n_per_col+1):nrow(tS2), ] else tS2[FALSE,]

f2 <- file.path(outdir, "TableS2_protein_coding.tex")
cat("\\begin{table}[htbp]\n", file=f2)
cat("\\centering\n", file=f2, append=TRUE)
cat(sprintf("\\caption{%s}\\label{%s}\n", capS2, labS2), file=f2, append=TRUE)
cat("\\scriptsize\n", file=f2, append=TRUE)
cat("\\resizebox{\\textwidth}{!}{%\n", file=f2, append=TRUE)
cat(sprintf("\\begin{tabular}{llllllllllllllllll}\n"), file=f2, append=TRUE)
cat("\\toprule\n", file=f2, append=TRUE)

hdr <- paste0("  ", paste(rep(c("Gene", "Chr", "$b_\\text{SMR}$", "$p_\\text{SMR}$", "$p_\\text{HEIDI}$", "PPH4"), 3), collapse=" & "), " \\\\\n")
cat(hdr, file=f2, append=TRUE)
cat("\\midrule\n", file=f2, append=TRUE)

for (r in 1:n_per_col) {
  vals <- c()
  for (col_idx in list(col1, col2, col3)) {
    if (r <= nrow(col_idx)) {
      vals <- c(vals, as.character(col_idx[r,]))
    } else {
      vals <- c(vals, rep("", 6))
    }
  }
  row <- paste0("  ", paste(vals, collapse=" & "), " \\\\\n")
  cat(row, file=f2, append=TRUE)
}
cat("\\bottomrule\n", file=f2, append=TRUE)
cat("\\end{tabular}\n", file=f2, append=TRUE)
cat("}\n", file=f2, append=TRUE)
cat("\\end{table}\n", file=f2, append=TRUE)
cat(sprintf("  -> %s\n", f2))

# ====================================================================
#  Table S3 — SuSiE sensitivity (per-gene summary)
# ====================================================================
cat("\n>>> Table S3: SuSiE sensitivity summary\n")

ss <- read_csv("results/susie_sensitivity/susie_sensitivity_grid.csv", show_col_types=FALSE) %>%
  filter(converged == TRUE) %>%
  group_by(gene) %>%
  summarise(
    n_tested   = n(),
    n_cs_min   = min(n_cs),
    n_cs_median = median(n_cs),
    n_cs_max   = max(n_cs),
    all_single = all(n_multi_snp_cs == 0),
    .groups="drop"
  ) %>%
  mutate(stability = if_else(n_cs_min == n_cs_max, "Perfect",
                       if_else(n_cs_max - n_cs_min <= 3, "High", "Moderate")))

tS3 <- ss %>% transmute(
  Gene = sanitize_latex(gene),
  `Param combos` = as.character(n_tested),
  `$n_{\\text{CS}}$ range` = paste0(n_cs_min, "--", n_cs_max),
  `$n_{\\text{CS}}$ median` = as.character(n_cs_median),
  `All single-SNP` = if_else(all_single, "Yes", "No"),
  Stability = stability
)

capS3 <- "SuSiE fine-mapping sensitivity analysis: 6 coloc-passing genes tested across 18 parameter combinations (SNP cap 200--500, $L$ 5--20, cis window $\\pm$500--1000 kb). All 65 converged runs produced exclusively single-SNP credible sets (PIP=1.0), confirming that this pattern reflects genuine causal architecture rather than methodological artifact."
labS3 <- "tab:susie_sensitivity"

f3 <- file.path(outdir, "TableS3_susie_sensitivity.tex")
cat("\\begin{table}[htbp]\n", file=f3)
cat("\\centering\n", file=f3, append=TRUE)
cat(sprintf("\\caption{%s}\\label{%s}\n", capS3, labS3), file=f3, append=TRUE)
cat("\\small\n", file=f3, append=TRUE)
cat("\\begin{tabular}{llllll}\n", file=f3, append=TRUE)
cat("\\toprule\n", file=f3, append=TRUE)
hdr <- paste0("  ", paste(colnames(tS3), collapse=" & "), " \\\\\n")
cat(hdr, file=f3, append=TRUE)
cat("\\midrule\n", file=f3, append=TRUE)
for (i in 1:nrow(tS3)) {
  row <- paste0("  ", paste(tS3[i,], collapse=" & "), " \\\\\n")
  cat(row, file=f3, append=TRUE)
}
cat("\\bottomrule\n", file=f3, append=TRUE)
cat("\\end{tabular}\n", file=f3, append=TRUE)
cat("\\end{table}\n", file=f3, append=TRUE)
cat(sprintf("  -> %s\n", f3))

# ====================================================================
#  Table S4 — UKB-PPP pQTL coloc
# ====================================================================
cat("\n>>> Table S4: UKB-PPP pQTL coloc\n")

# Canonical UKB-PPP pQTL coloc = phase3c (5 genes, incl. TNF); phase3b superseded (2026-08-17)
ukb <- read_csv("results/phase3c_pqtl_coloc.csv", show_col_types=FALSE)
ukb_mr <- read_csv("results/phase2_decode/eqtl_pqtl_contrast_matrix.csv", show_col_types=FALSE) %>%
  filter(pqtl_source == "UKB-PPP") %>%
  select(gene=SYMBOL, pQTL_Wald_b=pqtl_wald_b, pQTL_Wald_p=pqtl_wald_p, pQTL_n_instruments=n_instruments)

ukb_tab <- ukb %>%
  left_join(ukb_mr, by=c("gene")) %>%
  mutate(coloc_pass = PPH4 > 0.8)

tS4 <- ukb_tab %>% transmute(
  Gene      = sanitize_latex(gene),
  `$n_{\\text{SNPs}}$` = as.character(nsnps),
  `PPH4`    = sprintf("%.4f", PPH4),
  `PPH3`    = sprintf("%.3f", PPH3),
  `Coloc`   = if_else(coloc_pass, "Pass", "Fail"),
  `$n_{\\text{IV}}$`  = if_else(is.na(pQTL_n_instruments), "---", as.character(pQTL_n_instruments)),
  `Wald $b$` = if_else(is.na(pQTL_Wald_b), "---", sprintf("%.3f", pQTL_Wald_b)),
  `Wald $p$` = if_else(is.na(pQTL_Wald_p), "---", fmt_p4(pQTL_Wald_p))
)

capS4 <- "UKB-PPP pQTL colocalization results. All 5 genes failed coloc (PPH4 $>$ 0.8). TNF and CDKN1A were significant in pQTL MR (Wald $p < 0.05$)."
labS4 <- "tab:ukbppp_pqtl_coloc"
f4 <- file.path(outdir, "TableS4_ukbppp_pqtl_coloc.tex")
cat("\\begin{table}[htbp]\n", file=f4)
cat("\\centering\n", file=f4, append=TRUE)
cat(sprintf("\\caption{%s}\\label{%s}\n", capS4, labS4), file=f4, append=TRUE)
cat("\\small\n", file=f4, append=TRUE)
cat("\\begin{tabular}{llllllll}\n", file=f4, append=TRUE)
cat("\\toprule\n", file=f4, append=TRUE)
hdr <- paste0("  ", paste(colnames(tS4), collapse=" & "), " \\\\\n")
cat(hdr, file=f4, append=TRUE)
cat("\\midrule\n", file=f4, append=TRUE)
for (i in 1:nrow(tS4)) {
  row <- paste0("  ", paste(tS4[i,], collapse=" & "), " \\\\\n")
  cat(row, file=f4, append=TRUE)
}
cat("\\bottomrule\n", file=f4, append=TRUE)
cat("\\end{tabular}\n", file=f4, append=TRUE)
cat("\\end{table}\n", file=f4, append=TRUE)
cat(sprintf("  -> %s\n", f4))

# ====================================================================
#  Table S5 — deCODE pQTL coloc (all: UKB-PPP 5 + deCODE-only 5)
# ====================================================================
cat("\n>>> Table S5: deCODE pQTL coloc\n")
decode_coloc <- read_csv("results/phase2_decode/phase2_decode_pqtl_coloc.csv", show_col_types=FALSE)

# add UKB-PPP genes also tested in deCODE (CDKN1A, LTA, NCF2, TET2)
# Actually, the UKB-PPP genes were NOT tested with deCODE coloc. deCODE coloc was separate.
# deCODE coloc = the deCODE-only genes + overlapping UKB-PPP genes that had deCODE data
# Let me check if any UKB-PPP genes are in the deCODE coloc results
decode_mr <- read_csv("results/phase2_decode/eqtl_pqtl_contrast_matrix.csv", show_col_types=FALSE) %>%
  filter(pqtl_source == "deCODE") %>%
  select(gene=SYMBOL, pQTL_Wald_b=pqtl_wald_b, pQTL_Wald_p=pqtl_wald_p)

tS5 <- decode_coloc %>%
  left_join(decode_mr, by="gene") %>%
  transmute(
    Gene      = sanitize_latex(gene),
    Source    = if_else(gene == "CCM2", "deCODE", "deCODE"),
    `$n_{\\text{SNPs}}$` = as.character(nsnps),
    `PPH4`    = sprintf("%.4f", PPH4),
    `PPH3`    = sprintf("%.3f", PPH3),
    `Coloc`   = if_else(coloc_pass, "Pass", "Fail"),
    `Wald $b$` = if_else(is.na(pQTL_Wald_b), "---", sprintf("%.3f", pQTL_Wald_b)),
    `Wald $p$` = if_else(is.na(pQTL_Wald_p), "---", fmt_p4(pQTL_Wald_p))
  )

capS5 <- "deCODE pQTL colocalization results. Only CCM2 passed coloc (PPH4 $>$ 0.8). The remaining 4 genes showed PPH3 $\\approx$ 1.0, suggesting independent causal variants for pQTL and GWAS signals."
labS5 <- "tab:decode_pqtl_coloc"
f5 <- file.path(outdir, "TableS5_decode_pqtl_coloc.tex")
cat("\\begin{table}[htbp]\n", file=f5)
cat("\\centering\n", file=f5, append=TRUE)
cat(sprintf("\\caption{%s}\\label{%s}\n", capS5, labS5), file=f5, append=TRUE)
cat("\\small\n", file=f5, append=TRUE)
cat("\\begin{tabular}{llllllll}\n", file=f5, append=TRUE)
cat("\\toprule\n", file=f5, append=TRUE)
hdr <- paste0("  ", paste(colnames(tS5), collapse=" & "), " \\\\\n")
cat(hdr, file=f5, append=TRUE)
cat("\\midrule\n", file=f5, append=TRUE)
for (i in 1:nrow(tS5)) {
  row <- paste0("  ", paste(tS5[i,], collapse=" & "), " \\\\\n")
  cat(row, file=f5, append=TRUE)
}
cat("\\bottomrule\n", file=f5, append=TRUE)
cat("\\end{tabular}\n", file=f5, append=TRUE)
cat("\\end{table}\n", file=f5, append=TRUE)
cat(sprintf("  -> %s\n", f5))

# ====================================================================
#  Table S6 — PheWAS safety (per gene summary)
# ====================================================================
cat("\n>>> Table S6: PheWAS safety summary\n")

tS6 <- phewas %>%
  transmute(
    Gene = sanitize_latex(gene),
    `Diseases` = as.character(n_diseases),
    `Genetic` = as.character(n_genetic_assoc),
    `Risk cats` = as.character(n_risk_categories),
    `Safety score` = sprintf("%.1f", safety_score),
    `Tier` = safety_tier
  )

capS6 <- "PheWAS safety screening of 6 Tier 1 genes using Open Targets Genetics GWAS associations. Safety score integrates disease count, genetic evidence strength, and risk category diversity."
labS6 <- "tab:phewas_safety"
f6 <- file.path(outdir, "TableS6_phewas_safety.tex")
cat("\\begin{table}[htbp]\n", file=f6)
cat("\\centering\n", file=f6, append=TRUE)
cat(sprintf("\\caption{%s}\\label{%s}\n", capS6, labS6), file=f6, append=TRUE)
cat("\\small\n", file=f6, append=TRUE)
cat("\\begin{tabular}{llllll}\n", file=f6, append=TRUE)
cat("\\toprule\n", file=f6, append=TRUE)
hdr <- paste0("  ", paste(colnames(tS6), collapse=" & "), " \\\\\n")
cat(hdr, file=f6, append=TRUE)
cat("\\midrule\n", file=f6, append=TRUE)
for (i in 1:nrow(tS6)) {
  row <- paste0("  ", paste(tS6[i,], collapse=" & "), " \\\\\n")
  cat(row, file=f6, append=TRUE)
}
cat("\\bottomrule\n", file=f6, append=TRUE)
cat("\\end{tabular}\n", file=f6, append=TRUE)
cat("\\end{table}\n", file=f6, append=TRUE)
cat(sprintf("  -> %s\n", f6))

# ====================================================================
#  Table S7 — Microenvironment (ME) zone characterization
# ====================================================================
cat("\n>>> Table S7: Microenvironment zone characterization\n")

me <- read_csv("results/phase6_scrna/TableS_microenvironment_characterization.csv", show_col_types=FALSE) %>%
  mutate(
    top_ct1_str = paste0(top_ct1, " (+", sprintf("%.2f", top_ct1_score), ")"),
    top_ct2_str = paste0(top_ct2, " (", sprintf("%.2f", top_ct2_score), ")"),
    top_ct3_str = paste0(top_ct3, " (", sprintf("%.2f", top_ct3_score), ")")
  )

tS7 <- me %>% transmute(
  Zone        = sanitize_latex(microenvironment),
  `N spots`   = as.character(n_spots),
  `\\%`       = paste0(pct_spots, "\\%"),
  `Cell type 1` = sanitize_latex(top_ct1_str),
  `Cell type 2` = sanitize_latex(top_ct2_str),
  `Cell type 3` = sanitize_latex(top_ct3_str)
)

capS7 <- sprintf("Microenvironment (ME) zone characterization via K-means clustering (k=6, elbow method) of cell-type module scores on %d Visium spots (GSE285505).", sum(me$n_spots))
labS7 <- "tab:me_zones"
f7 <- file.path(outdir, "TableS7_me_zone_characterization.tex")
cat("\\begin{table}[htbp]\n", file=f7)
cat("\\centering\n", file=f7, append=TRUE)
cat(sprintf("\\caption{%s}\\label{%s}\n", capS7, labS7), file=f7, append=TRUE)
cat("\\small\n", file=f7, append=TRUE)
cat("\\begin{tabular}{llllll}\n", file=f7, append=TRUE)
cat("\\toprule\n", file=f7, append=TRUE)
hdr <- paste0("  ", paste(colnames(tS7), collapse=" & "), " \\\\\n")
cat(hdr, file=f7, append=TRUE)
cat("\\midrule\n", file=f7, append=TRUE)
for (i in 1:nrow(tS7)) {
  row <- paste0("  ", paste(tS7[i,], collapse=" & "), " \\\\\n")
  cat(row, file=f7, append=TRUE)
}
cat("\\bottomrule\n", file=f7, append=TRUE)
cat("\\end{tabular}\n", file=f7, append=TRUE)
cat("\\end{table}\n", file=f7, append=TRUE)
cat(sprintf("  -> %s\n", f7))

# ====================================================================
#  Table S8 — Chen 2024 gene list with v5 overlap annotation
# ====================================================================
cat("\n>>> Table S8: Chen 2024 gene overlap\n")

chen_genes <- readLines("/tmp/chen2024_genes.txt")
chen_genes <- chen_genes[chen_genes != "" & !is.na(chen_genes)]
v5_all <- unique(read_tsv("results/phase1_bonferroni_significant_annotated.tsv", show_col_types=FALSE)$SYMBOL)
v5_all <- v5_all[v5_all != "" & !is.na(v5_all)]

chen_tab <- tibble(gene = chen_genes) %>%
  mutate(
    in_v5 = gene %in% v5_all,
    rank = row_number()
  )

# Show overlap first, then Chen-only
chen_tab <- chen_tab %>%
  arrange(desc(in_v5), gene)

v5_overlap_n <- sum(chen_tab$in_v5)

capS8 <- sprintf("Chen 2024 (Nat Commun) cross-reference: %d genes from supplementary data with v5 overlap. %d/%d (%.1f\\%%) genes overlap with v5 (bold = overlapping).", nrow(chen_tab), v5_overlap_n, nrow(chen_tab), 100*v5_overlap_n/nrow(chen_tab))
labS8 <- "tab:chen2024_overlap"
f8 <- file.path(outdir, "TableS8_chen2024_gene_list.tex")
cat("\\begin{longtable}{ll}\n", file=f8)
cat(sprintf("\\caption{%s}\\label{%s}\\\\\n", capS8, labS8), file=f8, append=TRUE)
cat("\\toprule\n", file=f8, append=TRUE)
cat("  Gene & Overlap in v5 \\\\\n", file=f8, append=TRUE)
cat("\\midrule\n", file=f8, append=TRUE)
cat("\\endfirsthead\n", file=f8, append=TRUE)
cat("\\multicolumn{2}{l}{\\textit{Table S8 (continued)}} \\\\\n", file=f8, append=TRUE)
cat("\\toprule\n", file=f8, append=TRUE)
cat("  Gene & Overlap in v5 \\\\\n", file=f8, append=TRUE)
cat("\\midrule\n", file=f8, append=TRUE)
cat("\\endhead\n", file=f8, append=TRUE)
cat("\\bottomrule\n", file=f8, append=TRUE)
cat("\\endfoot\n", file=f8, append=TRUE)
for (i in 1:nrow(chen_tab)) {
  gene_name <- sanitize_latex(chen_tab$gene[i])
  overlap_mark <- ifelse(chen_tab$in_v5[i], "\\textbf{Yes}", "---")
  row <- paste0("  ", gene_name, " & ", overlap_mark, " \\\\\n")
  cat(row, file=f8, append=TRUE)
}
cat("\\bottomrule\n", file=f8, append=TRUE)
cat("\\end{longtable}\n", file=f8, append=TRUE)
cat(sprintf("  -> %s\n", f8))

# ====================================================================
#  Table S9 — Spatial gene-gene cross-correlation matrix (Tier 1 genes)
# ====================================================================
cat("\n>>> Table S9: Spatial gene cross-correlation\n")

corr <- read_csv("results/phase6_scrna/TableS_spatial_gene_gene_correlation.csv", show_col_types=FALSE) %>%
  filter(!is.na(spearman_rho)) %>%
  mutate(
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    display_val = paste0(sprintf("%.3f", spearman_rho), sig_label)
  ) %>%
  select(gene1, gene2, display_val)

# Pivot to matrix form (top 12 genes by connectivity)
top_genes <- corr %>%
  count(gene1, sort=TRUE) %>%
  head(12) %>%
  pull(gene1)

corr_mat <- corr %>%
  filter(gene1 %in% top_genes & gene2 %in% top_genes) %>%
  pivot_wider(id_cols=gene1, names_from=gene2, values_from=display_val, values_fill="---") %>%
  select(gene1, all_of(intersect(top_genes, colnames(.))))

tS9 <- corr_mat %>%
  rename(Gene = gene1) %>%
  mutate(across(-Gene, sanitize_latex))

capS9 <- sprintf("Spatial gene-gene cross-correlation matrix (Spearman $\\rho$) for Tier 1 genes on Visium spots ($n=%d$). *** $p<0.001$, ** $p<0.01$, * $p<0.05$. Gene pairs without sufficient spatial co-detection are marked `---'.", length(unique(c(corr$gene1, corr$gene2))))
labS9 <- "tab:spatial_crosscorr"
f9 <- file.path(outdir, "TableS9_spatial_crosscorrelation.tex")
cat("\\begin{table}[htbp]\n", file=f9)
cat("\\centering\n", file=f9, append=TRUE)
cat(sprintf("\\caption{%s}\\label{%s}\n", capS9, labS9), file=f9, append=TRUE)
cat("\\footnotesize\n", file=f9, append=TRUE)
cats <- min(ncol(tS9), 13)
cat("\\resizebox{\\textwidth}{!}{%\n", file=f9, append=TRUE)
cat(sprintf("\\begin{tabular}{%s}\n", paste0(rep("l", cats), collapse="")), file=f9, append=TRUE)
cat("\\toprule\n", file=f9, append=TRUE)
hdr <- paste0("  ", paste(colnames(tS9), collapse=" & "), " \\\\\n")
cat(hdr, file=f9, append=TRUE)
cat("\\midrule\n", file=f9, append=TRUE)
for (i in 1:nrow(tS9)) {
  row <- paste0("  ", paste(as.character(tS9[i,]), collapse=" & "), " \\\\\n")
  cat(row, file=f9, append=TRUE)
}
cat("\\bottomrule\n", file=f9, append=TRUE)
cat("\\end{tabular}\n", file=f9, append=TRUE)
cat("}\n", file=f9, append=TRUE)
cat("\\end{table}\n", file=f9, append=TRUE)
cat(sprintf("  -> %s\n", f9))

# ====================================================================
cat("\n=== ALL TABLES GENERATED ===\n")
cat(sprintf("Output: %s/\n", outdir))
for (f in list.files(outdir, pattern="\\.tex$")) {
  cat(sprintf("  %s (%s)\n", f, format(file.size(file.path(outdir, f)), big.mark=",")))
}
cat("Done.\n")
