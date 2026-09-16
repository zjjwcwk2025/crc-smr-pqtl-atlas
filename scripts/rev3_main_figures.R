#!/usr/bin/env Rscript
# rev3 main figures - regenerated at final print size (6.30 in = textwidth)
suppressMessages({
  library(data.table); library(ggplot2); library(dplyr); library(readr)
  library(tidyr); library(tibble); library(forcats); library(stringr)
  library(ggrepel); library(scales); library(cowplot); library(ComplexUpset)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "results/figures_rev3"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
TEXTW <- 6.30

source("scripts/theme_pub.R")
theme_nc <- theme_pub()
palette_coloc <- c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#CA0020")
palette_class <- c("A: pQTL concordant" = "#2166AC", "B: pQTL discordant" = "#B2182B", "C: pQTL not significant" = "#999999")
palette_direction <- c("concordant" = "#2166AC", "discordant" = "#B2182B", "not_significant" = "#999999")
palette_safety <- c("High alert" = "#B2182B", "Moderate" = "#F4A582", "Low" = "#92C5DE")
# ==============================================================
#  Fig 2a — coloc PPH4 barplot (phase3_coloc_all75.csv)
# ==============================================================
cat("\n>>> Figure 1c: coloc PPH4 barplot (all 76 probe-gene pairs)\n")
coloc <- read_csv("results/phase3_coloc_all75.csv", show_col_types=FALSE)
# Handle column name variants
if("SYMBOL" %in% names(coloc)) {
  coloc <- coloc %>% rename(gene = SYMBOL)
}
if("pph4" %in% names(coloc)) {
  coloc <- coloc %>% rename(PPH4 = pph4)
}
# Map bare ENSG IDs to SYMBOL + biotype (铁律 #8: no silent bare-ENSG labels)
biotype_map <- read_tsv("results/phase1_ensg_biotype.tsv", show_col_types=FALSE) %>%
  select(ensg = ENSG, gene_type, gene_name)
coloc <- coloc %>%
  left_join(biotype_map, by = "ensg") %>%
  mutate(
    gene_label = case_when(
      !grepl("^ENSG", gene) ~ gene,
      !is.na(gene_name) & gene_name != "" & !grepl("^ENSG", gene_name) ~ gene_name,
      !is.na(gene_type) & gene_type == "not_in_GENCODE_v47" ~ paste0(ensg, " (not in GENCODE)"),
      !is.na(gene_type) ~ paste0(ensg, " (", gene_type, ")"),
      TRUE ~ ensg
    ),
    coloc_status = case_when(
      PPH4 > 0.8 ~ "PASS (PPH4 > 0.8)",
      PPH4 >= 0.5 ~ "Marginal (0.5-0.8)",
      PPH4 >= 0.01 ~ "Weak (0.01-0.5)",
      TRUE ~ "FAIL (PPH4 < 0.01)"
    ),
    is_mhc = grepl("HLA-|HCG27|TNF$|LTA$|ZFP57", gene)
  )

# ── Fig 2a: coloc PPH4 dot/lollipop plot (curated subset, CNS-compliant) ──
# Show the colocalization-supported subset (PASS + MARGINAL, PPH4 > 0.5) on a
# log10 PPH4 axis as a dot/lollipop plot, and collapse the low-PPH4 tail
# (WEAK + FAIL, PPH4 <= 0.5) into a single data-driven text statement.
# A log-scale geom_col is a category error (col maps 0 -> -Inf), which clamps
# the 21 FAIL genes to the panel floor and renders them as an unreadable wall.

# Data-driven counts (铁律 #14: no hand-typed 43/76/49/27)
n_total  <- nrow(coloc)
n_pass   <- sum(coloc$status == "PASS")
n_margin <- sum(coloc$status == "MARGINAL")
n_weak   <- sum(coloc$status == "WEAK")
n_fail   <- sum(coloc$status == "FAIL")
n_hidden <- n_weak + n_fail
n_show   <- n_pass + n_margin

status_lvl <- c(
  sprintf("PASS (PPH4 > 0.8; n = %d)", n_pass),
  sprintf("Marginal (0.5\u20130.8; n = %d)", n_margin)
)
status_vals <- setNames(c("#2166AC", "#92C5DE"), status_lvl)   # colourblind-safe blue pair

# Curated subset = pairs with coloc-supporting PPH4 (PASS + MARGINAL).
coloc_show <- coloc %>%
  filter(status %in% c("PASS", "MARGINAL")) %>%
  mutate(
    coloc_status = factor(status, levels = c("PASS", "MARGINAL"), labels = status_lvl),
    PPH4_plot    = PPH4,
    gene_label   = fct_reorder(gene_label, PPH4_plot)   # highest PPH4 at top
  )

p2a <- ggplot(coloc_show, aes(x = PPH4_plot, y = gene_label)) +
  # lollipop stem from a baseline just below the marginal threshold
  geom_segment(aes(x = 0.3, xend = PPH4_plot, y = gene_label, yend = gene_label),
               colour = "#BDBDBD", linewidth = 0.5) +
  geom_point(aes(colour = coloc_status), size = 2.5) +
  scale_colour_manual(name = "Coloc status", values = status_vals, drop = FALSE) +
  # reference thresholds
  geom_vline(xintercept = 0.8, linetype = "dashed", colour = "#B2182B", linewidth = 0.6) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "#6B6B6B", linewidth = 0.5) +
  scale_x_log10(
    breaks = c(0.3, 0.5, 0.8, 1.0),
    labels = c("0.3", "0.5", "0.8", "1.0"),
    limits = c(0.28, 1.03),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Posterior Probability of Colocalization (PPH4, log10)",
    y = NULL,
    title = "eQTL\u2013GWAS colocalization (coloc)",
    subtitle = sprintf("%d/%d pairs (%.0f%%) show coloc support (%d PASS + %d marginal); %d pairs omitted",
                       n_show, n_total, 100 * n_show / n_total, n_pass, n_margin, n_hidden),
    caption = paste0("coloc.abf, 1000G EUR MAF. Red dashed, PPH4 = 0.8 (PASS); grey dotted, PPH4 = 0.5 (marginal). Omitted: ",
                     n_hidden, " pairs (PPH4 \u2264 0.5).")
  ) +
  theme_nc +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.caption = element_blank()
  )

ggsave(file.path(outdir, "Fig1c_coloc_pph4.pdf"), p2a, width = 6.30, height = 7.0, device = cairo_pdf)
cat("  -> Fig1c_coloc_pph4.pdf\n")
# ==============================================================
#  Fig 2b — eQTL/pQTL discordance classification bar
# ==============================================================
cat("\n>>> Figure 2b: eQTL/pQTL classification bar\n")
contrast <- read_csv("results/phase2_decode/eqtl_pqtl_contrast_matrix.csv", show_col_types=FALSE)
# Remove TNFRSF1A, TNFRSF1B (no eQTL SMR data in this table)
contrast <- contrast %>% filter(!is.na(SYMBOL) & SYMBOL != "")
class_counts <- contrast %>%
  count(classification) %>%
  mutate(pct = n/sum(n))

# Simplify labels for plotting
class_counts <- class_counts %>%
  mutate(class_label = case_when(
    grepl("A:", classification) ~ "A: pQTL concordant",
    grepl("B:", classification) ~ "B: pQTL discordant",
    grepl("C:", classification) ~ "C: pQTL not significant",
    TRUE ~ classification
  ))

p2b <- ggplot(class_counts, aes(x=reorder(class_label, n), y=n, fill=class_label)) +
  geom_col(width=0.6, color="black", linewidth=0.3) +
  geom_text(aes(label=paste0(n," (",round(pct*100,0),"%)")), hjust=-0.2, size=4, fontface="bold") +
  scale_fill_manual(values=palette_class, guide="none") +
  scale_y_continuous(limits=c(0, max(class_counts$n)*1.25), expand=c(0,0)) +
  labs(x="", y="Number of genes",
       title="eQTL vs pQTL Direction Concordance",
       subtitle=paste0("n=", sum(class_counts$n), " genes with both eQTL and pQTL results"),
       caption="Classification: A/B/C by pQTL Wald p < 0.05 and direction relative to eQTL SMR b.") +
  coord_flip() +
  theme_nc

ggsave(file.path(outdir, "Fig2b_eqtl_pqtl_classification.pdf"), p2b, width = 6.30, height = 2.9, device = cairo_pdf)
cat("  -> Fig2b_eqtl_pqtl_classification.pdf\n")
# ==============================================================
#  Fig 2c — eQTL/pQTL effect size scatter
# ==============================================================
cat("\n>>> Figure 2c: eQTL vs pQTL effect size scatter\n")
# Use genes that have both eQTL and pQTL betas
contrast_both <- contrast %>%
  filter(!is.na(b_SMR) & !is.na(pqtl_wald_b) & !is.na(pqtl_wald_p))

# Calculate Spearman correlation (铁律 6: 相关用 Spearman, 报 rho + exact p)
r_test <- if(nrow(contrast_both) >= 3) {
  tryCatch(cor.test(contrast_both$b_SMR, contrast_both$pqtl_wald_b,
                    method = "spearman", use = "complete.obs"),
           error = function(e) NULL)
} else NULL
r_val <- if(!is.null(r_test)) r_test$estimate else NA
r_p   <- if(!is.null(r_test)) r_test$p.value else NA
stats_label <- if(!is.na(r_val)) {
  paste0("Spearman rho = ", formatC(r_val, format = "f", digits = 3),
         ", p = ", formatC(r_p, format = "f", digits = 2), " (two-sided)")
} else {
  "Spearman correlation not estimable."
}

# Classify direction
contrast_both <- contrast_both %>%
  mutate(
    dir_status = case_when(
      direction_match == "TRUE" ~ "concordant",
      direction_match == "FALSE" ~ "discordant",
      TRUE ~ "not_significant"
    ),
    point_size = ifelse(pqtl_wald_p < 0.05, 3.5, 2)
  )

p2c <- ggplot(contrast_both, aes(x=b_SMR, y=pqtl_wald_b)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  geom_vline(xintercept=0, linetype="dashed", color="grey60") +
  geom_point(aes(fill=dir_status, size=point_size), shape=21, alpha=0.85) +
  geom_text_repel(aes(label=SYMBOL), size=2.9, max.overlaps=Inf, force=8, box.padding=0.7,
                  point.padding=0.35, min.segment.length=0, seed=5,
                  segment.size=0.25, segment.colour="grey40") +
  scale_fill_manual(values=palette_direction, name="Direction") +
  scale_size_identity() +
  labs(x="eQTL SMR b (blood)", y="pQTL Wald b",
       title="eQTL vs pQTL Effect Size Concordance",
       subtitle=paste0("n=", nrow(contrast_both), " genes with both eQTL and pQTL effect estimates; filled circles = pQTL p<0.05"),
       caption=stats_label) +
  theme_nc

ggsave(file.path(outdir, "Fig2c_eqtl_pqtl_scatter.pdf"), p2c, width = 6.30, height = 5.4, device = cairo_pdf)
cat("  -> Fig2c_eqtl_pqtl_scatter.pdf\n")
# ==============================================================
#  Fig 3 — GTEx Colon vs Blood SMR scatter
# ==============================================================
cat("\n>>> Figure 3b/3c: GTEx Colon SMR sensitivity\n")
gtex <- read_tsv("results/gtex_colon_transverse_smr.smr", show_col_types=FALSE)
# Read blood SMR for matched genes (Gene column = Ensembl ID here)
e1 <- read_tsv("results/phase1_eqtl_smr.smr", show_col_types=FALSE)

# GTEx has probeID=ENSG (same as eQTLGen), Gene=symbol
# eQTLGen has probeID=ENSG, Gene=ENSG also
# Merge on probeID (Ensembl gene ID)
e1_clean <- e1 %>%
  select(probeID, b_SMR_blood=b_SMR, se_SMR_blood=se_SMR, p_SMR_blood=p_SMR)

# Get gene symbols from eQTLGen if GTEx doesn't have them
merged <- gtex %>%
  inner_join(e1_clean, by="probeID") %>%
  mutate(
    sig_both = p_SMR < 3.21e-6 & p_SMR_blood < 3.21e-6,
    sig_status = case_when(
      sig_both ~ "Sig in both",
      p_SMR < 3.21e-6 ~ "Sig in GTEx only",
      p_SMR_blood < 3.21e-6 ~ "Sig in eQTLGen only",
      TRUE ~ "Neither sig"
    )
  )
cat(sprintf("  Merged: %d genes with both GTEx and eQTLGen SMR results\n", nrow(merged)))

ct <- cor.test(merged$b_SMR, merged$b_SMR_blood, method = "spearman", use="complete.obs")
r_rho   <- ct$estimate
r_pval  <- ct$p.value
concordant_n <- sum(sign(merged$b_SMR) == sign(merged$b_SMR_blood))
concordant_pct <- 100 * concordant_n / nrow(merged)

p3 <- ggplot(merged, aes(x=b_SMR_blood, y=b_SMR)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  geom_vline(xintercept=0, linetype="dashed", color="grey60") +
  geom_smooth(method="lm", se=TRUE, color="grey40", linewidth=0.8, alpha=0.15) +
  geom_point(aes(fill=sig_status), shape=21, size=3.5, alpha=0.85) +
  geom_text_repel(aes(label=Gene), size=2.7, max.overlaps=Inf, force=8, box.padding=0.7,
                  point.padding=0.35, min.segment.length=0, seed=7,
                  segment.size=0.25, segment.colour="grey40") +
  scale_x_continuous(expand = expansion(mult = 0.10)) +
  scale_fill_manual(values=c("Sig in both"="#2166AC", "Sig in GTEx only"="#F4A582",
                               "Sig in eQTLGen only"="#92C5DE", "Neither sig"="grey80"),
                    name="Significance") +
  labs(x="SMR b (eQTLGen blood)", y="SMR b (GTEx Colon Transverse)",
       title="Tissue Sensitivity: GTEx Colon vs eQTLGen Blood",
       subtitle=paste0("Spearman rho = ",round(r_rho,2),", p = ",formatC(r_pval,format="f",digits=4),
                        "; direction concordant: ",concordant_n,"/",nrow(merged)," (",round(concordant_pct,1),"%)"),
       caption=paste0("n = ", nrow(merged), " genes with SMR estimates in both tissues; Spearman correlation, two-sided.")) +
  theme_nc

ggsave(file.path(outdir, "Fig3b_gtex_colon_scatter.pdf"), p3, width = 6.30, height = 5.0, device = cairo_pdf)
cat("  -> Fig3b_gtex_colon_scatter.pdf\n")
# ==============================================================
#  Fig 3 — forest plot (companion panel)
# ==============================================================
merged_sorted <- merged %>%
  filter(sig_both) %>%
  arrange(b_SMR) %>%
  mutate(Gene = fct_reorder(Gene, b_SMR))

if(nrow(merged_sorted) > 0) {
  forest_data <- bind_rows(
    merged_sorted %>% mutate(source="GTEx Colon", b=b_SMR, se=se_SMR, p=p_SMR,
                             ci_low=b_SMR-1.96*se_SMR, ci_high=b_SMR+1.96*se_SMR),
    merged_sorted %>% mutate(source="eQTLGen Blood", b=b_SMR_blood, se=se_SMR_blood, p=p_SMR_blood,
                             ci_low=b_SMR_blood-1.96*se_SMR_blood, ci_high=b_SMR_blood+1.96*se_SMR_blood)
  ) %>% mutate(Gene = factor(Gene, levels=levels(merged_sorted$Gene)))

  # Prepare data with CI for both tissues
  forest_ci <- forest_data %>% filter(!is.na(se))

  p3f <- ggplot(forest_data, aes(x=b, y=Gene, color=source)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey60") +
    geom_errorbar(data=forest_ci,
                   aes(xmin=ci_low, xmax=ci_high), width=0.12, linewidth=0.6) +
    geom_point(aes(shape=source), size=3) +
    scale_shape_manual(values=c("GTEx Colon"=16, "eQTLGen Blood"=17)) +
    scale_color_manual(values=c("GTEx Colon"="#CA0020", "eQTLGen Blood"="#2166AC")) +
    labs(x="SMR effect size (b)", y="",
         title="Dual-significant Genes: GTEx Colon vs eQTLGen",
         subtitle=paste0(nrow(merged_sorted)," genes significant in both GTEx Colon and eQTLGen; error bars = 95% CI"),
         caption=paste0("n = ", nrow(merged_sorted), " dual-significant genes; point estimates and 95% CI.")) +
    theme_nc + theme(legend.position=c(0.8,0.15),
                     panel.grid.major.x=element_line(color="grey90", linewidth=0.3),
                     panel.grid.major.y=element_line(color="grey90", linewidth=0.3))

  ggsave(file.path(outdir, "Fig3c_gtex_colon_forest.pdf"), p3f, width = 6.30, height = 4.2, device = cairo_pdf)
  cat("  -> Fig3c_gtex_colon_forest.pdf\n")
}
# ==============================================================
#  Fig 4a — Drug annotation dot plot  (phase5b_deep_drug_annotation.csv)
# ==============================================================
cat("\n>>> Figure 4b: Drug annotation summary\n")
drug <- read_csv("results/phase5b_deep_drug_annotation.csv", show_col_types=FALSE) %>%
  filter(!is.na(gene)) %>%
  mutate(
    has_drug = n_drugs > 0,
    has_approved = n_approved > 0,
    tier_label = case_when(
      has_approved ~ "Approved drug",
      has_drug ~ "Investigational",
      TRUE ~ "No drugs"
    )
  )

# Show top 25 by repurposing_score
drug_top <- drug %>%
  arrange(desc(repurposing_score)) %>%
  slice_head(n=25) %>%
  mutate(gene = fct_reorder(gene, repurposing_score))

p5a <- ggplot(drug_top, aes(x=gene, y=repurposing_score)) +
  geom_col(aes(fill=tier_label), width=0.7, color="black", linewidth=0.2) +
  geom_text(aes(label=ifelse(n_drugs>0, paste0(n_drugs," drug",ifelse(n_drugs>1,"s","")), ""),
                y=repurposing_score+0.03), size=2.8, fontface="italic") +
  scale_fill_manual(values=c("Approved drug"="#2166AC", "Investigational"="#F4A582", "No drugs"="grey85"),
                    name="Drug status") +
  scale_y_continuous(limits=c(0, max(drug_top$repurposing_score)*1.15), expand=c(0,0)) +
  labs(x="", y="Drug Repurposing Score",
       title="Drug Repurposing Potential",
       subtitle=paste0("Top ", nrow(drug_top), " of ", nrow(drug), " genes by integrated repurposing score (SMR+coloc+SuSiE+TCGA+scRNA+drugs)"),
       caption="Score components and weights are defined in Methods; n = 25 genes shown.") +
  coord_flip() +
  theme_nc

ggsave(file.path(outdir, "Fig4b_drug_repurposing.pdf"), p5a, width = 6.30, height = 4.9, device = cairo_pdf)
cat("  -> Fig4b_drug_repurposing.pdf\n")
# ==============================================================
#  Fig 4b — PheWAS safety heatmap  (phase5c_phewas_safety.csv)
# ==============================================================
cat("\n>>> Figure 4a: PheWAS safety heatmap\n")
phewas <- read_csv("results/phase5c_phewas_safety.csv", show_col_types=FALSE)

# Build risk category matrix
risk_cat_long <- phewas %>%
  separate_rows(category_breakdown, sep="\\s*\\|\\s*") %>%
  mutate(
    category_name = str_trim(str_extract(category_breakdown, "^[^=]+")),
    category_count = as.numeric(str_extract(category_breakdown, "\\d+$"))
  ) %>%
  select(gene, category_name, category_count) %>%
  filter(!is.na(category_count)) %>%
  distinct(gene, category_name, .keep_all=TRUE)

p5b <- ggplot(risk_cat_long, aes(x=category_name, y=gene, fill=category_count)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=category_count), size=3.5, fontface="bold") +
  scale_fill_gradient(low="#F7FBFF", high="#08306B", name="N diseases") +
  labs(x="Disease Category", y="",
       title="PheWAS Disease Association Profile",
       subtitle="Tier 1 genes: Open Targets Platform disease associations by category",
       caption="n = 6 Tier 1 genes; tile value = number of associated diseases per category.") +
  theme_nc + theme(
    axis.text.x = element_text(angle=45, hjust=1, size=8.5),
    plot.margin = margin(5, 6, 3, 30),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(outdir, "Fig4a_phewas_safety.pdf"), p5b, width = 6.30, height = 3.2, device = cairo_pdf)
cat("  -> Fig4a_phewas_safety.pdf\n")
# ==============================================================
#  Fig 4c — Competitor overlap  (competitor_overlap.tsv + UpSet)
# ==============================================================
cat("\n>>> Figure 5a: Competitor overlap\n")
comp <- read_csv("results/competitor_overlap.tsv", show_col_types=FALSE)

# Chen 2024 is absent from competitor_overlap.tsv. Add it from Table S8 (the
# authoritative cross-reference) so the bar matches the 44.9% quoted in the text.
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
cat(sprintf("  Chen 2024 added: %d/%d overlap\n", n_chen_ov, length(v5_sym)))

# Barplot of overlap rates
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

p5c_bar <- ggplot(comp_plot, aes(y=study_short, x=overlap_n)) +
  geom_col(aes(fill=journal_if_num), width=0.6, color="black", linewidth=0.3) +
  geom_text(aes(label=paste0(overlap_n, " (", overlap_rate, ")")),
            hjust=-0.12, size=3.5, fontface="bold") +
  scale_fill_gradient(low="#FDDBC7", high="#2166AC", name="Journal IF") +
  scale_x_continuous(limits=c(0, max(comp_plot$overlap_n)*1.35), expand=c(0,0)) +
  labs(x="Overlapping genes with this study", y="",
       title="Overlap with prior CRC target-discovery studies",
       subtitle=paste0("This study: 69 unique genes. ",comp$overlap_n[comp$study=="ALL_COMBINED"],
                       " overlap (", comp$overlap_rate[comp$study=="ALL_COMBINED"],")"),
       caption="Overlap counts are gene-symbol based; journal impact factor shown by fill.") +
  theme_nc +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3))

ggsave(file.path(outdir, "Fig5a_competitor_overlap.pdf"), p5c_bar, width = 6.30, height = 3.5, device = cairo_pdf)
cat("  -> Fig5a_competitor_overlap.pdf\n")

# Fig 5b (UpSet) is generated by scripts/63_fig5b_upset.R, which reads Table S8
# directly. The previous ComplexUpset version parsed only the highlighted rows and
# drew background ribbons that did not line up with the dot matrix.

# ==============================================================
#  Fig 1 - study design (print size, no in-figure title)
# ==============================================================
cat("\n>>> Figure 1b: study design\n")
rows <- c(4.85, 2.95, 1.05)
bands <- data.frame(
  ymin = rows - 0.85, ymax = rows + 0.85,
  fill = c("#EEF4FA", "#EEF4FA", "#F2F2F2"),
  label = c("Discovery", "Replication &\ntissue validation", "Translation")
)
N_CIS_PQTL <- nrow(fread("results/phase2_decode/phase2_decode_pqtl_mr_combined.csv")) +
               nrow(fread("results/phase2_pqtl/phase2_pqtl_mr_combined.csv"))
boxes <- data.frame(
  x = rep(c(4.3, 7.9, 11.5), 3),
  y = rep(rows, each = 3),
  fill = c("#D6E4F2", "#D6E4F2", "#D6E4F2", "#E6EEF7", "#E6EEF7", "#E6EEF7", "#EDEDED", "#EDEDED", "#EDEDED"),
  label = c(
    "SMR: eQTL \u2192 GWAS\n75 genes pass Bonferroni\n(62 protein-coding)",
    "Colocalization + SuSiE\n43 of 76 probe-gene pairs\ncolocalize",
    sprintf("pQTL MR + colocalization\n1 of %d genes with cis-pQTL\ninstruments is convergent", N_CIS_PQTL),
    "FinnGen R13 replication\n1 of 8 Tier 1 genes\nreplicated",
    "GTEx colon eQTL\ntissue-of-action\nsensitivity",
    "scRNA-seq + Visium spatial\n6 microenvironment\nzones",
    "Drug repurposing\n8 of 62 genes with\nknown drugs",
    "IDG druggability\n19 Tier 1 genes\ngraded",
    "PheWAS safety screen\n6 core genes on two\ntransparent dimensions"
  ),
  stringsAsFactors = FALSE
)
boxes$w <- 3.3; boxes$h <- 1.2
boxes$xmin <- boxes$x - boxes$w/2; boxes$xmax <- boxes$x + boxes$w/2
boxes$ymin <- boxes$y - boxes$h/2; boxes$ymax <- boxes$y + boxes$h/2
arrows <- data.frame(x = rep(c(5.95, 9.55), 3),
                     y = rep(rows, each = 2),
                     xend = rep(c(6.25, 9.85), 3),
                     yend = rep(rows, each = 2))
p1 <- ggplot() +
  geom_rect(data = bands, aes(xmin = 0.10, xmax = 13.70, ymin = ymin, ymax = ymax, fill = fill), color = NA) +
  scale_fill_identity() +
  geom_text(data = bands, aes(x = 1.00, y = (ymin + ymax)/2, label = label),
            size = 3.0, fontface = "bold", color = "grey25", lineheight = 1.0) +
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "grey30", linewidth = 0.35) +
  geom_text(data = boxes, aes(x = x, y = y, label = label),
            size = 2.75, color = "grey10", fontface = "bold", lineheight = 1.15) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.10, "cm"), type = "closed"), color = "grey45", linewidth = 0.5) +
  annotate("text", x = 6.90, y = -0.22, hjust = 0.5, size = 2.6, color = "grey35",
           label = "Data: 78,473-case CRC GWAS (GCST90255675); eQTLGen (n = 31,684); deCODE and UKB-PPP plasma proteomes;\nFinnGen R13; TCGA-COAD and READ; Visium spatial transcriptomics.") +
  coord_cartesian(xlim = c(0, 13.8), ylim = c(-0.78, 5.7), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(4, 4, 4, 4))
ggsave(file.path(outdir, "Fig1b_study_design.pdf"), p1, width = TEXTW, height = 3.56, device = cairo_pdf)

# ==============================================================
#  Fig 2a - evidence funnel: 62 -> 15 (measured) -> 7 of 15 / 9 of 17 (cis-pQTL instruments) -> 1
# ==============================================================
cat("\n>>> Figure 2a: evidence funnel\n")
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
funnel <- data.frame(
  y = c(4, 3, 2, 1),
  n = c(n_coding, length(measured), sum(with_inst %in% measured), 1),
  head = c("Protein-coding genes",
           "With plasma protein measurement",
           "With genome-wide significant cis-pQTL instruments",
           "With pQTL-GWAS colocalization"),
  detail = c("of 75 genes passing Bonferroni-corrected SMR (p < 3.21e-6)",
             sprintf("deCODE %d + UKB-PPP %d, %d in both; +2 TNF-pathway genes outside the Bonferroni set (TNFRSF1A, TNFRSF1B) give 17 genes analysed",
                     length(intersect(dec, coding)), length(ukb), length(intersect(dec, ukb))),
             sprintf("%d of these 15; %d of the 17-gene analysis set", n_inst_meas, sum(with_inst %in% all_meas)),
             "CCM2 only (PPH4 = 0.95)"),
  stringsAsFactors = FALSE
)
# Wrap long provenance strings so that nothing runs off the right edge of the panel
wrap_lab <- function(x, width = 66) paste(strwrap(x, width = width), collapse = "\n")
funnel$detail <- vapply(funnel$detail, wrap_lab, character(1))
note_2a <- wrap_lab("Of the six core Tier 1 genes, five have no protein measurement in either panel; BMP2 is measured but carries no cis-pQTL instruments.", width = 100)
funnel$bar <- pmax(0.38 * funnel$n / 62, 0.020)
funnel$fill <- c("#08519C", "#2171B5", "#4393C3", "#B2182B")
p2a <- ggplot(funnel) +
  geom_rect(aes(xmin = 0, xmax = bar, ymin = y - 0.30, ymax = y + 0.30, fill = fill), color = NA) +
  scale_fill_identity() +
  geom_text(aes(x = bar + 0.015, y = y, label = n), hjust = 0, size = 3.6, fontface = "bold", color = "grey20") +
  geom_text(aes(x = 0.46, y = y + 0.14, label = head), hjust = 0, size = 2.9, fontface = "bold", color = "grey20") +
  geom_text(aes(x = 0.46, y = y - 0.245, label = detail), hjust = 0, size = 2.45, color = "grey45", lineheight = 0.95) +
  annotate("text", x = 0, y = 0.22, hjust = 0, size = 2.45, color = "grey35", lineheight = 0.95,
           label = note_2a) +
  coord_cartesian(xlim = c(0, 1.0), ylim = c(0, 4.75), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(4, 4, 4, 4))
ggsave(file.path(outdir, "Fig2a_coverage_funnel.pdf"), p2a, width = TEXTW, height = 2.5, device = cairo_pdf)

# ==============================================================
#  Fig 3a - pQTL convergence: MR effect and colocalization, aligned
# ==============================================================
cat("\n>>> Figure 3a: pQTL convergence panel\n")
mr <- fread("results/phase2_decode/phase2_decode_pqtl_mr_combined.csv")
co <- fread("results/phase2_decode/phase2_decode_pqtl_coloc.csv")
dd <- merge(mr[, .(gene, wald_b, wald_se, wald_pval, n_instruments)], co[, .(gene, PPH4)], by = "gene")
dd[, lo := wald_b - 1.96 * wald_se]
dd[, hi := wald_b + 1.96 * wald_se]
dd[, hl := ifelse(gene == "CCM2", "CCM2", "other")]
dd[, gene := factor(gene, levels = rev(c("CCM2","LIMA1","STAT6","TNFRSF1A","TNFRSF1B")))]
cols <- c("CCM2" = "#B2182B", "other" = "grey55")
pL <- ggplot(dd, aes(x = wald_b, y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.55, colour = "grey45") +
  geom_point(aes(colour = hl), size = 2.2) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = "pQTL effect on CRC (Wald ratio b, 95% CI)", y = NULL) +
  theme_nc + theme(legend.position = "none", plot.margin = margin(4, 6, 2, 16))
pR <- ggplot(dd, aes(x = PPH4, y = gene)) +
  geom_vline(xintercept = 0.8, linetype = "dashed", colour = "#B2182B", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = PPH4, y = gene, yend = gene), colour = "grey85", linewidth = 0.4) +
  geom_point(aes(colour = hl), size = 2.2) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "pQTL-GWAS colocalization (PPH4)", y = NULL) +
  theme_nc + theme(legend.position = "none", axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(), axis.line.y = element_blank())
p3a <- plot_grid(pL, pR, nrow = 1, rel_widths = c(1.25, 1))
ggsave(file.path(outdir, "Fig3a_ccm2_convergence.pdf"), p3a, width = TEXTW, height = 2.4, device = cairo_pdf)


# ==============================================================
#  Final panel file names.
#  Every delivered panel is named after the figure number it carries in
#  the manuscript (Figure 1b -> Fig1b_*.pdf ... Figure S19 -> FigS19_*.pdf),
#  so the file name can never drift away from the printed figure number.
#  The main panels are written directly under their final names above;
#  the supplementary panels are produced by other scripts and are
#  collected into results/figures_rev3/ by the table below.
# ==============================================================
panel_names <- c(
  "../figures/FigS1a_manhattan.pdf"                    = "FigS1a_manhattan.pdf",
  "../figures/FigS1b_qq_plot.pdf"                      = "FigS1b_qq_plot.pdf",
  "../figures/FigS23_tier1_regional_plots.pdf"         = "FigS2_tier1_regional_plots.pdf",
  "../figures/FigS3_susie_sensitivity.pdf"             = "FigS3_susie_sensitivity.pdf",
  "../figures/FigS2_smr_or_forest.pdf"                 = "FigS4_smr_or_forest.pdf",
  "../figures/FigS6_mr_sensitivity.pdf"                = "FigS5_mr_sensitivity.pdf",
  "../figures/FigS4_ukbppp_pqtl_coloc.pdf"             = "FigS6_ukbppp_pqtl_coloc.pdf",
  "../figures/FigS5_decode_pqtl_coloc.pdf"             = "FigS7_decode_pqtl_coloc.pdf",
  "../figures/FigS20_pqtl_mr_decode_forest.pdf"        = "FigS8_pqtl_mr_decode_forest.pdf",
  "../figures/FigS21_pqtl_mr_ukbppp_forest.pdf"        = "FigS9_pqtl_mr_ukbppp_forest.pdf",
  "../figures/FigS24_pqtl_mr_scatter.pdf"              = "FigS10_pqtl_mr_scatter.pdf",
  "../figures/FigS25_pqtl_mr_loo.pdf"                  = "FigS11_pqtl_mr_loo.pdf",
  "../figures/FigS26_pqtl_mr_methods_forest.pdf"       = "FigS12_pqtl_mr_methods_forest.pdf",
  "../figures/FigS27_pqtl_mr_funnel_pleiotropy.pdf"    = "FigS13_pqtl_mr_funnel_pleiotropy.pdf",
  "../figures/FigS7_celltype_dotplot.pdf"              = "FigS14_celltype_dotplot.pdf",
  "../figures/FigS8_spatial_tier1_genes.pdf"           = "FigS15_spatial_tier1_genes.pdf",
  "../phase6_scrna/FigS14_microenvironment_elbow_plot.pdf" = "FigS16a_microenvironment_elbow_plot.pdf",
  "../figures/FigS15_me_zone_barplot.pdf"              = "FigS16b_me_zone_barplot.pdf",
  "../figures/FigS16_me_tier1_heatmap.pdf"             = "FigS16c_me_tier1_heatmap.pdf",
  "../figures/FigS22_finngen_replication_forest.pdf"   = "FigS17_finngen_replication_forest.pdf",
  "../figures/FigS18_idg_druggability.pdf"             = "FigS18_idg_druggability.pdf",
  "../figures/FigS19_phewas_safety_scores.pdf"         = "FigS19_phewas_safety_scores.pdf"
)

main_panels <- c(
  "Fig1b_study_design.pdf", "Fig1c_coloc_pph4.pdf",
  "Fig2a_coverage_funnel.pdf", "Fig2b_eqtl_pqtl_classification.pdf", "Fig2c_eqtl_pqtl_scatter.pdf",
  "Fig3a_ccm2_convergence.pdf", "Fig3b_gtex_colon_scatter.pdf", "Fig3c_gtex_colon_forest.pdf",
  "Fig4a_phewas_safety.pdf", "Fig4b_drug_repurposing.pdf",
  "Fig5a_competitor_overlap.pdf", "Fig5b_gene_overlap_venn.pdf"
)
miss <- main_panels[!file.exists(file.path(outdir, main_panels))]
if (length(miss)) stop("main panel not produced: ", paste(miss, collapse = ", "))

cat("\n-> final panel names\n")
for (nm in names(panel_names)) {
  a <- file.path(outdir, nm); b <- file.path(outdir, panel_names[[nm]])
  if (identical(normalizePath(a, mustWork = FALSE), normalizePath(b, mustWork = FALSE))) next
  if (!file.exists(a)) {
    if (file.exists(b)) { cat(sprintf("   %-36s (source missing, keeping existing %s)\n", nm, panel_names[[nm]])); next }
    stop("figure not produced: ", a)
  }
  if (!file.copy(a, b, overwrite = TRUE)) stop("copy failed: ", a, " -> ", b)
  cat(sprintf("   %-36s -> %s\n", nm, panel_names[[nm]]))
}
cat("\n== rev3 main figures done ==\n")
