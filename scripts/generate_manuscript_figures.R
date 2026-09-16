#!/usr/bin/env Rscript
# generate_manuscript_figures.R
# Generate ALL individual panel figures for the CRC SMR Atlas manuscript.
# Output: results/figures/Fig*_*.pdf — one PDF per panel, publication quality.
# Does NOT assemble composite figures (user will do that manually).

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(forcats)
  library(stringr)
  library(ggrepel)
  library(scales)
  library(ComplexUpset)
  library(cowplot)
})

setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "results/figures"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# ── ggplot2 theme ──────────────────────────────────────────────
# Figures are generated at their final printed width (see manuscript), so these
# sizes are the on-page point sizes: nothing falls below ~8 pt after inclusion.
theme_nc <- theme_classic(base_size=9, base_family="Liberation Sans") +
  theme(
    axis.title = element_text(size=10, face="bold"),
    axis.text  = element_text(size=8.5, color="black"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.title = element_text(size=9, face="bold"),
    legend.text = element_text(size=8),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color="grey92", linewidth=0.25),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(5,6,4,4)
  )

# ==== Colours ====
palette_coloc <- c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#CA0020")
palette_class <- c(
  "A: pQTL concordant" = "#2166AC",
  "B: pQTL discordant" = "#B2182B",
  "C: pQTL not significant" = "#999999"
)
palette_direction <- c("concordant" = "#2166AC", "discordant" = "#B2182B", "not_significant" = "#999999")
palette_safety <- c("High alert" = "#B2182B", "Moderate" = "#F4A582", "Low" = "#92C5DE")

# ==============================================================
#  Fig 1 — Study Design Overview (8-phase pipeline schematic)
# ==============================================================
cat("\n>>> Fig 1: Study Design Overview\n")
# Build an 8-phase flow chart as a ggplot with rectangles and arrows
phases <- tribble(
  ~phase, ~label, ~x, ~y, ~width, ~height, ~fill,
  1L, "Phase 1\nSMR eQTL→GWAS\n75 Bonf. significant", 1, 1, 1.6, 0.75, "#2166AC",
  2L, "Phase 2\ncoloc + SuSiE\n43 PASS, 19 Tier 1", 3, 1, 1.6, 0.75, "#4393C3",
  3L, "Phase 3\npQTL MR + coloc\nCCM2 triple-convergent", 5, 1, 1.6, 0.75, "#92C5DE",
  4L, "Phase 4\nFinnGen replication\n1/8 LIMA1 replicated", 1, 0, 1.6, 0.75, "#D1E5F0",
  5L, "Phase 5\nGTEx Colon\ntissue sensitivity", 3, 0, 1.6, 0.75, "#D1E5F0",
  6L, "Phase 6\nscRNA + Visium + ME\nk=6 zones, cross-correlation", 5, 0, 1.6, 0.75, "#D1E5F0",
  7L, "Phase 7\nDrug repurposing\n+Tier 1 druggability", 1, -1, 1.6, 0.75, "#FDDBC7",
  8L, "Phase 8\nCompetitor landscape\nChen 2024, Hazelwood 2025", 3, -1, 1.6, 0.75, "#FDDBC7"
) %>%
  mutate(
    phase_label = paste0("Phase ", phase),
    xmin = x - width/2, xmax = x + width/2,
    ymin = y - height/2, ymax = y + height/2
  )

# Arrows between phases
arrows <- tribble(
  ~x, ~y, ~xend, ~yend,
  2, 1, 2.2, 1,  # 1→2
  4, 1, 4.2, 1,  # 2→3
  1.8, 0.625, 1.8, 0.375,  # 1→4 (down)
  1, 0, 2.2, 0,   # 4→5
  4, 0, 4.2, 0,   # 5→6
  3.8, -0.375, 3.8, -0.625,  # 5→8 (down)
  1, -1, 2.2, -1 # 7→8
)

# Border box
border_box <- data.frame(
  xmin = rep(c(0.1, 0.1, 5.9), 2),
  xmax = rep(c(5.9, 0.1, 5.9), 2),
  ymin = c(-1.8, -1.8, -1.8, 1.8, 1.8, -1.8),
  ymax = c(-1.8, 1.8, 1.8, 1.8, -1.8, -1.8),
  group = rep(1:2, each = 3)
)

# Data sources annotation
data_sources <- tribble(
  ~x, ~y, ~label,
  6.45, 0.55, "eQTLGen (31,684 blood)\nGCST90255675 (78K CRC)",
  6.45, 0.05, "deCODE (4,907 proteins)\nUKB-PPP (2,923 proteins)",
  6.45, -0.45, "FinnGen R13 (10.5K CRC)",
  6.45, -1.45, "Pharos/DrugCentral\nOpen Targets Platform"
)

# Layer labels
layer_labels <- tribble(
  ~x, ~y, ~label,
  -0.6, 1, "Discovery",
  -0.6, 0, "Replication &\nValidation",
  -0.6, -1, "Translation"
)

p1 <- ggplot() +
  # Layer separators
  geom_hline(yintercept = c(0.55, -0.55), linetype = "dashed", color = "grey60", linewidth = 0.4) +
  # Layer labels
  geom_text(data = layer_labels, aes(x = x, y = y, label = label),
            size = 3.8, fontface = "bold", color = "grey30", lineheight = 0.9) +
  # Phase boxes
  geom_rect(data = phases,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "grey30", linewidth = 0.5, alpha = 0.92) +
  # Phase labels inside boxes
  geom_text(data = phases,
            aes(x = x, y = y, label = label),
            size = 2.6, color = "white", fontface = "bold", lineheight = 0.9) +
  # Phase numbers on top
  geom_text(data = phases,
            aes(x = x, y = ymax + 0.1, label = phase_label),
            size = 2.8, fontface = "bold", color = "grey20") +
  # Arrows
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = "grey40", linewidth = 0.6) +
  # Vertical connector arrow (Phase 3→6)
  geom_segment(aes(x = 5, y = 0.625, xend = 5, yend = 0.375),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = "grey40", linewidth = 0.6) +
  # Data source annotations
  geom_text(data = data_sources, aes(x = x, y = y, label = label),
            size = 2.2, hjust = 0, color = "grey40", lineheight = 0.9) +
  # Global arrow (bottom)
  geom_segment(aes(x = 5.9, y = -1.6, xend = 5.9, yend = 1.6),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed", ends = "both"),
               color = "grey50", linewidth = 0.8) +
  annotate("text", x = 5.9, y = 0, label = "8-phase\nanalytical\npipeline",
           angle = 90, size = 3, hjust = 0.5, color = "grey40", vjust = -1.5) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-1.3, 7.8), ylim = c(-1.9, 1.9)) +
  labs(title = "CRC SMR Atlas: Multi-Omics Drug-Target Discovery Pipeline",
       subtitle = "78K CRC GWAS + eQTL + pQTL + spatial transcriptomics + drug annotation") +
  theme_void() +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(outdir, "Fig1_study_design.pdf"), p1, width = 12, height = 7, device = cairo_pdf)
cat("  -> Fig1_study_design.pdf\n")

# ==============================================================
#  Fig 2a — coloc PPH4 barplot (phase3_coloc_all75.csv)
# ==============================================================
cat("\n>>> Fig 2a: coloc PPH4 barplot (all 76 probe-gene pairs)\n")
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
    plot.caption = element_text(hjust = 0.5, size = 8)
  )

ggsave(file.path(outdir, "Fig2a_coloc_pph4_barplot.pdf"), p2a,
       width = 7, height = 7.5, device = cairo_pdf)
cat("  -> Fig2a_coloc_pph4_barplot.pdf\n")

# ==============================================================
#  Fig 2a supplement — SMR→OR forest plot for 6 non-MHC Tier 1 genes
# ==============================================================
cat("\n>>> Fig 2a-supp: SMR→OR forest plot (6 non-MHC Tier 1 genes)\n")
smr_data <- read_tsv("results/phase1_bonferroni_significant_annotated.tsv", show_col_types=FALSE)

tier1_non_mhc <- c("UTP11", "BMP2", "PNKD", "LAMC1", "TMBIM1", "GPBAR1")

smr_tier1 <- smr_data %>%
  filter(SYMBOL %in% tier1_non_mhc) %>%
  mutate(
    OR = exp(b_SMR),
    OR_lower = exp(b_SMR - 1.96 * se_SMR),
    OR_upper = exp(b_SMR + 1.96 * se_SMR),
    direction = ifelse(b_SMR > 0, "Risk", "Protective"),
    SYMBOL = factor(SYMBOL, levels = rev(tier1_non_mhc))
  ) %>%
  arrange(SYMBOL)

p2a_supp <- ggplot(smr_tier1, aes(x = OR, y = SYMBOL, color = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_errorbar(aes(xmin = OR_lower, xmax = OR_upper), width = 0.25, linewidth = 1.2) +
  geom_point(size = 4.5, shape = 18) +
  geom_text(aes(label = sprintf("OR=%.2f [%.2f, %.2f]", OR, OR_lower, OR_upper),
                x = OR_upper + 0.04), size = 3.0, hjust = 0, color = "black") +
  scale_color_manual(values = c("Risk" = "#B2182B", "Protective" = "#2166AC"), guide = "none") +
  scale_x_continuous(limits = c(min(smr_tier1$OR_lower) * 0.9, max(smr_tier1$OR_upper) * 1.75)) +
  labs(x = "Odds Ratio (per 1 SD increase in gene expression)", y = "",
       title = "Tier 1 Gene Effects on CRC Risk",
       subtitle = paste0("SMR-derived OR with 95% CI; 6 non-MHC coloc+SuSiE-confirmed genes"),
       caption = "n = 6 Tier 1 genes; error bars = 95% CI from SMR standard error.") +
  annotate("text", x = 1.32, y = 0.62, label = "OR = 1 (no effect)",
           color = "grey50", size = 3.0, hjust = 1) +
  theme_nc + theme(panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
                   panel.grid.major.y = element_blank())

ggsave(file.path(outdir, "FigS2_smr_or_forest.pdf"), p2a_supp, width = 5.04, height = 2.52, device = cairo_pdf)
cat("  -> FigS2_smr_or_forest.pdf\n")

# ==============================================================
#  Fig 2b — eQTL/pQTL discordance classification bar
# ==============================================================
cat("\n>>> Fig 2b: eQTL/pQTL classification bar\n")
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

ggsave(file.path(outdir,"Fig2b_eqtl_pqtl_classification.pdf"), p2b, width=8, height=3.5, device=cairo_pdf)
cat("  -> Fig2b_eqtl_pqtl_classification.pdf\n")

# ==============================================================
#  Fig 2c — eQTL/pQTL effect size scatter
# ==============================================================
cat("\n>>> Fig 2c: eQTL vs pQTL effect size scatter\n")
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
  geom_text_repel(aes(label=SYMBOL), size=3.2, max.overlaps=20, force=2, box.padding=0.5) +
  scale_fill_manual(values=palette_direction, name="Direction") +
  scale_size_identity() +
  labs(x="eQTL SMR b (blood)", y="pQTL Wald b",
       title="eQTL vs pQTL Effect Size Concordance",
       subtitle=paste0("n=", nrow(contrast_both), " genes with both eQTL and pQTL effect estimates; filled circles = pQTL p<0.05"),
       caption=stats_label) +
  theme_nc

ggsave(file.path(outdir,"Fig2c_eqtl_pqtl_scatter.pdf"), p2c, width=7, height=6, device=cairo_pdf)
cat("  -> Fig2c_eqtl_pqtl_scatter.pdf\n")

# ==============================================================
#  Fig 3 — GTEx Colon vs Blood SMR scatter
# ==============================================================
cat("\n>>> Fig 3: GTEx Colon SMR sensitivity\n")
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
  geom_text_repel(aes(label=Gene), size=2.8, max.overlaps=15, force=2, box.padding=0.3) +
  scale_fill_manual(values=c("Sig in both"="#2166AC", "Sig in GTEx only"="#F4A582",
                               "Sig in eQTLGen only"="#92C5DE", "Neither sig"="grey80"),
                    name="Significance") +
  labs(x="SMR b (eQTLGen blood)", y="SMR b (GTEx Colon Transverse)",
       title="Tissue Sensitivity: GTEx Colon vs eQTLGen Blood",
       subtitle=paste0("Spearman rho = ",round(r_rho,2),", p = ",formatC(r_pval,format="f",digits=4),
                        "; direction concordant: ",concordant_n,"/",nrow(merged)," (",round(concordant_pct,1),"%)"),
       caption=paste0("n = ", nrow(merged), " genes with SMR estimates in both tissues; Spearman correlation, two-sided.")) +
  theme_nc

ggsave(file.path(outdir,"Fig3_gtex_colon_scatter.pdf"), p3, width=7.5, height=6.5, device=cairo_pdf)
cat("  -> Fig3_gtex_colon_scatter.pdf\n")

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
                   aes(xmin=ci_low, xmax=ci_high), width=0.3, linewidth=1.0) +
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

  ggsave(file.path(outdir,"Fig3_gtex_colon_forest.pdf"), p3f, width=7, height=5, device=cairo_pdf)
  cat("  -> Fig3_gtex_colon_forest.pdf\n")
}

# ==============================================================
#  Fig 4a — Drug annotation dot plot  (phase5b_deep_drug_annotation.csv)
# ==============================================================
cat("\n>>> Fig 4a: Drug annotation summary\n")
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

ggsave(file.path(outdir,"Fig4a_drug_repurposing.pdf"), p5a, width=9, height=7, device=cairo_pdf)
cat("  -> Fig4a_drug_repurposing.pdf\n")

# ==============================================================
#  Fig 5a supplement — IDG druggability gradient for Tier 1 genes
# ==============================================================
cat("\n>>> Fig 5a-supp: IDG druggability gradient (Tier 1 genes)\n")
idg_data <- tibble(
  gene = c("LAMC1", "GPBAR1", "BMP2", "UTP11", "TMBIM1", "PNKD"),
  tdl = c("Tclin", "Tchem", "Tchem", "Tchem*", "Tbio", "Tbio"),
  tdl_score = c(4, 3, 3, 2.5, 2, 1),  # Tclin > Tchem > Tchem*(no AB) > Tbio
  has_drug = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  antibody = c(3, 3, 3, 0, 3, 3),  # AB accessibility score
  sm = c("None", "Structure+Ligand\nDruggable Family", "Structure with Ligand",
         "Structure with Ligand", "None", "None")
) %>%
  mutate(gene = factor(gene, levels = rev(c("LAMC1", "GPBAR1", "BMP2", "UTP11", "TMBIM1", "PNKD"))))

# TDL color scale
tdl_colors <- c("Tclin" = "#2166AC", "Tchem" = "#4393C3", "Tchem*" = "#92C5DE", "Tbio" = "#D1E5F0")

p5a_idg <- ggplot(idg_data, aes(x = gene, y = tdl_score, fill = tdl)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = tdl), hjust = -0.2, size = 4.5, fontface = "bold") +
  geom_hline(yintercept = c(1, 2, 3, 4), linetype = "dotted", color = "grey80", linewidth = 0.3) +
  scale_fill_manual(values = tdl_colors, name = "IDG TDL") +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0), breaks = 1:4,
                     labels = c("Tbio", "Tbio+", "Tchem", "Tclin")) +
  labs(x = "", y = "Druggability Level",
       title = "Tier 1 Gene Druggability Gradient (IDG Pharos/DrugCentral)",
       subtitle = "Tclin = approved drug target; Tchem = chemical probe available; Tbio = biological only; *UTP11 has 0/3 AB accessibility",
       caption = "n = 6 Tier 1 genes; TDL and antibody accessibility from Pharos/DrugCentral.") +
  coord_flip() +
  theme_nc

ggsave(file.path(outdir, "FigS18_idg_druggability.pdf"), p5a_idg, width = 5.04, height = 2.52, device = cairo_pdf)
cat("  -> FigS18_idg_druggability.pdf\n")

# ==============================================================
#  Fig 4b — PheWAS safety heatmap  (phase5c_phewas_safety.csv)
# ==============================================================
cat("\n>>> Fig 4b: PheWAS safety heatmap\n")
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
    axis.text.x = element_text(angle=45, hjust=1, size=10),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(outdir,"Fig4b_phewas_safety_heatmap.pdf"), p5b, width=9, height=4.5, device=cairo_pdf)
cat("  -> Fig4b_phewas_safety_heatmap.pdf\n")

# Safety score barplot (companion)
p5b2 <- ggplot(phewas %>% mutate(gene=fct_reorder(gene, safety_concern_score)),
       aes(x=gene, y=safety_concern_score, fill=safety_tier)) +
  geom_col(width=0.6, color="black", linewidth=0.3) +
    geom_text(aes(label=paste0(sprintf("%.1f", safety_concern_score),
                             " (", n_genetic_diseases, " genetic)")),
            hjust=-0.05, size=2.6) +
  scale_fill_manual(values=palette_safety, guide="none") +
  scale_y_continuous(limits=c(0, max(phewas$safety_concern_score)*1.45), expand=c(0,0)) +
  labs(x="", y="Safety Concern Score",
       title="PheWAS Safety Scores for Tier 1 Genes",
       subtitle=paste0(nrow(phewas)," Tier 1 genes assessed via Open Targets Platform; score integrates disease count, genetic evidence, and risk category diversity"),
       caption="Safety score = 2.0×category diversity + 1.5×genetic evidence + 0.5×high-risk count.") +
  coord_flip() +
  theme_nc

ggsave(file.path(outdir,"FigS19_phewas_safety_scores.pdf"), p5b2, width = 4.41, height = 2.52, device=cairo_pdf)
cat("  -> FigS19_phewas_safety_scores.pdf\n")

# ==============================================================
#  Fig 4c — Competitor overlap  (competitor_overlap.tsv + UpSet)
# ==============================================================
cat("\n>>> Fig 4c: Competitor overlap\n")
comp <- read_csv("results/competitor_overlap.tsv", show_col_types=FALSE)

# Barplot of overlap rates
comp_plot <- comp %>%
  filter(study != "ALL_COMBINED") %>%
  mutate(study_short = case_when(
    grepl("Hong", study) ~ "Hong 2025\nFront Immunol",
    grepl("Wang_2025", study) ~ "Wang 2025\nDisc Oncol",
    grepl("Wang_suppl", study) ~ "Wang suppl\nBMC Cancer",
    grepl("Tian", study) ~ "Tian 2025\nBiomedicines",
    grepl("Hazelwood", study) ~ "Hazelwood 2025\nNat Commun"
  ),
  study_short = fct_reorder(study_short, overlap_n),
  journal_if_num = as.numeric(journal_IF))

p5c_bar <- ggplot(comp_plot, aes(x=study_short, y=overlap_n)) +
  geom_col(aes(fill=journal_if_num), width=0.6, color="black", linewidth=0.3) +
  geom_text(aes(label=paste0(overlap_n, " (", overlap_rate, ")")),
            vjust=-0.5, size=3.5, fontface="bold") +
  scale_fill_gradient(low="#FDDBC7", high="#2166AC", name="Journal IF") +
  scale_y_continuous(limits=c(0, max(comp_plot$overlap_n)*1.25), expand=c(0,0)) +
  labs(x="", y="Overlapping genes with v5",
       title="Competitor Gene Overlap Analysis",
       subtitle=paste0("v5 total: 69 unique genes. ",comp$overlap_n[comp$study=="ALL_COMBINED"],
                       " overlap (", comp$overlap_rate[comp$study=="ALL_COMBINED"],")"),
       caption="Overlap counts are gene-symbol based; journal impact factor shown by fill.") +
  coord_flip() +
  theme_nc

ggsave(file.path(outdir,"Fig4c_competitor_overlap.pdf"), p5c_bar, width=8, height=4.5, device=cairo_pdf)
cat("  -> Fig4c_competitor_overlap.pdf\n")

# UpSet plot: v5 vs Chen 2024 vs Hazelwood 2025 (3 sets, intersections computed programmatically)
cat("  Building UpSet plot...\n")
v5_genes <- unique(read_tsv("results/phase1_bonferroni_significant_annotated.tsv", show_col_types=FALSE)$SYMBOL)
v5_genes <- v5_genes[v5_genes != "" & !is.na(v5_genes)]

# Chen 2024 gene list parsed from Table S8 (exclude the literal "Gene" artifact row)
chen_lines <- readLines("manuscript/tables/TableS8_chen2024_gene_list.tex")
chen_rows <- grepl("^\\s*\\S+\\s*&\\s*(\\\\textbf\\{Yes\\}|---)\\s*\\\\\\\\", chen_lines)
chen_genes <- sub("^\\s*(\\S+)\\s*&.*$", "\\1", chen_lines[chen_rows])
chen_genes <- chen_genes[chen_genes != "Gene"]

hazelwood_all <- c("AAMP","ABCC2","AC087392.1","AC144831.1","AL121832.3","ARPC2","ATF1","CCM2","CENPBD2P","COLCA1","COX14","COX15","DACT1","EPM2AIP1","FADS1","FEN1","GPATCH1","KLF5","LAMC1","LAMC1-AS1","LIMA1","LRRFIP2","METRNL","MLH1","MZF1","MZF1-AS1","PLEKHG6","PNKD","POU2AF2","POU2AF3","POU5F1B","RBBP8NL","RP11-129K12.1","RPS21-DT","SEMA4D","TCF19","TMBIM1","GPBAR1","LTBR","PDCD1","PTGER3")

# Programmatic intersection counts (canonical: v5∩Chen=31, v5∩Hazelwood=12, triple=8)
v5_chen <- sum(v5_genes %in% chen_genes)
v5_haz  <- sum(v5_genes %in% hazelwood_all)
triple  <- sum(v5_genes %in% chen_genes & v5_genes %in% hazelwood_all)
v5_only <- sum(!(v5_genes %in% c(chen_genes, hazelwood_all)))
cat(sprintf("    v5=%d  v5&Chen=%d  v5&Hazelwood=%d  triple=%d  v5-only=%d\n",
            length(v5_genes), v5_chen, v5_haz, triple, v5_only))

all_genes <- unique(c(v5_genes, chen_genes, hazelwood_all))
upset_df <- data.frame(
  gene = all_genes,
  "v5 (this study)" = all_genes %in% v5_genes,
  "Chen 2024" = all_genes %in% chen_genes,
  "Hazelwood 2025" = all_genes %in% hazelwood_all,
  check.names = FALSE
)

p4c_upset <- upset(
  upset_df,
  intersect = c("v5 (this study)", "Chen 2024", "Hazelwood 2025"),
  name = "Gene set membership",
  width_ratio = 0.3,
  sort_intersections_by = "cardinality"
) +
  labs(title = "Gene Overlap: v5 vs Chen 2024 vs Hazelwood 2025",
       subtitle = sprintf("v5 = %d  |  v5 & Chen = %d  |  v5 & Hazelwood = %d  |  triple = %d",
                          length(v5_genes), v5_chen, v5_haz, triple)) +
  theme(plot.title = element_text(hjust=0.5, face="bold", size=14),
        plot.subtitle = element_text(hjust=0.5, size=10))

ggsave(file.path(outdir,"Fig4c_competitor_upset.pdf"), p4c_upset, width=8, height=7, device=cairo_pdf)
cat("  -> Fig4c_competitor_upset.pdf\n")

# ==============================================================
#  Supp Fig S2 — Suggestive genes (not done; skip)
# ==============================================================
cat("\n>>> Fig S2: Skipping (suggestive manhattan — data already in phase1 outputs)\n")

# ==============================================================
#  Supp Fig S3 — SuSiE sensitivity grid
# ==============================================================
cat("\n>>> Fig S3: SuSiE sensitivity\n")
susie_sens <- read_csv("results/susie_sensitivity/susie_sensitivity_grid.csv", show_col_types=FALSE)
if(nrow(susie_sens) > 0) {
  susie_sens <- susie_sens %>%
    mutate(param_label = paste0("SNPs=", snp_cap, ", L=", L_param, ", ", window_kb, "kb"))

  # Transposed grid: the parameter combinations sit on the y axis and the six genes
  # on the x axis, so the long parameter labels stay legible instead of overlapping
  # as rotated x tick labels. Cells with no converged run are drawn grey.
  susie_grid <- expand.grid(gene = unique(susie_sens$gene),
                            param_label = unique(susie_sens$param_label),
                            stringsAsFactors = FALSE)
  pS3 <- ggplot() +
    geom_tile(data = susie_grid, aes(x = gene, y = param_label),
              fill = "grey90", color = "white", linewidth = 0.5) +
    geom_tile(data = susie_sens, aes(x = gene, y = param_label, fill = factor(n_cs)),
              color = "white", linewidth = 0.5) +
    geom_text(data = susie_sens, aes(x = gene, y = param_label, label = n_cs),
              size = 3.0, fontface = "bold") +
    scale_fill_brewer(palette = "Blues", name = "N CS") +
    labs(x = "", y = "") +
    theme_nc + theme(
      axis.text.x = element_text(size = 8.5),
      axis.text.y = element_text(size = 7.5),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  ggsave(file.path(outdir,"FigS3_susie_sensitivity.pdf"), pS3, width = 6.3, height = 2.95, device=cairo_pdf)
  cat("  -> FigS3_susie_sensitivity.pdf\n")
}

# ==============================================================
#  Supp Fig S4 — UKB-PPP pQTL coloc
# ==============================================================
cat("\n>>> Fig S4: UKB-PPP pQTL coloc\n")
# Canonical UKB-PPP pQTL coloc = phase3c (5 genes, incl. TNF); phase3b superseded (2026-08-17)
ukb_coloc <- read_csv("results/phase3c_pqtl_coloc.csv", show_col_types=FALSE) %>%
  select(gene, PPH4)

if(nrow(ukb_coloc) > 0) {
  ukb_coloc <- ukb_coloc %>%
    mutate(
      coloc_pass = PPH4 > 0.8,
      pph4_label = case_when(
        PPH4 < 0.001 ~ "<0.001",
        PPH4 > 0.999 ~ ">0.999",
        TRUE ~ sprintf("%.4f", PPH4)
      )
    )

  pS4 <- ggplot(ukb_coloc, aes(x=reorder(gene, PPH4), y=PPH4, fill=coloc_pass)) +
    geom_col(width=0.5, color="black", linewidth=0.3) +
    geom_text(aes(label=pph4_label), hjust=-0.12, vjust=0.5, size=3.0, fontface="bold") +
    geom_hline(yintercept=0.8, linetype="dashed", color="#B2182B") +
    scale_fill_manual(values=c("TRUE"="#2166AC", "FALSE"="grey80"), guide="none") +
    scale_y_continuous(limits=c(0,1.30), expand=c(0,0)) +
    labs(x="", y="PPH4",
         title="UKB-PPP pQTL–GWAS Colocalization",
         subtitle=paste0(nrow(ukb_coloc)," UKB-PPP genes tested; 0/", nrow(ukb_coloc), " pass PPH4 > 0.8"),
         caption=paste0("n = ", nrow(ukb_coloc), " UKB-PPP genes; dashed line, PPH4 = 0.8.")) +
    coord_flip() + theme_nc

  ggsave(file.path(outdir,"FigS4_ukbppp_pqtl_coloc.pdf"), pS4, width = 3.78, height = 2.21, device=cairo_pdf)
  cat("  -> FigS4_ukbppp_pqtl_coloc.pdf\n")
}

# ==============================================================
#  Supp Fig S5 — deCODE pQTL coloc
# ==============================================================
cat("\n>>> Fig S5: deCODE pQTL coloc\n")
decode_coloc <- read_csv("results/phase2_decode/phase2_decode_pqtl_coloc.csv", show_col_types=FALSE)
if(nrow(decode_coloc) > 0) {
  decode_coloc <- decode_coloc %>%
    mutate(
      coloc_pass = PPH4 > 0.8,
      pph4_label = case_when(
        PPH4 < 0.001 ~ "<0.001",
        PPH4 > 0.999 ~ ">0.999",
        TRUE ~ sprintf("%.4f", PPH4)
      )
    )

  pS5 <- ggplot(decode_coloc, aes(x=reorder(gene, PPH4), y=PPH4, fill=coloc_pass)) +
    geom_col(width=0.5, color="black", linewidth=0.3) +
    geom_text(aes(label=pph4_label), hjust=-0.12, vjust=0.5, size=3.0) +
    geom_hline(yintercept=0.8, linetype="dashed", color="#B2182B") +
    scale_fill_manual(values=c("TRUE"="#2166AC", "FALSE"="grey60"), guide="none") +
    scale_y_continuous(limits=c(0,1.32), expand=c(0,0)) +
    labs(x="", y="PPH4",
         title="deCODE pQTL–GWAS Colocalization",
         subtitle=paste0(nrow(decode_coloc)," deCODE genes with cis-pQTL instruments; only CCM2 passes coloc"),
         caption=paste0("n = ", nrow(decode_coloc), " deCODE genes; CCM2 PPH4 = 0.9514.")) +
    coord_flip() + theme_nc

  ggsave(file.path(outdir,"FigS5_decode_pqtl_coloc.pdf"), pS5, width = 3.78, height = 2.21, device=cairo_pdf)
  cat("  -> FigS5_decode_pqtl_coloc.pdf\n")
}

# ==============================================================
#  Supp Fig S6 — ME zone spatial distribution + Tier 1 expression
#  (TableS_microenvironment_characterization.csv + tier1_microenvironment_expression.csv)
# ==============================================================
cat("\n>>> Fig S6: Microenvironment (ME) zone distribution\n")
me_chars <- read_csv("results/phase6_scrna/TableS_microenvironment_characterization.csv", show_col_types=FALSE)
me_chars <- me_chars %>%
  mutate(
    # Compact axis labels: the dominant cell type is carried by the fill legend,
    # and the zone sizes and descriptions are given in the figure caption.
    me_label = factor(microenvironment, levels = paste0("ME", 1:6))
  )

# Panel A: spot count barplot per ME zone
# Rendered at the final printed width (0.44 x 6.30 in = 2.77 in). Titles, the
# per-zone spot counts and the zone descriptions are carried by the figure
# caption, so the panel itself only carries the bars and the fill legend.
pS6a <- ggplot(me_chars, aes(x = me_label, y = n_spots, fill = top_ct1)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set2", name = "Dominant cell type") +
  scale_y_continuous(limits = c(0, max(me_chars$n_spots) * 1.06), expand = c(0, 0)) +
  labs(x = "", y = "Visium spots") +
  theme_nc +
  theme(axis.text.x = element_text(size = 7.5),
        axis.text.y = element_text(size = 7.0),
        axis.title.y = element_text(size = 8, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7.0),
        legend.key.size = unit(0.22, "cm"),
        legend.spacing.x = unit(0.08, "cm"),
        legend.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5))

ggsave(file.path(outdir, "FigS15_me_zone_barplot.pdf"), pS6a, width = 2.90, height = 2.60, device = cairo_pdf)
cat("  -> FigS15_me_zone_barplot.pdf\n")

# Panel B: Tier 1 gene expression heatmap across ME zones
me_expr <- read_csv("results/phase6_scrna/TableS_tier1_microenvironment_expression.csv", show_col_types=FALSE)

# Select key Tier 1 genes (6 non-MHC + a few MHC for context)
tier1_genes_show <- c("UTP11", "BMP2", "PNKD", "LAMC1", "TMBIM1", "GPBAR1",
                      "TNF", "LTA", "SF3A3", "CLDN4", "LEMD3", "CCM2", "GDPGP1", "ABHD12B", "METTL27")
me_expr_plot <- me_expr %>%
  filter(gene %in% tier1_genes_show) %>%
  mutate(
    gene = factor(gene, levels = rev(tier1_genes_show)),
    microenvironment = factor(microenvironment, levels = paste0("ME", 1:6))
  )

pS6b <- ggplot(me_expr_plot, aes(x = microenvironment, y = gene, fill = mean_expr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", mean_expr)), size = 3, fontface = "bold") +
  scale_fill_gradient(low = "#F7FBFF", high = "#08306B", name = "Mean SCT\nexpression") +
  labs(x = "Microenvironment Zone", y = "",
       title = "Tier 1 Gene Expression Across Microenvironment Zones",
       subtitle = "Mean SCT-normalized expression per Visium spot",
       caption = "n = 15 genes across 6 zones; mean of SCT-normalized counts.") +
  theme_nc + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1)) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

ggsave(file.path(outdir, "FigS16_me_tier1_heatmap.pdf"), pS6b, width = 5.54, height = 3.60, device = cairo_pdf)
cat("  -> FigS16_me_tier1_heatmap.pdf\n")

# ==== Summary ====
cat("\n=== ALL FIGURES GENERATED ===\n")
cat("Output directory:", file.path(getwd(), outdir), "\n")
cat("Files:\n")
for(f in sort(list.files(outdir, pattern="\\.pdf$"))) {
  cat(sprintf("  %s  (%s)\n", f,
              format(file.size(file.path(outdir,f)), big.mark=",")))
}
cat("\nDone.\n")
