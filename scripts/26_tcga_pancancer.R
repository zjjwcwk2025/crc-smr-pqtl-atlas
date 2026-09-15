#!/usr/bin/env Rscript
#'
#' Phase 5c-7: TCGA Expression Analysis of 6 Tier 1 genes
#' Uses existing COAD+READ data (pan-cancer download too slow from China)
#' Plots tumor vs normal expression for Tier 1 genes
#'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

PROJ <- "/ifs1/User/zhouman/project9-v5-crc-atlas"
setwd(PROJ)
dir.create("results/phase5c_pancancer", showWarnings = FALSE, recursive = TRUE)

# ── 1. Load Tier 1 genes ──
susie <- fread("results/phase3_susie_finemap.csv")
tier1_genes <- susie[tier == 1, gene]
cat("Tier 1 genes:", paste(tier1_genes, collapse = ", "), "\n")

# ── 2. Load TCGA COAD+READ data ──
cat("\n[1/2] Loading TCGA COAD+READ expression...\n")
tcga <- fread("data/tcga_coad/TCGA_COADREAD_HiSeqV2.gz")

# First column is gene symbol
setnames(tcga, names(tcga)[1], "gene")

# Check available genes
avail_genes <- tier1_genes[tier1_genes %in% tcga$gene]
missing_genes <- tier1_genes[!(tier1_genes %in% tcga$gene)]
cat(sprintf("  %d/%d Tier 1 genes found in TCGA COAD+READ\n", length(avail_genes), length(tier1_genes)))
if (length(missing_genes) > 0) cat(sprintf("  Missing: %s\n", paste(missing_genes, collapse = ", ")))

tcga_sub <- tcga[gene %in% avail_genes]
cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(tcga_sub), ncol(tcga_sub) - 1))

# ── 3. Parse sample types ──
# TCGA barcode: TCGA-XX-YYYY-ZZ
# ZZ = 01: primary tumor, 11: solid normal
sample_ids <- names(tcga_sub)[-1]
sample_type <- ifelse(grepl("-01[A-Z]?$", sample_ids), "Tumor",
                      ifelse(grepl("-11[A-Z]?$", sample_ids), "Normal", "Other"))

# Melt to long
tcga_long <- melt(tcga_sub, id.vars = "gene", variable.name = "sample", value.name = "expression")
tcga_long$type <- sample_type[match(tcga_long$sample, sample_ids)]
tcga_long <- tcga_long[type %in% c("Tumor", "Normal")]

# Expression is log2(RSEM+1) — already normalized
cat(sprintf("  %d Tumor samples, %d Normal samples\n",
            sum(tcga_long$type == "Tumor") / length(avail_genes),
            sum(tcga_long$type == "Normal") / length(avail_genes)))

# ── 4. Boxplot ──
cat("\n[2/2] Generating plots...\n")

p <- ggplot(tcga_long, aes(x = type, y = expression, fill = type)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.8, width = 0.5) +
  geom_jitter(width = 0.1, size = 0.2, alpha = 0.15) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Tumor" = "#E41A1C", "Normal" = "#377EB8"),
                    name = "") +
  labs(x = "", y = "Expression (log2 RSEM+1)",
       title = "Tier 1 Gene Expression: TCGA COAD+READ",
       subtitle = "380 Tumor vs 51 Normal samples") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave("results/phase5c_pancancer/tcga_coadread_expression.pdf", p, width = 10, height = 7)
ggsave("results/phase5c_pancancer/tcga_coadread_expression.png", p, width = 10, height = 7, dpi = 150)

# ── 5. Differential expression tests ──
cat("\n── Differential Expression (Wilcoxon) ──\n")

de_results <- data.frame(
  gene = character(), tumor_mean = numeric(), normal_mean = numeric(),
  log2FC = numeric(), p_value = numeric(), significance = character(),
  stringsAsFactors = FALSE
)

for (g in avail_genes) {
  gdata <- tcga_long[gene == g]
  tumor_vals <- gdata[type == "Tumor", expression]
  normal_vals <- gdata[type == "Normal", expression]
  if (length(normal_vals) < 3) next

  wt <- wilcox.test(tumor_vals, normal_vals)
  fc <- mean(tumor_vals) - mean(normal_vals)  # log2 fold change

  p_label <- sprintf("p = %.2e", wt$p.value)
  sig <- ifelse(wt$p.value < 0.001, "p < 0.001",
                ifelse(wt$p.value < 0.01, "p < 0.01",
                       ifelse(wt$p.value < 0.05, "p < 0.05", "ns")))

  de_results <- rbind(de_results, data.frame(
    gene = g,
    tumor_mean = mean(tumor_vals),
    normal_mean = mean(normal_vals),
    log2FC = fc,
    p_value = wt$p.value,
    significance = sig,
    p_label = p_label,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %s: Tumor=%.2f vs Normal=%.2f, log2FC=%+.2f, %s (%s)\n",
              g, mean(tumor_vals), mean(normal_vals), fc, p_label, sig))
}

fwrite(de_results, "results/phase5c_pancancer/tcga_coadread_de.csv")

# ── 6. Highlight plot with significance ──
de_dt <- as.data.table(de_results)
tcga_long_sig <- merge(as.data.table(tcga_long), de_dt[, .(gene, significance)], by = "gene")

p2 <- ggplot(de_results, aes(x = reorder(gene, -log2FC), y = log2FC,
                              fill = significance)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = p_label, vjust = ifelse(log2FC >= 0, -0.4, 1.4)),
            size = 3, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("p < 0.001" = "#E41A1C", "p < 0.01" = "#FF7F00",
                                "p < 0.05" = "#377EB8", "ns" = "grey70"),
                    name = "Two-sided Wilcoxon p") +
  labs(x = "", y = expression(log[2] ~ "Fold Change (Tumor / Normal)"),
       title = "Tier 1 Genes: Tumor vs Normal Expression",
       subtitle = "TCGA COAD+READ (380 tumor, 51 normal)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))

ggsave("results/phase5c_pancancer/tcga_log2fc_barplot.pdf", p2, width = 8, height = 5)
ggsave("results/phase5c_pancancer/tcga_log2fc_barplot.png", p2, width = 8, height = 5, dpi = 150)

cat("\nDone. Results saved to results/phase5c_pancancer/\n")
