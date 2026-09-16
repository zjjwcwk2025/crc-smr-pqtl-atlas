#!/usr/bin/env Rscript
# =====================================================================
#  Per-panel regeneration: figures whose text was below the 7 pt floor
#  at final print size.
#    FigS1a_manhattan        gene labels were drawn at 0.55 cex (~3.6 pt)
#    FigS16b_me_zone_barplot rendered at 4.41 in but printed at 2.90 in
#  Both are written to results/figures/ (the canonical generator output)
#  and copied into results/figures_rev3/ under their final figure names.
# =====================================================================
suppressMessages({
  library(data.table); library(ggplot2); library(dplyr); library(readr)
  library(CMplot)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
outdir <- "results/figures_rev3"
srcdir <- "results/figures"
source("scripts/theme_pub.R")

# ---------------- Fig S1a : Manhattan, readable gene labels ----------------
smr <- fread("results/phase1_eqtl_smr.smr")
setnames(smr, c("probeID","ProbeChr","Gene","Probe_bp","topSNP","topSNP_chr","topSNP_bp",
                "A1","A2","Freq","b_GWAS","se_GWAS","p_GWAS","b_eQTL","se_eQTL",
                "p_eQTL","b_SMR","se_SMR","p_SMR","p_HEIDI","nsnp_HEIDI"))
bio <- fread("results/phase1_ensg_biotype.tsv"); setnames(bio, c("ENSG","gene_type","gene_name"))
d <- data.frame(SNP = smr$Gene, Chromosome = as.integer(smr$ProbeChr),
                Position = as.integer(smr$Probe_bp), P.value = as.numeric(smr$p_SMR),
                stringsAsFactors = FALSE)
d <- d[!is.na(d$Chromosome) & !is.na(d$Position) & !is.na(d$P.value) & d$P.value > 0, ]
d$gene_name <- bio$gene_name[match(d$SNP, bio$ENSG)]
d$gene_type <- bio$gene_type[match(d$SNP, bio$ENSG)]
mk <- function(ensg, sym, bt) if (is.na(sym) || sym == "" || sym == ensg) paste0(ensg, " (", bt, ")") else sym
d$label <- mapply(mk, d$SNP, d$gene_name, d$gene_type)
threshold <- 0.05 / nrow(d)
sig <- d[d$P.value < threshold, ]
top_hits <- head(sig[order(sig$P.value), ], min(12, nrow(sig)))  # top 12 keeps labels legible
d_cm <- d[, c("SNP","Chromosome","Position","P.value")]

old <- getwd(); setwd(srcdir)
CMplot(d_cm, plot.type = "m", pch = 19,
       threshold = threshold, threshold.col = "red", threshold.lwd = 1, threshold.lty = 2,
       amplify = TRUE, signal.col = "red3", signal.pch = 19, signal.cex = 0.8,
       highlight = top_hits$SNP, highlight.col = "red3", highlight.pch = 19,
       highlight.text = top_hits$label,
       highlight.text.col = "black",
       highlight.text.cex = 0.78,          # ~8 pt effective at 6.3 in print width
       highlight.text.font = 3,
       axis.cex = 0.7, lab.cex = 0.8, cex = 0.45,
       file.output = TRUE, file = "pdf", file.name = "FigS1a_manhattan",
       dpi = 300, width = 6.3, height = 3.4, verbose = FALSE)
file.rename("Rect_Manhtn.FigS1a_manhattan.pdf", "FigS1a_manhattan.pdf")
setwd(old)
cat("  -> results/figures/FigS1a_manhattan.pdf (gene labels ~7.2 pt)\n")

# ---------------- Fig S16b : ME zone barplot at final print width ----------
me_chars <- read_csv("results/phase6_scrna/TableS_microenvironment_characterization.csv",
                     show_col_types = FALSE) %>%
  mutate(me_label = factor(microenvironment, levels = paste0("ME", 1:6)))
pS6a <- ggplot(me_chars, aes(x = me_label, y = n_spots, fill = top_ct1)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set2", name = "Dominant cell type") +
  scale_y_continuous(limits = c(0, max(me_chars$n_spots) * 1.06), expand = c(0, 0)) +
  labs(x = "", y = "Visium spots") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 7.5),
        axis.text.y = element_text(size = 7.0),
        axis.title.y = element_text(size = 8, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7.0),
        legend.key.size = unit(0.20, "cm"),
        legend.spacing.x = unit(0.06, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt")) +
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5))
ggsave(file.path(srcdir, "FigS15_me_zone_barplot.pdf"), pS6a,
       width = 2.90, height = 2.60, device = cairo_pdf)
cat("  -> results/figures/FigS15_me_zone_barplot.pdf (2.90 in print width)\n")

# ---------------- deliver into results/figures_rev3 -----------------------
for (p in list(c("FigS1a_manhattan.pdf","FigS1a_manhattan.pdf"),
               c("FigS15_me_zone_barplot.pdf","FigS16b_me_zone_barplot.pdf"))) {
  ok <- file.copy(file.path(srcdir, p[1]), file.path(outdir, p[2]), overwrite = TRUE)
  cat(sprintf("   %-34s -> %s  (%s)\n", p[1], p[2], if (ok) "ok" else "FAILED"))
}
cat("== 71_regen_small_text done ==\n")
