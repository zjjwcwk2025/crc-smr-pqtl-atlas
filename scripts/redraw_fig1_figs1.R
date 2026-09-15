#!/usr/bin/env Rscript
# Redraw Fig1 (fixed text-overlap block) + FigS1a Manhattan + FigS1b QQ (CMplot)
suppressMessages({
  library(CMplot)
  library(data.table)
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
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================
#  Fig 1 — Study Design Overview (3-layer x 7-phase, clean layout)
# ==============================================================
cat("\n>>> Fig 1: Study Design Overview\n")

# Layer bands (subtle grouping; blue = discovery/replication, green = translation)
layer_bands <- tribble(
  ~xmin, ~xmax, ~ymin, ~ymax, ~fill,
  0.15, 5.85,  0.55,  1.55, "#EDF3FA",
  0.15, 5.85, -0.45,  0.45, "#EDF3FA",
  0.15, 5.85, -1.55, -0.55, "#EAF4EC"
)

phases <- tribble(
  ~phase, ~x, ~y, ~width, ~height, ~fill, ~label,
  1L, 1, 1, 1.7, 0.8, "#08519C", "SMR eQTL\u2192GWAS\n75 Bonf. significant",
  2L, 3, 1, 1.7, 0.8, "#08519C", "coloc + SuSiE\n43 PASS, 19 Tier 1",
  3L, 5, 1, 1.7, 0.8, "#08519C", "pQTL MR + coloc\nCCM2 triple-convergent",
  4L, 1, 0, 1.7, 0.8, "#2171B5", "FinnGen replication\n1/8 LIMA1 replicated",
  5L, 3, 0, 1.7, 0.8, "#2171B5", "GTEx Colon\ntissue sensitivity",
  6L, 5, 0, 1.7, 0.8, "#2171B5", "scRNA + Visium + ME\nk=6 zones, cross-correlation",
  7L, 3, -1, 2.7, 0.8, "#1B7837", "Drug repurposing\n+ IDG druggability (Tier 1)"
) %>%
  mutate(
    phase_label = paste0("Phase ", phase),
    xmin = x - width/2, xmax = x + width/2,
    ymin = y - height/2, ymax = y + height/2
  )

# Within-row flow arrows only (left -> right)
arrows <- tribble(
  ~x, ~y, ~xend, ~yend,
  1.85, 1, 2.15, 1,
  3.85, 1, 4.15, 1,
  1.85, 0, 2.15, 0,
  3.85, 0, 4.15, 0
)

data_sources <- tribble(
  ~x, ~y, ~label,
  6.25, 1.20, "eQTLGen (31,684 blood)\nGCST90255675 (78K CRC)",
  6.25, 0.80, "deCODE (4,907 proteins)\nUKB-PPP (2,923 proteins)",
  6.25, 0.00, "FinnGen R13 (10.5K CRC)",
  6.25, -1.00, "Pharos/DrugCentral\nOpen Targets Platform"
)

layer_labels <- tribble(
  ~x, ~y, ~label,
  -1.10, 1, "Discovery",
  -1.10, 0, "Replication &\nValidation",
  -1.10, -1, "Translation"
)

p1 <- ggplot() +
  geom_rect(data = layer_bands,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = NA) +
  geom_text(data = layer_labels, aes(x = x, y = y, label = label),
            size = 3.8, fontface = "bold", color = "grey25", lineheight = 0.9) +
  geom_rect(data = phases,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "grey25", linewidth = 0.5) +
  geom_text(data = phases,
            aes(x = x, y = y + 0.22, label = phase_label),
            size = 3.3, color = "white", fontface = "bold") +
  geom_text(data = phases,
            aes(x = x, y = y - 0.08, label = label),
            size = 2.9, color = "white", fontface = "bold", lineheight = 0.95) +
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = "grey45", linewidth = 0.7) +
  geom_text(data = data_sources, aes(x = x, y = y, label = label),
            size = 2.9, hjust = 0, color = "grey40", lineheight = 0.9) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-1.7, 8.7), ylim = c(-1.75, 1.75)) +
  labs(title = "CRC SMR Atlas: Multi-Omics Drug-Target Discovery Pipeline",
       subtitle = "78K CRC GWAS + eQTL + pQTL + spatial transcriptomics + drug annotation") +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(outdir, "Fig1_study_design.pdf"), p1, width = 12, height = 7, device = cairo_pdf)
cat("  -> Fig1_study_design.pdf\n")

# ==============================================================
#  Fig S1a/S1b — CMplot Manhattan + QQ (gene symbols, not ENSG IDs)
# ==============================================================
cat("\n>>> Fig S1a/S1b: CMplot Manhattan + QQ\n")
smr <- fread("results/phase1_eqtl_smr.smr")
setnames(smr, c("probeID","ProbeChr","Gene","Probe_bp","topSNP","topSNP_chr","topSNP_bp",
                "A1","A2","Freq","b_GWAS","se_GWAS","p_GWAS","b_eQTL","se_eQTL",
                "p_eQTL","b_SMR","se_SMR","p_SMR","p_HEIDI","nsnp_HEIDI"))

bio <- fread("results/phase1_ensg_biotype.tsv")
setnames(bio, c("ENSG","gene_type","gene_name"))

d <- data.frame(
  SNP        = smr$Gene,
  Chromosome = as.integer(smr$ProbeChr),
  Position   = as.integer(smr$Probe_bp),
  P.value    = as.numeric(smr$p_SMR),
  stringsAsFactors = FALSE
)
d <- d[!is.na(d$Chromosome) & !is.na(d$Position) & !is.na(d$P.value) & d$P.value > 0, ]

# genomic inflation factor (lambda)
chi2 <- qchisq(d$P.value, 1, lower.tail = FALSE)
lam  <- median(chi2, na.rm = TRUE) / qchisq(0.5, 1)
cat(sprintf("  lambda = %.6f  (display %.3f)\n", lam, lam))

# gene symbol / biotype mapping (iron rule #8: no silent bare-ENSG labels)
bio_map <- bio[, .(ENSG, gene_type, gene_name)]
d$gene_name <- bio_map$gene_name[match(d$SNP, bio_map$ENSG)]
d$gene_type <- bio_map$gene_type[match(d$SNP, bio_map$ENSG)]

make_label <- function(ensg, sym, biotype) {
  if (is.na(sym) || sym == "" || sym == ensg) {
    paste0(ensg, " (", biotype, ")")
  } else {
    sym
  }
}
d$label <- mapply(make_label, d$SNP, d$gene_name, d$gene_type)

threshold <- 0.05 / nrow(d)
sig <- d[d$P.value < threshold, ]
cat(sprintf("  n significant (p < %.3e) = %d\n", threshold, nrow(sig)))

topN <- min(30, nrow(sig))
top_hits <- head(sig[order(sig$P.value), ], topN)
cat("  Top-hit labels (gene symbol):\n")
print(top_hits[, c("SNP","label","P.value")], row.names = FALSE)

# CMplot requires a 4-column data.frame (SNP, Chromosome, Position, P.value);
# extra columns (gene_name/gene_type/label) make it mis-parse as multi-trait input.
d_cm <- d[, c("SNP","Chromosome","Position","P.value")]

# CMplot writes to getwd() and prepends a fixed prefix ("Rect_Manhtn." for a
# single-trait Manhattan, "QQplot." for a single-trait QQ). file.name must be a
# bare basename (no directory), so generate inside outdir and rename afterwards.
oldwd <- getwd()
setwd(outdir)

# --- Manhattan ---
# 省略 col 参数 -> 使用 CMplot 默认 10 色板（按染色体循环着色，横坐标彩色）
CMplot(d_cm, plot.type = "m",
       pch = 19,
       threshold = threshold, threshold.col = "red", threshold.lwd = 1, threshold.lty = 2,
       amplify = TRUE, signal.col = "red3", signal.pch = 19, signal.cex = 0.8,
       highlight = top_hits$SNP, highlight.col = "red3", highlight.pch = 19,
       highlight.text = top_hits$label,
       highlight.text.col = "black", highlight.text.cex = 0.65, highlight.text.font = 3,
       main = "SMR eQTL to CRC Manhattan plot",
       main.cex = 0.9, axis.cex = 0.7, lab.cex = 0.8, cex = 0.45,
       file.output = TRUE, file = "pdf",
       file.name = "FigS1a_manhattan",
       dpi = 300, width = 6.3, height = 2.7,
       verbose = FALSE)

# --- QQ ---
CMplot(d_cm, plot.type = "q",
       conf.int = TRUE, conf.int.col = "grey80", box = FALSE,
       main = paste0("QQ plot of SMR p-values (lambda = ", sprintf("%.3f", lam), ")"),
       main.cex = 0.9, axis.cex = 0.7, lab.cex = 0.8,
       file.output = TRUE, file = "pdf",
       file.name = "FigS1b_qq_plot",
       dpi = 300, width = 4.54, height = 4.54,
       verbose = FALSE)

# rename CMplot-prefixed files to clean names
file.rename("Rect_Manhtn.FigS1a_manhattan.pdf", "FigS1a_manhattan.pdf")
file.rename("QQplot.FigS1b_qq_plot.pdf", "FigS1b_qq_plot.pdf")
setwd(oldwd)

cat("  -> FigS1a_manhattan.pdf + FigS1b_qq_plot.pdf\n")

# --- postcondition checks (iron rule #5) ---
targets <- c(
  "Fig1_study_design.pdf",
  "FigS1a_manhattan.pdf",
  "FigS1b_qq_plot.pdf"
)
cat("\n>>> Postcondition check\n")
ok <- TRUE
for (f in targets) {
  p <- file.path(outdir, f)
  exists <- file.exists(p)
  size <- if (exists) file.info(p)$size else 0
  cat(sprintf("  %-32s exists=%s size=%d\n", f, exists, size))
  if (!exists || size == 0) ok <- FALSE
}
cat(if (ok) "  ALL OK\n" else "  FAILED\n")
