#!/usr/bin/env Rscript
# ============================================================
#  rev8_fig5.R -- Figure 5, redrawn to CNS art standard (rev8).
#  Panel sizes identical to the previous round.
#  5a: the fill no longer encodes journal impact factor (which the caption
#      never mentioned). It now encodes the two study tiers the caption does
#      name, so the legend and the caption use the same words. The stale
#      label "Wang suppl / BMC Cancer" is corrected to the name the
#      manuscript uses ("Xia & Wang 2025, BMC Cancer").
#  5b: same Venn geometry and the same numbers as before; palette unified
#      with the rest of the paper and type set from one theme.
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2); library(dplyr); library(readr)
  library(tidyr); library(tibble); library(forcats); library(stringr)
  library(ggforce)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
source("scripts/theme_pub2.R")
outdir <- "results/figures_rev8_cns"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
th <- theme_pub2()

# ============================================================
#  Figure 5a -- overlap with six prior CRC target-discovery studies
# ============================================================
cat("\n>>> Figure 5a: overlap with prior studies\n")
comp <- read_csv("results/competitor_overlap.tsv", show_col_types = FALSE)
v5_sym <- unique(read_tsv("results/phase1_bonferroni_significant_annotated.tsv",
                          show_col_types = FALSE)$SYMBOL)
v5_sym <- sort(unique(v5_sym[!is.na(v5_sym) & v5_sym != ""]))
stopifnot(length(v5_sym) == 69)

t9 <- readLines("manuscript/tables/TableS9_chen2024_gene_list.tex")
t9 <- t9[grepl("^\\s*[A-Za-z0-9][A-Za-z0-9._-]*\\s*&", t9)]
chen <- sort(unique(sub("^\\s*([A-Za-z0-9._-]+)\\s*&.*$", "\\1", t9)))
chen <- chen[chen != "Gene"]
n_chen_ov <- sum(v5_sym %in% chen)
stopifnot(length(chen) == 289, n_chen_ov == 31)
comp <- bind_rows(comp, tibble(
  study = "Chen_2024_NatCommun", journal_IF = "16",
  total_genes = length(chen), overlap_n = n_chen_ov,
  overlap_rate = sprintf("%.1f%%", 100 * n_chen_ov / length(v5_sym)),
  overlapping_genes = ""))
comp <- comp %>% filter(study != "ALL_COMBINED") %>% filter(!is.na(overlap_n))
stopifnot(nrow(comp) == 6)

HI <- c("Chen_2024_NatCommun", "Hazelwood_2025_NatCommun")
comp <- comp %>%
  mutate(
    short = case_when(
      grepl("^Hong", study)       ~ "Hong 2025, Front Immunol",
      grepl("^Wang_2025", study)  ~ "Wang 2025, Disc Oncol",
      grepl("^Wang_suppl", study) ~ "Xia & Wang 2025, BMC Cancer",
      grepl("^Tian", study)       ~ "Tian 2025, Biomedicines",
      grepl("^Chen", study)       ~ "Chen 2024, Nat Commun",
      grepl("^Hazelwood", study)  ~ "Hazelwood 2025, Nat Commun",
      TRUE ~ study),
    tier = factor(ifelse(study %in% HI, "High-tier studies", "Low-tier studies"),
                  levels = c("Low-tier studies", "High-tier studies")),
    short = fct_reorder(short, overlap_n))
stopifnot(n_distinct(comp$short) == 6, sum(comp$tier == "High-tier studies") == 2)

p5a <- ggplot(comp, aes(y = short)) +
  geom_rect(aes(xmin = 0, xmax = overlap_n, ymin = as.numeric(short) - 0.30,
                ymax = as.numeric(short) + 0.30, fill = tier), colour = NA) +
  scale_fill_manual(name = "Prior study",
                    values = c(`Low-tier studies` = "#9DC3E6",
                               `High-tier studies` = COL$genetic)) +
  geom_text(aes(x = overlap_n, label = sprintf("  %d  (%s)", overlap_n, overlap_rate)),
            hjust = 0, size = sz(6.6), colour = COL$ink, fontface = "bold") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = 0.65)) +
  coord_cartesian(xlim = c(0, max(comp$overlap_n) * 1.35), clip = "off") +
  labs(x = "Overlapping gene symbols (this study = 69 unique symbols, denominator of every %)",
       y = NULL) +
  th + theme(panel.grid.major.y = element_blank(),
             legend.position = "bottom", legend.box = "horizontal",
             plot.margin = margin(4, 8, 3, 6))
save_panel(p5a, file.path(outdir, "Fig5a_competitor_overlap.pdf"),
           mm_snap(169.69), mm_snap(94.19))

# ============================================================
#  Figure 5b -- shared / unique gene symbols (Venn)
# ============================================================
cat("\n>>> Figure 5b: gene overlap Venn\n")
haz <- sort(unique(c("AAMP","ABCC2","AC087392.1","AC144831.1","AL121832.3","ARPC2","ATF1",
  "CCM2","CENPBD2P","COLCA1","COX14","COX15","DACT1","EPM2AIP1","FADS1","FEN1","GPATCH1",
  "KLF5","LAMC1","LAMC1-AS1","LIMA1","LRRFIP2","METRNL","MLH1","MZF1","MZF1-AS1","PLEKHG6",
  "PNKD","POU2AF2","POU2AF3","POU5F1B","RBBP8NL","RP11-129K12.1","RPS21-DT","SEMA4D","TCF19",
  "TMBIM1","GPBAR1","LTBR","PDCD1","PTGER3")))
stopifnot(length(haz) == 41, length(chen) == 289, length(v5_sym) == 69)

allg <- unique(c(v5_sym, chen, haz))
mem <- data.frame(gene = allg, Chen = allg %in% chen, v5 = allg %in% v5_sym, Haz = allg %in% haz)
cnt <- function(a, b, cc) sum(mem$Chen == a & mem$v5 == b & mem$Haz == cc)
n_v5 <- length(v5_sym); n_chen <- length(chen); n_haz <- length(haz)

BLUE <- COL$genetic; CYAN <- "#6BAED6"; GREY <- COL$mute
SELF <- "This study"
DNAM <- c("Chen 2024", SELF, "Hazelwood 2025")
DSET <- c(n_chen, n_v5, n_haz)

up <- data.frame(
  lab = c("Chen 2024 only", "This study only", "This study + Chen 2024",
          "Chen 2024 + Hazelwood 2025", "Hazelwood 2025 only", "All three studies",
          "This study + Hazelwood 2025"),
  n   = c(cnt(TRUE, FALSE, FALSE), cnt(FALSE, TRUE, FALSE), cnt(TRUE, TRUE, FALSE),
          cnt(TRUE, FALSE, TRUE), cnt(FALSE, FALSE, TRUE), cnt(TRUE, TRUE, TRUE),
          cnt(FALSE, TRUE, TRUE)),
  Chen = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE),
  v5   = c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
  Haz  = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE)
stopifnot(sum(up$n) == length(allg),
          cnt(TRUE, FALSE, FALSE) == 242, cnt(FALSE, TRUE, FALSE) == 34,
          cnt(TRUE, TRUE, FALSE) == 23, cnt(TRUE, FALSE, TRUE) == 16,
          cnt(FALSE, FALSE, TRUE) == 13, cnt(TRUE, TRUE, TRUE) == 8,
          cnt(FALSE, TRUE, TRUE) == 4)
up <- up[order(-up$n), ]
up$yp <- rev(seq_len(nrow(up)))
up$xend <- up$n / max(up$n) * 0.50

DX <- c(0.620, 0.765, 0.910)
YTOP <- nrow(up)

p5b <- ggplot(up) +
  geom_segment(aes(x = DX[1], xend = DX[3], y = yp, yend = yp),
               colour = "#DCDCDC", linewidth = 0.35) +
  geom_rect(aes(xmin = 0, xmax = xend, ymin = yp - 0.29, ymax = yp + 0.29),
            fill = BLUE, colour = NA) +
  geom_text(aes(x = xend + 0.014, y = yp, label = n), hjust = 0,
            size = sz(6.2), colour = COL$ink, fontface = "bold") +
  geom_text(aes(x = -0.016, y = yp, label = lab), hjust = 1,
            size = sz(6.0), colour = COL$ink2) +
  geom_point(aes(x = DX[1], y = yp, colour = Chen), size = 2.6, stroke = 0.7) +
  geom_point(aes(x = DX[2], y = yp, colour = v5),   size = 2.6, stroke = 0.7) +
  geom_point(aes(x = DX[3], y = yp, colour = Haz),  size = 2.6, stroke = 0.7) +
  scale_colour_manual(values = c(`TRUE` = BLUE, `FALSE` = "#DCDCDC"), guide = "none") +
  annotate("text", x = DX, y = YTOP + 0.92, hjust = 0.5, vjust = 1,
           size = sz(6.0), fontface = "bold", colour = COL$ink, lineheight = 1.10,
           label = sprintf("%s\nn = %d", DNAM, DSET)) +
  annotate("text", x = 0, y = YTOP + 0.92, hjust = 0, vjust = 1,
           size = sz(6.0), fontface = "bold", colour = COL$ink,
           label = "Number of genes") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.34, 1.00), ylim = c(0.35, YTOP + 1.55), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "Liberation Sans") +
  theme(plot.margin = margin(3, 4, 3, 4))
save_panel(p5b, file.path(outdir, "Fig5b_gene_overlap_upset.pdf"),
           mm_snap(169.69), mm_snap(68))
cat("\n== rev8 Figure 5 done ==\n")
