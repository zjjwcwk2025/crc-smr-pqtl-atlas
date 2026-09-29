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
CX <- 2.75; CY <- 2.40; R <- 1.28; D <- 1.15
off <- data.frame(
  set = c(SELF, "Chen 2024", "Hazelwood 2025"),
  x0  = CX + c(0, -0.866, 0.866) * D,
  y0  = CY + c(1, -0.5, -0.5) * D,
  col = c(BLUE, CYAN, GREY))

gx <- seq(CX - (D + R) * 1.02, CX + (D + R) * 1.02, length.out = 700)
gr <- expand.grid(x = gx, y = seq(CY - (D + R) * 1.02, CY + (D + R) * 1.02, length.out = 700))
ins <- function(i) (gr$x - off$x0[i])^2 + (gr$y - off$y0[i])^2 <= R^2
key <- paste(ins(1), ins(2), ins(3))
cent <- function(k) { i <- which(key == k); c(x = mean(gr$x[i]), y = mean(gr$y[i])) }
K <- c(v5_only = "TRUE FALSE FALSE",   chen_only = "FALSE TRUE FALSE",
       haz_only = "FALSE FALSE TRUE",  v5_chen = "TRUE TRUE FALSE",
       v5_haz = "TRUE FALSE TRUE",     chen_haz = "FALSE TRUE TRUE",
       triple = "TRUE TRUE TRUE")
vals <- c(v5_only = cnt(FALSE, TRUE, FALSE), chen_only = cnt(TRUE, FALSE, FALSE),
          haz_only = cnt(FALSE, FALSE, TRUE), v5_chen = cnt(TRUE, TRUE, FALSE),
          v5_haz = cnt(FALSE, TRUE, TRUE),   chen_haz = cnt(TRUE, FALSE, TRUE),
          triple = cnt(TRUE, TRUE, TRUE))
cat("   regions:", paste(sprintf("%s=%d", names(vals), vals), collapse = " "), "\n")
stopifnot(sum(vals) == length(allg), vals[["v5_only"]] == 34)

lab <- do.call(rbind, lapply(names(K), function(nm) {
  cc <- cent(K[[nm]]); data.frame(x = cc[["x"]], y = cc[["y"]], txt = unname(vals[nm]), nm = nm)
}))
lab$col <- ifelse(lab$nm %in% c("v5_only", "v5_chen", "v5_haz", "triple"), BLUE, COL$ink2)

nms <- data.frame(x = c(CX, 1.22, 4.28), y = c(4.99, 0.20, 0.20),
                  txt = c(SELF, "Chen 2024", "Hazelwood 2025"),
                  col = c(BLUE, CYAN, COL$ink2))
szdf <- data.frame(set = c("Chen 2024", SELF, "Hazelwood 2025"),
                   n = c(n_chen, n_v5, n_haz), col = c(CYAN, BLUE, GREY),
                   y = c(3.92, 2.55, 1.18))
BMAX <- 2.95; BX <- 6.90
szdf$len <- szdf$n / max(szdf$n) * BMAX

p5b <- ggplot() +
  geom_circle(data = off, aes(x0 = x0, y0 = y0, r = R, fill = col, colour = col),
              alpha = 0.20, linewidth = 0.55) +
  geom_text(data = lab, aes(x = x, y = y, label = txt, colour = col), size = sz(7.2)) +
  geom_text(data = nms, aes(x = x, y = y, label = txt, colour = col),
            size = sz(7.2), fontface = "bold") +
  geom_text(aes(x = BX - 0.16, y = szdf$y, label = szdf$set), hjust = 1,
            size = sz(6.6), colour = COL$ink) +
  geom_rect(aes(xmin = BX, xmax = BX + szdf$len, ymin = szdf$y - 0.25,
                ymax = szdf$y + 0.25), fill = szdf$col, colour = NA) +
  geom_text(aes(x = BX + szdf$len + 0.16, y = szdf$y, label = szdf$n), hjust = 0,
            size = sz(6.6), colour = COL$ink, fontface = "bold") +
  annotate("text", x = BX, y = 4.72, label = "Set size", hjust = 0,
           size = sz(6.6), fontface = "bold", colour = COL$ink) +
  scale_colour_identity() + scale_fill_identity() +
  coord_fixed(ratio = 1, xlim = c(0.20, 11.00), ylim = c(0.00, 5.22), expand = FALSE) +
  theme_void(base_family = "Liberation Sans") +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(3, 4, 3, 4))
save_panel(p5b, file.path(outdir, "Fig5b_gene_overlap_venn.pdf"),
           mm_snap(169.69), mm_snap(82.20))

cat("\n== rev8 Figure 5 done ==\n")
