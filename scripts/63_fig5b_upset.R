#!/usr/bin/env Rscript
# Fig 5b - shared / unique genes: this study vs Chen 2024 vs Hazelwood 2025.
# 3 sets -> Venn diagram. Region counts sit at the computed centroid of each region.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(ggforce)
})
setwd("/ifs1/User/zhouman/project9-v5-crc-atlas")
OUT <- "results/figures_rev3/Fig5b_gene_overlap_venn.pdf"

## ---------- data ----------
v5 <- unique(read_tsv("results/phase1_bonferroni_significant_annotated.tsv",
                      show_col_types = FALSE)$SYMBOL)
v5 <- sort(unique(v5[!is.na(v5) & v5 != ""]))

t8 <- readLines("manuscript/tables/TableS8_chen2024_gene_list.tex")
keep <- grepl("^\\s*[A-Za-z0-9][A-Za-z0-9._-]*\\s*&", t8)
chen <- sub("^\\s*([A-Za-z0-9._-]+)\\s*&.*$", "\\1", t8[keep])
chen <- sort(unique(chen[chen != "Gene"]))

haz <- sort(unique(c("AAMP","ABCC2","AC087392.1","AC144831.1","AL121832.3","ARPC2","ATF1",
  "CCM2","CENPBD2P","COLCA1","COX14","COX15","DACT1","EPM2AIP1","FADS1","FEN1","GPATCH1",
  "KLF5","LAMC1","LAMC1-AS1","LIMA1","LRRFIP2","METRNL","MLH1","MZF1","MZF1-AS1","PLEKHG6",
  "PNKD","POU2AF2","POU2AF3","POU5F1B","RBBP8NL","RP11-129K12.1","RPS21-DT","SEMA4D","TCF19",
  "TMBIM1","GPBAR1","LTBR","PDCD1","PTGER3")))

allg <- unique(c(v5, chen, haz))
mem <- data.frame(gene = allg, Chen = allg %in% chen, v5 = allg %in% v5, Haz = allg %in% haz)
cnt <- function(a, b, c) sum(mem$Chen == a & mem$v5 == b & mem$Haz == c)
n_v5 <- length(v5); n_chen <- length(chen); n_haz <- length(haz)

BLUE <- "#2166AC"; ORANGE <- "#E08214"; GREY <- "#6E6E6E"
SELF <- "This study"

## ---------- circle geometry ----------
CX <- 2.75; CY <- 2.40; R <- 1.28; D <- 1.15
off <- data.frame(
  set = c(SELF, "Chen 2024", "Hazelwood 2025"),
  x0  = CX + c(0, -0.866, 0.866) * D,
  y0  = CY + c(1, -0.5, -0.5) * D,
  col = c(BLUE, ORANGE, GREY))

## ---------- region centroids (numeric) ----------
gx <- seq(CX - (D + R) * 1.02, CX + (D + R) * 1.02, length.out = 600)
gr <- expand.grid(x = gx, y = seq(CY - (D + R) * 1.02, CY + (D + R) * 1.02, length.out = 600))
ins <- function(i) (gr$x - off$x0[i])^2 + (gr$y - off$y0[i])^2 <= R^2
key <- paste(ins(1), ins(2), ins(3))
cent <- function(k) { i <- which(key == k); c(x = mean(gr$x[i]), y = mean(gr$y[i])) }
K <- c(v5_only = "TRUE FALSE FALSE",   chen_only = "FALSE TRUE FALSE",
       haz_only = "FALSE FALSE TRUE",  v5_chen = "TRUE TRUE FALSE",
       v5_haz = "TRUE FALSE TRUE",     chen_haz = "FALSE TRUE TRUE",
       triple = "TRUE TRUE TRUE")
vals <- c(v5_only = cnt(FALSE,TRUE,FALSE), chen_only = cnt(TRUE,FALSE,FALSE),
          haz_only = cnt(FALSE,FALSE,TRUE), v5_chen = cnt(TRUE,TRUE,FALSE),
          v5_haz = cnt(FALSE,TRUE,TRUE),   chen_haz = cnt(TRUE,FALSE,TRUE),
          triple = cnt(TRUE,TRUE,TRUE))
print(vals)

lab <- do.call(rbind, lapply(names(K), function(nm) {
  cc <- cent(K[[nm]]); data.frame(x = cc[["x"]], y = cc[["y"]], txt = unname(vals[nm]), nm = nm)
}))
lab$col <- ifelse(lab$nm %in% c("v5_only", "v5_chen", "v5_haz", "triple"), BLUE, "grey20")

nms <- data.frame(x = c(CX, 1.25, 4.25), y = c(4.99, 0.26, 0.26),
                  txt = c(SELF, "Chen 2024", "Hazelwood 2025"),
                  col = c(BLUE, ORANGE, GREY))

szdf <- data.frame(set = c("Chen 2024", SELF, "Hazelwood 2025"),
                   n = c(n_chen, n_v5, n_haz), col = c(ORANGE, BLUE, GREY),
                   y = c(3.90, 2.55, 1.20))
BMAX <- 2.95; BX <- 6.85
szdf$len <- szdf$n / max(szdf$n) * BMAX

th <- theme_void(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(3, 4, 3, 4))

p <- ggplot() +
  geom_circle(data = off, aes(x0 = x0, y0 = y0, r = R, fill = col, colour = col),
              alpha = 0.20, linewidth = 0.55) +
  geom_text(data = lab, aes(x = x, y = y, label = txt, colour = col), size = 2.9) +
  geom_text(data = nms, aes(x = x, y = y, label = txt, colour = col),
            size = 2.9, fontface = "bold") +
  geom_text(aes(x = BX - 0.14, y = szdf$y, label = szdf$set),
            hjust = 1, size = 2.7, colour = "grey15") +
  geom_rect(aes(xmin = BX, xmax = BX + szdf$len, ymin = szdf$y - 0.26, ymax = szdf$y + 0.26),
            fill = szdf$col, colour = NA) +
  geom_text(aes(x = BX + szdf$len + 0.14, y = szdf$y, label = szdf$n),
            hjust = 0, size = 2.7, colour = "grey15") +
  annotate("text", x = BX, y = 4.72, label = "Set size", hjust = 0,
           size = 2.7, colour = "grey35") +
  scale_colour_identity() + scale_fill_identity() +
  coord_fixed(ratio = 1, xlim = c(0.25, 10.90), ylim = c(0.02, 5.20), expand = FALSE) +
  th

ggsave(OUT, p, width = 6.30, height = 3.06, device = cairo_pdf)
cat(sprintf("  -> %s\n", OUT))