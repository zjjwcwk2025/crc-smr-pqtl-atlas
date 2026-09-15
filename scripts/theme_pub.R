# Unified publication theme for project9 (main and supplementary figures).
# Figures are generated at final print size (single column = 6.30 in) so that no
# text falls below 7 pt after inclusion in the manuscript.
suppressMessages({library(ggplot2)})

FONT <- "Liberation Sans"

PAL <- list(
  coloc     = c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#CA0020"),
  class     = c("A: pQTL concordant"="#2166AC","B: pQTL discordant"="#B2182B","C: pQTL not significant"="#999999"),
  direction = c("concordant"="#2166AC","discordant"="#B2182B","not_significant"="#999999"),
  safety    = c("High alert"="#B2182B","Moderate"="#F4A582","Low"="#92C5DE"),
  accent    = c("#2166AC","#B2182B","#1B7837","#F4A582","#999999"),
  tint      = c("#D6E4F2","#E6EEF7","#EDEDED")
)

theme_pub <- function(base_size = 9, base_family = FONT) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text               = element_text(colour = "grey10"),
      axis.title         = element_text(size = base_size, face = "bold", colour = "grey10"),
      axis.text          = element_text(size = base_size - 1, colour = "black"),
      axis.line          = element_line(linewidth = 0.35, colour = "grey20"),
      axis.ticks         = element_line(linewidth = 0.35, colour = "grey20"),
      axis.ticks.length  = unit(2, "pt"),
      plot.title         = element_blank(),
      plot.subtitle      = element_blank(),
      plot.caption       = element_blank(),
      legend.title       = element_text(size = base_size - 1, face = "bold"),
      legend.text        = element_text(size = base_size - 1.5),
      legend.key.size    = unit(9, "pt"),
      legend.position    = "bottom",
      legend.background  = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(size = base_size - 1, face = "bold"),
      plot.margin        = margin(5, 6, 3, 4)
    )
}

save_pub <- function(plot, file, width, height, ...) {
  ggsave(file, plot, width = width, height = height, device = cairo_pdf, bg = "white", ...)
  cat(sprintf("  -> %s (%.2f x %.2f in)\n", basename(file), width, height))
}
