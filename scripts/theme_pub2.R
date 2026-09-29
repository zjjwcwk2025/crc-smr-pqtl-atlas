# ============================================================
#  theme_pub2.R -- publication design system (rev8, 2026-09-28)
#  One theme file shared by every panel, main and supplementary.
#  Implements figure-spec.md section 8.4 (single shared theme) and
#  section 11 (colour rules: colour-blind safe, print-safe, one
#  meaning per colour, max 6-8 colours).
# ============================================================
suppressMessages({library(ggplot2); library(grDevices)})

FONT <- "Liberation Sans"          # metric clone of Arial/Helvetica

## geom_text size -> points: ggplot multiplies size by .pt (72.27/25.4)
.PT <- 25.4 / 72.27
sz  <- function(pt) pt * .PT

## --- semantic palette: one meaning per colour, fixed for the whole paper ---
COL <- list(
  ink       = "#1A1A1A",
  ink2      = "#4D4D4D",
  ink3      = "#6E6E6E",
  grid      = "#EAEAEA",
  axis      = "#5A5A5A",
  genetic   = "#1F5FA8",   # genetic layer / eQTL / blood / PASS
  genetic2  = "#9DC3E6",   # same layer, secondary (marginal PPH4)
  protein   = "#D1783C",   # protein layer / pQTL / plasma / GTEx colon
  alert     = "#B22222",   # discordant / high-risk  (always + 2nd encoding)
  mute      = "#A6A6A6",   # not significant / no measurement
  tint_blue = "#EAF1F9",
  tint_amb  = "#FBEFE6",
  tint_grey = "#F2F2F2"
)

ramp_blue <- grDevices::colorRampPalette(
  c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5", "#08306B"))

theme_pub2 <- function(base_pt = 7.4, base_family = FONT) {
  theme_classic(base_size = base_pt, base_family = base_family) +
    theme(
      text               = element_text(colour = COL$ink),
      axis.title         = element_text(size = base_pt, face = "bold", colour = COL$ink),
      axis.text          = element_text(size = base_pt * 0.90, colour = COL$ink2),
      axis.line          = element_line(linewidth = 0.40, colour = COL$axis),
      axis.ticks         = element_line(linewidth = 0.40, colour = COL$axis),
      axis.ticks.length  = unit(2.2, "pt"),
      plot.title         = element_blank(),
      plot.subtitle      = element_blank(),
      plot.caption       = element_blank(),
      legend.title       = element_text(size = base_pt * 0.95, face = "bold"),
      legend.text        = element_text(size = base_pt * 0.87),
      legend.key.size    = unit(8, "pt"),
      legend.position    = "bottom",
      legend.background  = element_blank(),
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      panel.grid.major.y = element_line(colour = COL$grid, linewidth = 0.30),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(size = base_pt * 0.95, face = "bold"),
      plot.margin        = margin(4, 6, 3, 4)
    )
}

## panel writer: always cairo_pdf (embeds fonts) at the final printed size
save_panel <- function(plot, file, w, h) {
  ggsave(file, plot, width = w, height = h, device = cairo_pdf, bg = "white")
  cat(sprintf("   -> %s  (%.2f x %.2f in)\n", basename(file), w, h))
}

## millimetres -> inches, so every panel is declared in journal units
mm <- function(x) x / 25.4

## millimetres -> inches, snapped to a whole number of points.
## cairo_pdf floors the page size to whole points: 7.43 in -> 534 pt, but
## 188.38 mm -> 533.99 pt -> 533 pt, one point short. Requesting the middle of
## the (n, n+1) point interval makes the saved page land on exactly
## round(mm/25.4*72) points, so a panel declared at 169.69 x 134.76 mm really
## is 481 x 382 pt.
mm_snap <- function(x) (round(x / 25.4 * 72) + 0.5) / 72