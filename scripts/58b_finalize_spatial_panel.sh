#!/bin/bash
# 58b_finalize_spatial_panel.sh -- reproduce the tier-1 spatial panel end to end.
#
#   58_spatial_tier1_panel.R          -> results/figures/FigS8_spatial_tier1_genes_raw.pdf  (28 MB)
#   79_compress_figure_pdf.py         -> results/figures/FigS8_spatial_tier1_genes.pdf      (1.5 MB)
#
# The second file is what rev3_main_figures.R renames to FigS15_spatial_tier1_genes.pdf.
# Before 2026-09-29 the compression step was an undocumented manual step, so the
# shipped FigS15 could not be regenerated from the repository; it is scripted now.
set -e
cd "$(dirname "$0")/.."
Rscript scripts/58_spatial_tier1_panel.R
python3 scripts/79_compress_figure_pdf.py \
    results/figures/FigS8_spatial_tier1_genes_raw.pdf \
    results/figures/FigS8_spatial_tier1_genes.pdf
