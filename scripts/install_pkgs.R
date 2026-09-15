#!/usr/bin/env Rscript
# Install missing R packages for enhancement analyses
options(repos = c(CRAN = "https://cloud.r-project.org"))

# remotes for GitHub installs
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# BiocManager for Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# locuszoomr from GitHub (regional association plots)
if (!requireNamespace("locuszoomr", quietly = TRUE)) {
  remotes::install_github("mrcieu/locuszoomr", upgrade = "never")
}

# STRINGdb from Bioconductor
if (!requireNamespace("STRINGdb", quietly = TRUE)) {
  BiocManager::install("STRINGdb", update = FALSE, ask = FALSE)
}

cat("Installation complete, checking packages:\n")
for (p in c("locuszoomr", "STRINGdb", "clusterProfiler", "org.Hs.eg.db",
            "TwoSampleMR", "MendelianRandomization", "MRPRESSO", "ggplot2", "dplyr")) {
  cat(sprintf("  %-25s: %s\n", p, requireNamespace(p, quietly = TRUE)))
}
