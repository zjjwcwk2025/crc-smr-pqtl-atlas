#!/usr/bin/env Rscript
#'
#' Phase 5c-6: Pathway Enrichment for 75 SMR genes
#' GO + KEGG enrichment analysis using clusterProfiler
#'
suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

PROJ <- "/ifs1/User/zhouman/project9-v5-crc-atlas"
setwd(PROJ)
dir.create("results/phase5c_pathway", showWarnings = FALSE, recursive = TRUE)

# ── 1. Load 75 SMR genes ──
cat("[1/4] Loading SMR Bonferroni-significant genes...\n")
smr <- fread("results/phase1_bonferroni_significant_annotated.tsv")

genes <- unique(smr$SYMBOL[smr$SYMBOL != "" & !is.na(smr$SYMBOL)])
cat(sprintf("  %d unique gene symbols\n", length(genes)))

# Convert to Entrez IDs
entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
cat(sprintf("  %d genes mapped to Entrez IDs\n", nrow(entrez)))

# ── 2. GO Enrichment ──
cat("\n[2/4] Running GO enrichment...\n")

go_bp <- enrichGO(gene = entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

go_mf <- enrichGO(gene = entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

go_cc <- enrichGO(gene = entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# ── 3. KEGG enrichment ──
cat("[3/4] Running KEGG enrichment...\n")
kegg <- enrichKEGG(gene = entrez$ENTREZID,
                   organism = "hsa",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)
if (!is.null(kegg) && nrow(kegg) > 0) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# ── 4. Generate plots ──
cat("[4/4] Generating plots...\n")

plot_and_save <- function(enrich_result, name, title) {
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    cat(sprintf("  %s: no significant terms\n", name))
    return()
  }

  # Top 15 terms
  top_terms <- enrich_result@result %>%
    arrange(p.adjust) %>%
    head(15) %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  # Dotplot
  p_dot <- ggplot(top_terms, aes(x = Count, y = Description, color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "FDR") +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = title, x = "Gene Count", y = "") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))

  ggsave(sprintf("results/phase5c_pathway/%s_dotplot.pdf", name),
         p_dot, width = 10, height = 5)

  # Barplot
  p_bar <- ggplot(top_terms, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "red", high = "blue", name = "FDR") +
    labs(title = title, x = "Gene Count", y = "") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))

  ggsave(sprintf("results/phase5c_pathway/%s_barplot.pdf", name),
         p_bar, width = 10, height = 5)

  cat(sprintf("  %s: %d significant terms\n", name, nrow(enrich_result)))

  # Save full table
  fwrite(enrich_result@result,
         sprintf("results/phase5c_pathway/%s_table.csv", name))
}

plot_and_save(go_bp, "GO_BP", "GO Biological Process (75 SMR genes)")
plot_and_save(go_mf, "GO_MF", "GO Molecular Function (75 SMR genes)")
plot_and_save(go_cc, "GO_CC", "GO Cellular Component (75 SMR genes)")
plot_and_save(kegg, "KEGG", "KEGG Pathways (75 SMR genes)")

# ── Summary Table ──
cat("\n── Pathway Enrichment Summary ──\n")

summary_rows <- list()
for (name in c("GO_BP", "GO_MF", "GO_CC", "KEGG")) {
  table_path <- sprintf("results/phase5c_pathway/%s_table.csv", name)
  if (file.exists(table_path)) {
    tbl <- fread(table_path)
    if (nrow(tbl) > 0) {
      cat(sprintf("\n%s: Top 5 terms\n", name))
      for (j in seq_len(min(5, nrow(tbl)))) {
        cat(sprintf("  %d. %s (Count=%d, FDR=%.2e)\n",
                    j, tbl$Description[j], tbl$Count[j], tbl$p.adjust[j]))
      }
    }
  }
}

cat("\nDone. Results saved to results/phase5c_pathway/\n")
