#!/usr/bin/env Rscript
#'
#' Phase 5c-8: Protein-Protein Interaction Network for 6 Tier 1 genes
#' Uses STRING REST API v12.0
#'
suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(ggplot2)
  library(ggrepel)
  library(data.table)
  library(igraph)
})

PROJ <- "/ifs1/User/zhouman/project9-v5-crc-atlas"
setwd(PROJ)
dir.create("results/phase5c_ppi", showWarnings = FALSE, recursive = TRUE)

# ── 1. Load Tier 1 genes ──
susie <- fread("results/phase3_susie_finemap.csv")
tier1_genes <- susie[tier == 1, gene]
cat("Tier 1 genes:", paste(tier1_genes, collapse = ", "), "\n")

# ── 2. Query STRING API ──
cat("\n[1/3] Querying STRING API for PPI interactions...\n")

STRING_BASE <- "https://string-db.org/api/json"
gene_list <- paste(tier1_genes, collapse = "%0d")

# Get STRING identifiers
res <- GET(sprintf("%s/get_string_ids?identifiers=%s&species=9606&caller_identity=project9_crc",
                   STRING_BASE, gene_list))
if (status_code(res) != 200) {
  cat(sprintf("  ERROR: STATUS %d\n", status_code(res)))
  cat("  Response:", content(res, "text"), "\n")
  quit(status = 1)
}

string_ids <- fromJSON(content(res, "text", encoding = "UTF-8"))
cat(sprintf("  %d genes mapped to STRING IDs\n", nrow(string_ids)))

# Get PPI network
string_id_list <- paste(string_ids$stringId, collapse = "%0d")
res2 <- GET(sprintf("%s/network?identifiers=%s&species=9606&required_score=400&add_nodes=10&caller_identity=project9_crc",
                    STRING_BASE, string_id_list))

if (status_code(res2) != 200) {
  cat(sprintf("  ERROR: network API returned status %d\n", status_code(res2)))
  quit(status = 1)
}

ppi <- as.data.table(fromJSON(content(res2, "text", encoding = "UTF-8")))
cat(sprintf("  %d interactions, %d nodes\n", nrow(ppi),
            length(unique(c(ppi$preferredName_A, ppi$preferredName_B)))))
fwrite(ppi, "results/phase5c_ppi/string_interactions.csv")

# ── 3. Build network ──
cat("\n[2/3] Building network visualization...\n")

edges <- ppi[, .(from = preferredName_A, to = preferredName_B, score)]
node_names <- unique(c(edges$from, edges$to))

g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = node_names))

# Layout
set.seed(42)
coords <- layout_with_fr(g)

# Build node/edge data for ggplot
nodes_df <- data.frame(
  name = V(g)$name,
  x = coords[, 1],
  y = coords[, 2],
  is_tier1 = V(g)$name %in% tier1_genes,
  degree = degree(g),
  stringsAsFactors = FALSE
)

el <- as_edgelist(g, names = TRUE)
edges_df <- data.frame(
  from = el[, 1],
  to = el[, 2],
  score = E(g)$score,
  stringsAsFactors = FALSE
)
edges_df$x_from <- nodes_df$x[match(edges_df$from, nodes_df$name)]
edges_df$y_from <- nodes_df$y[match(edges_df$from, nodes_df$name)]
edges_df$x_to <- nodes_df$x[match(edges_df$to, nodes_df$name)]
edges_df$y_to <- nodes_df$y[match(edges_df$to, nodes_df$name)]

# Plot
p <- ggplot() +
  geom_segment(data = edges_df,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                   linewidth = score, alpha = score),
               color = "grey60") +
  scale_linewidth_continuous(range = c(0.2, 2), guide = "none") +
  scale_alpha_continuous(range = c(0.15, 0.7), guide = "none") +
  geom_point(data = nodes_df,
             aes(x = x, y = y, color = is_tier1, size = degree)) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                     labels = c("TRUE" = "Tier 1 Gene", "FALSE" = "Interaction Partner"),
                     name = "") +
  scale_size_continuous(range = c(2, 10), name = "Degree") +
  geom_text_repel(data = nodes_df[nodes_df$name %in% tier1_genes, ],
                  aes(x = x, y = y, label = name),
                  size = 4.5, fontface = "bold", color = "#E41A1C",
                  box.padding = 0.6, max.overlaps = 20) +
  geom_text_repel(data = nodes_df[!(nodes_df$name %in% tier1_genes) & nodes_df$degree >= 2, ],
                  aes(x = x, y = y, label = name),
                  size = 3, color = "grey30", box.padding = 0.3, max.overlaps = 30) +
  labs(title = "STRING PPI Network: Tier 1 CRC Drug Target Genes",
       subtitle = sprintf("%d nodes, %d edges (confidence > 0.4)",
                          nrow(nodes_df), nrow(edges_df))) +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("results/phase5c_ppi/string_ppi_network.pdf", p, width = 12, height = 10)
ggsave("results/phase5c_ppi/string_ppi_network.png", p, width = 12, height = 10, dpi = 150)

# ── 4. Functional enrichment via STRING ──
cat("\n[3/3] Querying STRING functional enrichment...\n")

res3 <- GET(sprintf("%s/enrichment?identifiers=%s&species=9606&caller_identity=project9_crc",
                    STRING_BASE, string_id_list))
if (status_code(res3) == 200) {
  enrich <- as.data.table(fromJSON(content(res3, "text", encoding = "UTF-8")))
  if (nrow(enrich) > 0) {
    fwrite(enrich, "results/phase5c_ppi/string_enrichment.csv")
    cat(sprintf("  %d enriched terms retrieved\n", nrow(enrich)))

    enrich_sig <- enrich[fdr < 0.05]
    if (nrow(enrich_sig) > 0) {
      cats <- unique(enrich_sig$category)
      cat("\n-- STRING Functional Enrichment (FDR < 0.05) --\n")
      for (cn in cats) {
        ct <- enrich_sig[category == cn][order(fdr)]
        if (nrow(ct) > 0) {
          cat(sprintf("\n%s:\n", cn))
          for (j in seq_len(min(5, nrow(ct)))) {
            cat(sprintf("  %d. %s (FDR=%.2e, %d genes)\n",
                        j, ct$term[j], ct$fdr[j], ct$number_of_genes[j]))
          }
        }
      }

      # Plot
      top_terms <- enrich_sig[order(fdr)][1:20]
      top_terms[, Description := factor(term, levels = rev(term))]
      p_enrich <- ggplot(top_terms, aes(x = -log10(fdr), y = Description, fill = category)) +
        geom_bar(stat = "identity") +
        labs(x = expression(-log[10](FDR)), y = "",
             title = "STRING Enrichment: Tier 1 Gene PPI Network") +
        theme_bw(base_size = 9) +
        theme(axis.text.y = element_text(size = 7))
      ggsave("results/phase5c_ppi/string_enrichment.pdf", p_enrich, width = 12, height = 6)
      ggsave("results/phase5c_ppi/string_enrichment.png", p_enrich, width = 12, height = 6, dpi = 150)
    } else {
      cat("  No terms with FDR < 0.05\n")
    }
  }
} else {
  cat(sprintf("  STRING enrichment API returned status %d\n", status_code(res3)))
}

cat("\nDone. Results saved to results/phase5c_ppi/\n")
