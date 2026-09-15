#!/usr/bin/env Rscript
# Phase 1 Step 3: Parse SMR results, Bonferroni filter, Manhattan & QQ plots
# Project 9 v5: CRC Multi-Omics Drug-Target Atlas

library(data.table)
library(ggplot2)

# --- Config ---
smr_file  <- "results/phase1_eqtl_smr.smr"
out_dir   <- "results"

# --- Read results ---
cat("[1/4] Reading SMR results...\n")
smr <- fread(smr_file, sep = "\t", header = TRUE)
cat(sprintf("  Total genes tested: %s\n", format(nrow(smr), big.mark=",")))

# Bonferroni threshold
n_tests <- nrow(smr)
bonf_threshold <- 0.05 / n_tests
cat(sprintf("  Bonferroni threshold: %.2e (0.05 / %s)\n", bonf_threshold, format(n_tests, big.mark=",")))

# Also note the standard genome-wide suggestive threshold
suggestive_threshold <- 1e-5

# --- Filter significant ---
bonf_sig <- smr[p_SMR < bonf_threshold]
suggestive <- smr[p_SMR >= bonf_threshold & p_SMR < suggestive_threshold]

cat(sprintf("\n[2/4] Results:\n"))
cat(sprintf("  Bonferroni significant (p < %.2e): %d genes\n", bonf_threshold, nrow(bonf_sig)))
cat(sprintf("  Suggestive (%.2e <= p < 1e-5): %d genes\n", bonf_threshold, nrow(suggestive)))
cat(sprintf("  Top 20 SMR p-values:\n"))

# Print top 20
top20 <- smr[order(p_SMR)][1:min(20, nrow(smr))]
for (i in 1:nrow(top20)) {
    cat(sprintf("    %2d. %s (Chr%s) p_SMR=%.2e b_SMR=%.4f topSNP=%s\n",
                i, top20$Gene[i], top20$ProbeChr[i],
                top20$p_SMR[i], top20$b_SMR[i], top20$topSNP[i]))
}

# --- Save significant results ---
if (nrow(bonf_sig) > 0) {
    bonf_sig <- bonf_sig[order(p_SMR)]
    fwrite(bonf_sig, file.path(out_dir, "phase1_bonferroni_significant.tsv"),
           sep = "\t", quote = FALSE)
    cat(sprintf("\n  Significant results saved to: phase1_bonferroni_significant.tsv\n"))
}

if (nrow(suggestive) > 0) {
    suggestive <- suggestive[order(p_SMR)]
    fwrite(suggestive, file.path(out_dir, "phase1_suggestive.tsv"),
           sep = "\t", quote = FALSE)
    cat(sprintf("  Suggestive results saved to: phase1_suggestive.tsv\n"))
}

# --- QQ Plot ---
cat("\n[3/4] Generating QQ plot...\n")

# Compute expected -log10(p)
smr_clean <- smr[!is.na(p_SMR) & p_SMR > 0 & p_SMR <= 1]
observed <- sort(smr_clean$p_SMR)
n <- length(observed)
expected <- (1:n) / (n + 1)

qq_data <- data.table(
    expected = -log10(expected),
    observed = -log10(observed)
)

# Genetic inflation factor lambda
chisq <- qchisq(1 - smr_clean$p_SMR, 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

pdf(file.path(out_dir, "phase1_qq_plot.pdf"), width = 6, height = 6)
par(mar = c(5, 5, 4, 2))
plot(qq_data$expected, qq_data$observed,
     pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3),
     xlab = expression(Expected ~ -log[10](italic(P))),
     ylab = expression(Observed ~ -log[10](italic(P))),
     main = "SMR eQTL → CRC (GCST90255675)",
     las = 1)
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c(
    sprintf("N genes = %s", format(nrow(smr_clean), big.mark = ",")),
    sprintf("λ = %.3f", lambda)
), bty = "n", cex = 0.9)
dev.off()
cat(sprintf("  QQ plot saved (λ = %.3f)\n", lambda))

# --- Manhattan Plot ---
cat("[4/4] Generating Manhattan plot...\n")

# Prepare manhattan data (use ProbeChr and Probe_bp)
manh <- copy(smr_clean)
manh[, CHR := as.integer(ProbeChr)]
manh <- manh[!is.na(CHR) & CHR >= 1 & CHR <= 22]

# Cumulative position
manh[, BP := Probe_bp]
setorder(manh, CHR, BP)

# Build cumulative axis
chr_info <- manh[, .(chr_len = max(BP)), by = CHR][order(CHR)]
chr_info[, cum_pos := cumsum(as.numeric(chr_len)) - chr_len]
manh <- merge(manh, chr_info[, .(CHR, cum_pos)], by = "CHR")
manh[, manh_pos := cum_pos + BP]

# Colors
manh[, color := ifelse(CHR %% 2 == 0, "#4DBBD5", "#E64B35")]

# Axis ticks at chromosome centers
chr_centers <- manh[, .(center = cum_pos[1] + max(BP) / 2), by = CHR]

pdf(file.path(out_dir, "phase1_manhattan.pdf"), width = 14, height = 5)
par(mar = c(4, 5, 3, 1))

# Only plot significant points for speed
plot(manh$manh_pos, -log10(manh$p_SMR),
     pch = 16, cex = 0.2,
     col = manh$color,
     xlab = "", ylab = expression(-log[10](italic(P)[SMR])),
     main = "SMR: eQTLGen → CRC (GCST90255675, 78K cases)",
     xaxt = "n", las = 1)
axis(1, at = chr_centers$center, labels = chr_centers$CHR, tick = FALSE, line = -0.5, cex.axis = 0.8)

# Bonferroni line
abline(h = -log10(bonf_threshold), col = "red", lwd = 2, lty = 2)
text(max(manh$manh_pos) * 0.98, -log10(bonf_threshold) + 0.3,
     sprintf("Bonferroni (%.1e)", bonf_threshold),
     col = "red", cex = 0.7, adj = 1)

# Suggestive line
abline(h = -log10(suggestive_threshold), col = "blue", lwd = 1, lty = 3)
text(max(manh$manh_pos) * 0.98, -log10(suggestive_threshold) + 0.2,
     "p = 1e-5", col = "blue", cex = 0.7, adj = 1)

# Build ENSG → gene symbol lookup from annotated file
gene_lookup <- NULL
annot_tsv <- file.path(out_dir, "phase1_bonferroni_significant_annotated.tsv")
if (file.exists(annot_tsv)) {
  annot <- fread(annot_tsv, sep = "\t", header = TRUE)
  if ("SYMBOL" %in% names(annot) && "Gene" %in% names(annot)) {
    gene_lookup <- setNames(annot$SYMBOL, annot$Gene)
    cat(sprintf("  Loaded %d gene symbol mappings from annotated TSV\n", nrow(annot)))
  }
}

# Label top hits with gene symbols (fall back to ENSG if no symbol)
if (nrow(bonf_sig) > 0) {
    top_hits <- bonf_sig[order(p_SMR)][1:min(10, nrow(bonf_sig))]
    for (i in 1:nrow(top_hits)) {
        chr <- top_hits$ProbeChr[i]
        bp <- top_hits$Probe_bp[i]
        pos <- chr_info[CHR == chr, cum_pos] + bp
        ensg_id <- top_hits$Gene[i]
        gene_label <- if (is.null(gene_lookup) || is.na(gene_lookup[ensg_id])) ensg_id else gene_lookup[ensg_id]
        text(pos, -log10(top_hits$p_SMR[i]),
             labels = gene_label,
             pos = 3, cex = 0.5, col = "black", offset = 0.3)
    }
}

dev.off()
cat("  Manhattan plot saved\n")

# --- Stop-loss decision ---
cat("\n========================================\n")
cat("PHASE 1 STOP-LOSS DECISION\n")
cat("========================================\n")
cat(sprintf("Bonferroni significant genes: %d\n", nrow(bonf_sig)))
cat(sprintf("Threshold: p < %.2e\n", bonf_threshold))

if (nrow(bonf_sig) < 3) {
    cat("\n🔴 STOP: < 3 Bonferroni significant genes. Archive and reconsider.\n")
} else if (nrow(bonf_sig) < 5) {
    cat(sprintf("\n⚠️  CAUTION: %d significant genes (3-4 range). Assess overlap with competitors before continuing.\n", nrow(bonf_sig)))
} else if (nrow(bonf_sig) < 10) {
    cat(sprintf("\n🟡 PASS with caution: %d significant genes. Continue to Phase 2.\n", nrow(bonf_sig)))
} else {
    cat(sprintf("\n🟢 PASS: %d significant genes. Strong signal. Continue to Phase 2.\n", nrow(bonf_sig)))
}

cat(sprintf("Suggestive genes (p < 1e-5): %d\n", nrow(suggestive)))
cat(sprintf("Total genes with p < 1e-5: %d\n", nrow(bonf_sig) + nrow(suggestive)))
cat("========================================\n")
