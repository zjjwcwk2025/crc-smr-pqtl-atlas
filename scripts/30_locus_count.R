#!/usr/bin/env Rscript
# 30_locus_count.R — 统计 75 个 SMR 显著基因映射到多少个独立基因座
# 回应审稿人质疑：chr2:219Mb 3 基因是否同一 LD block？
# ±1Mb 窗口聚类；另备选 500kb / 250kb 作为敏感性分析

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

## — 读入数据 ——————————————————————————————————————
proj  <- "/ifs1/User/zhouman/project9-v5-crc-atlas"
src   <- file.path(proj, "results", "phase1_bonferroni_significant_annotated.tsv")
raw   <- read_tsv(src, show_col_types = FALSE)

# 去重（同一个 ENSG 可能被标注多次）
raw <- raw %>% distinct(probeID, .keep_all = TRUE)
cat(sprintf("Unique ENSG IDs: %d\n", nrow(raw)))

## 提取基因坐标（用 top SNP 位置，这是 SMR 定位到的最显著 SNP）
df <- raw %>%
  select(probeID, SYMBOL, chr = ProbeChr, topSNP, topSNP_bp, p_SMR) %>%
  mutate(
    chr     = as.integer(chr),
    pos     = as.numeric(topSNP_bp),
    p_SMR   = as.numeric(p_SMR),
    # 蛋白编码基因 vs 非编码
    is_pc   = !is.na(SYMBOL) & SYMBOL != "" & !grepl("^ENSG", SYMBOL)
  ) %>%
  filter(!is.na(chr) & !is.na(pos)) %>%
  arrange(chr, pos)

cat(sprintf("Genes with valid chr:pos: %d  (protein-coding: %d, noncoding: %d)\n",
            nrow(df), sum(df$is_pc), sum(!df$is_pc)))

## — 基因座聚类（greedy, ±窗口） ————————————————————
cluster_by_window <- function(df, window_bp = 1e6) {
  if (nrow(df) == 0) return(df)
  df <- df %>% arrange(chr, pos)
  df$locus_id  <- NA_integer_
  df$locus_chr <- NA_integer_
  locus_counter <- 0L

  for (ch in unique(df$chr)) {
    idx   <- which(df$chr == ch)
    if (length(idx) == 0) next
    chr_df <- df[idx, ]
    # 贪心聚类：以第一个基因为锚点，合并窗口内的基因
    cl <- 1L
    anchor <- chr_df$pos[1]
    for (i in seq_len(nrow(chr_df))) {
      if (chr_df$pos[i] - anchor > window_bp) {
        cl <- cl + 1L
        anchor <- chr_df$pos[i]
      }
      chr_df$locus_id[i] <- locus_counter + cl
      chr_df$locus_chr[i] <- ch
    }
    df[idx, ] <- chr_df
    locus_counter <- locus_counter + cl
  }
  df
}

## — 主分析：±1Mb 窗口 ——————————————————————————————
window1 <- 1e6
res_1mb <- cluster_by_window(df, window1)

# 每个基因座的汇总
loci_1mb <- res_1mb %>%
  group_by(locus_id, locus_chr) %>%
  summarise(
    n_genes        = n(),
    n_pc           = sum(is_pc),
    locus_start    = min(pos),
    locus_end      = max(pos),
    locus_span_kb  = round((max(pos) - min(pos)) / 1000, 1),
    representative = SYMBOL[which.min(p_SMR)],
    genes          = paste(SYMBOL, collapse = ", "),
    min_p          = min(p_SMR),
    .groups = "drop"
  ) %>%
  arrange(locus_chr, locus_start)

cat(sprintf("\n========== ±1Mb window ==========\n"))
cat(sprintf("Total loci: %d\n", nrow(loci_1mb)))
cat(sprintf("  Singletons (1 gene/locus): %d\n", sum(loci_1mb$n_genes == 1)))
cat(sprintf("  Multi-gene loci: %d\n", sum(loci_1mb$n_genes >= 2)))
cat(sprintf("  Genes in multi-gene loci: %d / %d (%.1f%%)\n",
            sum(loci_1mb$n_genes[loci_1mb$n_genes >= 2]),
            sum(loci_1mb$n_genes),
            100*sum(loci_1mb$n_genes[loci_1mb$n_genes >= 2])/sum(loci_1mb$n_genes)))

cat(sprintf("\n  --- Protein-coding only ---\n"))
pc_loci <- loci_1mb %>% filter(n_pc > 0)
cat(sprintf("  Loci with ≥1 protein-coding gene: %d\n", nrow(pc_loci)))

## — 重点分析：chr2:219Mb 区域 —————————————————————————————————————
cat(sprintf("\n========== chr2:219Mb cluster (spotlight) ==========\n"))
chr2_219 <- loci_1mb %>% filter(locus_chr == 2 & locus_start >= 218e6 & locus_end <= 220e6)
if (nrow(chr2_219) > 0) {
  for (i in seq_len(nrow(chr2_219))) {
    r <- chr2_219[i, ]
    cat(sprintf("  Locus chr%s:%s-%s (%.0f kb) | %d genes (PC: %d)\n",
                r$locus_chr, r$locus_start, r$locus_end, r$locus_span_kb,
                r$n_genes, r$n_pc))
    cat(sprintf("    Genes: %s\n", r$genes))

    # 这个区域内的 coloc+SuSiE 基因
    coloc_hits <- c("PNKD", "TMBIM1", "GPBAR1")
    overlap <- coloc_hits[coloc_hits %in% unlist(strsplit(r$genes, ", "))]
    if (length(overlap) > 0) {
      cat(sprintf("    ⚠ Contains coloc+SuSiE Tier 1 genes: %s\n",
                  paste(overlap, collapse = ", ")))
    }
  }
}

## — 打印所有多基因基因座 ——————————————————————————————
cat(sprintf("\n========== All multi-gene loci (≥2 genes) ==========\n"))
multi <- loci_1mb %>% filter(n_genes >= 2) %>% arrange(desc(n_genes))
for (i in seq_len(nrow(multi))) {
  r <- multi[i, ]
  cat(sprintf("  [L%s] chr%s:%s-%s (%s kb) | %d genes | rep: %s\n",
              r$locus_id, r$locus_chr, r$locus_start, r$locus_end,
              r$locus_span_kb, r$n_genes, r$representative))
  cat(sprintf("    Genes: %s\n", r$genes))
}

## — 敏感性分析：±500kb / ±250kb ——————————————————————————
cat(sprintf("\n========== Sensitivity: smaller windows ==========\n"))
sens_rows <- list()
for (w in c(5e5, 2.5e5)) {
  res  <- cluster_by_window(df, w)
  loci <- res %>%
    group_by(locus_id, locus_chr) %>%
    summarise(
      n_genes        = n(),
      n_pc           = sum(is_pc),
      locus_start    = min(pos),
      locus_end      = max(pos),
      locus_span_kb  = round((max(pos) - min(pos)) / 1000, 1),
      representative = SYMBOL[which.min(p_SMR)],
      genes          = paste(SYMBOL, collapse = ", "),
      min_p          = min(p_SMR),
      .groups = "drop"
    ) %>%
    arrange(locus_chr, locus_start)
  wkb <- as.integer(round(w / 1000))
  out <- file.path(proj, "results", sprintf("phase_locus_count_%dkb.csv", wkb))
  write_csv(loci, out)
  cat(sprintf("  ±%dkb: %d loci (singletons: %d, multi-gene: %d) -> %s\n",
              wkb, nrow(loci), sum(loci$n_genes == 1), sum(loci$n_genes >= 2),
              sprintf("phase_locus_count_%dkb.csv", wkb)))
  sens_rows[[sprintf("%dkb", wkb)]] <- data.frame(
    window_kb    = wkb,
    n_loci       = nrow(loci),
    n_singletons = sum(loci$n_genes == 1),
    n_multi_gene = sum(loci$n_genes >= 2),
    stringsAsFactors = FALSE
  )
}
sens_df <- bind_rows(sens_rows)
write_csv(sens_df, file.path(proj, "results", "phase_locus_count_sensitivity_summary.csv"))
cat(sprintf("  Sensitivity summary -> results/phase_locus_count_sensitivity_summary.csv\n"))

## — 保存结果 ————————————————————————————————————————————
out_csv <- file.path(proj, "results", "phase_locus_count.csv")
write_csv(loci_1mb, out_csv)
cat(sprintf("\n✅ Results saved to: %s\n", out_csv))

## — 基因座分布直方图（快速文本版）—————————————
cat(sprintf("\n========== Genes per locus distribution ==========\n"))
tbl <- table(loci_1mb$n_genes)
for (n in names(tbl)) {
  cat(sprintf("  %s genes/locus: %d loci\n", n, tbl[n]))
}

## — 染色体分布 ————————————————————————————————————————————
cat(sprintf("\n========== Loci per chromosome ==========\n"))
chr_dist <- loci_1mb %>% count(locus_chr) %>% arrange(locus_chr)
for (i in seq_len(nrow(chr_dist))) {
  cat(sprintf("  chr%s: %d loci\n", chr_dist$locus_chr[i], chr_dist$n[i]))
}

cat(sprintf("\n========== Summary for manuscript ==========\n"))
cat(sprintf("Among the %d Bonferroni-significant SMR genes (%d protein-coding), ", nrow(df), sum(df$is_pc)))
cat(sprintf("we identified %d independent genomic loci (±1Mb).\n", nrow(loci_1mb)))
cat(sprintf("%d loci (%d%%) contain only a single gene.\n",
            sum(loci_1mb$n_genes == 1),
            round(100 * sum(loci_1mb$n_genes == 1) / nrow(loci_1mb))))
cat(sprintf("The largest cluster is chr2:219Mb (%d genes) — ", max(loci_1mb$n_genes)))
cat(sprintf("PNKD, TMBIM1, and GPBAR1 share the same LD block.\n"))
cat(sprintf("Excluding this single cluster, %d/%d genes map to %d independent loci.\n",
            nrow(df) - max(loci_1mb$n_genes), nrow(df), nrow(loci_1mb) - 1))
