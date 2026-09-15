#!/usr/bin/env Rscript
# 32_tier1_smr_or.R — 6 个非-MHC Tier 1 基因的 per-SD SMR 效应量转 OR + 95%CI 落盘
# 回应逻辑断点 L5：SMR b 的临床含义（OR = exp(b_SMR)；95%CI = exp(b_SMR ± 1.96·se_SMR)）
# 铁律 #14：OR/CI 全部由脚本从 phase1 SMR 输出计算，禁止手敲。

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

proj <- "/ifs1/User/zhouman/project9-v5-crc-atlas"

## — 读入数据 ——————————————————————————————————————
combined <- read_csv(file.path(proj, "results", "phase3_susie_combined.csv"),
                     show_col_types = FALSE)   # 19 Tier 1 基因
extended <- read_csv(file.path(proj, "results", "phase3c_susie_extended.csv"),
                     show_col_types = FALSE)   # 13 扩展基因

# 6 个非-MHC Tier 1 初始基因 = combined - extended（按 ensg anti_join）
# 注意：不能用 is_mhc==FALSE 过滤，因为 13 个扩展里也有非-MHC Tier 1 基因（会得到 13 个）
initial <- anti_join(combined, extended, by = "ensg")
cat(sprintf("Initial non-MHC Tier 1 genes (combined - extended): %d\n", nrow(initial)))

# SMR 效应量源
smr <- read_tsv(file.path(proj, "results", "phase1_bonferroni_significant_annotated.tsv"),
                show_col_types = FALSE)
smr <- smr %>% distinct(probeID, .keep_all = TRUE)

## — 连接 + 计算 OR/CI ————————————————————————————
df <- initial %>%
  select(gene, ensg, PPH4, n_cs) %>%
  left_join(
    smr %>% select(probeID, SYMBOL, b_SMR, se_SMR, p_SMR),
    by = c("ensg" = "probeID")
  ) %>%
  mutate(
    OR        = exp(b_SMR),
    ci_low    = exp(b_SMR - 1.96 * se_SMR),
    ci_high   = exp(b_SMR + 1.96 * se_SMR),
    direction = if_else(b_SMR > 0, "Risk", "Protective")
  ) %>%
  arrange(desc(OR))

## — 后置断言（铁律 #5 / #14） ————————————————————
stopifnot(nrow(df) == 6)

# 期望值（来自 CLAUDE.md §P0-1；用于交叉校验脚本计算，不手敲进产出表）
expected <- tibble(
  gene    = c("BMP2", "TMBIM1", "UTP11", "GPBAR1", "LAMC1", "PNKD"),
  OR_r2   = c(1.38, 1.35, 1.27, 0.78, 1.14, 1.12),
  lo_r2   = c(1.27, 1.25, 1.20, 0.73, 1.11, 1.09),
  hi_r2   = c(1.50, 1.46, 1.36, 0.83, 1.17, 1.15)
)
chk <- df %>%
  left_join(expected, by = "gene") %>%
  mutate(ok = round(OR, 2) == OR_r2 &
               round(ci_low, 2) == lo_r2 &
               round(ci_high, 2) == hi_r2)
stopifnot(all(chk$ok), all(df$gene %in% expected$gene))

cat(sprintf("\n%-8s %-12s %-9s %-9s %-9s %-10s %s\n",
            "Gene", "b_SMR", "OR", "CI_low", "CI_high", "Direction", "p_SMR"))
for (i in seq_len(nrow(df))) {
  r <- df[i, ]
  cat(sprintf("%-8s %-12.6f %-9.2f %-9.2f %-9.2f %-10s %.3e\n",
              r$gene, r$b_SMR, r$OR, r$ci_low, r$ci_high, r$direction, r$p_SMR))
}

## — 落盘 ————————————————————————————————————————
out <- file.path(proj, "results", "phase_tier1_smr_or.csv")
write_csv(df, out)
cat(sprintf("\n✅ Written %d rows -> %s\n", nrow(df), out))
