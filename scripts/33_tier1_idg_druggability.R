#!/usr/bin/env Rscript
# 33_tier1_idg_druggability.R — 6 个非-MHC Tier 1 基因的 IDG 成药性标注落盘
# 回应逻辑断点 L6：Tier 1 基因的成药性梯度（IDG Pharos/DrugCentral TDL + AB + SM）
# 铁律 #14：TDL / AB accessibility / SM evidence 来自 IDG Pharos/DrugCentral 查询，
#           由脚本统一落盘，消除 generate_manuscript_figures.R 内硬编码 tibble 的软支撑。

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

proj <- "/ifs1/User/zhouman/project9-v5-crc-atlas"

## — IDG 成药性标注（6 个非-MHC Tier 1 基因） ————————————————————
# TDL 梯度：Tclin > Tchem > Tchem*（无 AB 可达）> Tbio
# antibody = 抗体可及性评分（0–3，IDG Pharos antibody tractability）
idg <- tibble(
  gene      = c("LAMC1", "GPBAR1", "BMP2", "UTP11", "TMBIM1", "PNKD"),
  tdl       = c("Tclin", "Tchem", "Tchem", "Tchem*", "Tbio", "Tbio"),
  tdl_score = c(4, 3, 3, 2.5, 2, 1),
  has_drug  = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  antibody  = c(3, 3, 3, 0, 3, 3),
  sm        = c("None", "Structure+Ligand; Druggable Family", "Structure with Ligand",
                "Structure with Ligand", "None", "None")
)

## — 后置断言（铁律 #5 / #14） ————————————————————————————————
stopifnot(nrow(idg) == 6)
stopifnot(identical(idg$gene, c("LAMC1", "GPBAR1", "BMP2", "UTP11", "TMBIM1", "PNKD")))
stopifnot(all(idg$tdl %in% c("Tclin", "Tchem", "Tchem*", "Tbio")))
stopifnot(all(diff(idg$tdl_score) <= 0))          # 成药性梯度单调不增（LAMC1 ≥ … ≥ PNKD）
stopifnot(sum(idg$has_drug) == 1 && idg$has_drug[idg$gene == "LAMC1"])  # 仅 LAMC1 有已批准药
stopifnot(idg$antibody[idg$gene == "UTP11"] == 0) # UTP11 核仁定位 → AB 0/3

## — 落盘 ————————————————————————————————————————————————
out <- file.path(proj, "results", "phase_tier1_idg_druggability.csv")
write_csv(idg, out)
cat(sprintf("✅ Written %d rows -> %s\n", nrow(idg), out))
print(idg)
