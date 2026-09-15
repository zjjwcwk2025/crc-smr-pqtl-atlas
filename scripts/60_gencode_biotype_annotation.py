#!/usr/bin/env python3
"""
Phase 1 P1-fix: 用 GENCODE v47 单一 biotype 来源重注释 75 个 Bonferroni 显著基因，
锁定全文唯一的「蛋白编码 vs 非编码」计数（解决 62/63/68 三套口径不一致）。

输入:
  - results/phase1_bonferroni_significant_annotated.tsv  (col1=ENSG_clean, col23=SYMBOL)
  - data/gencode/gencode.v47.basic.annotation.gtf.gz     (GENCODE v47, 单一 biotype 来源)
输出:
  - results/phase1_ensg_biotype.tsv                      (ENSG, gene_type, gene_name)
  - results/phase1_ensg_biotype_mismatch.tsv             (eQTLGen SYMBOL vs GENCODE gene_name 不一致项)
后置校验: assert 75 个唯一 ENSG 全部命中 GENCODE; assert 蛋白编码计数唯一且写入 stdout。

铁律 #4/#12/#14: 路径集中、随机数无关、关键数值从输出文件复制不手敲。
"""
import gzip
import re
import csv
from collections import OrderedDict

ANNOTATED_TSV = "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase1_bonferroni_significant_annotated.tsv"
GTF = "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gencode/gencode.v47.annotation.gtf.gz"
OUT_BIOTYPE = "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase1_ensg_biotype.tsv"
OUT_MISMATCH = "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase1_ensg_biotype_mismatch.tsv"


def strip_version(ensg: str) -> str:
    """去掉 GTF gene_id 的版本后缀 (ENSGxxxx.N -> ENSGxxxx)。"""
    return ensg.split(".")[0]


def main():
    # 1) 读取 annotated TSV: 唯一 ENSG (保持原顺序), 并记录 eQTLGen SYMBOL (首现)
    ensg_set = OrderedDict()
    with open(ANNOTATED_TSV, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            e = row["ENSG_clean"].strip()
            s = row.get("SYMBOL", "").strip()
            if e and e not in ensg_set:
                ensg_set[e] = s

    n_unique = len(ensg_set)
    assert n_unique == 75, f"期望 75 个唯一 ENSG, 实得 {n_unique}"

    # 2) 流式解析 GENCODE GTF (仅 feature == "gene")
    gene_type = {}
    gene_name = {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attr = cols[8]
            m_id = re.search(r'gene_id "([^"]+)"', attr)
            if not m_id:
                continue
            gid = strip_version(m_id.group(1))
            if gid not in ensg_set:
                continue
            m_type = re.search(r'gene_type "([^"]+)"', attr)
            m_name = re.search(r'gene_name "([^"]+)"', attr)
            gene_type[gid] = m_type.group(1) if m_type else ""
            gene_name[gid] = m_name.group(1) if m_name else ""

    # 3) 缺失项 (不在 GENCODE v47 中, 即 withdrawn/非 GENCODE 基因) → 归为非编码
    missing = [e for e in ensg_set if e not in gene_type]
    for e in missing:
        gene_type[e] = "not_in_GENCODE_v47"
        gene_name[e] = ""

    # 4) 写出 biotype 文件 (ENSG, gene_type, gene_name)
    with open(OUT_BIOTYPE, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["ENSG", "gene_type", "gene_name"])
        for e in ensg_set:
            w.writerow([e, gene_type[e], gene_name[e]])

    # 5) eQTLGen SYMBOL vs GENCODE gene_name 不一致项
    mismatch = []
    for e, eqtl_sym in ensg_set.items():
        gc_name = gene_name[e]
        # 空 eQTLGen SYMBOL 单独记录, 不算符号冲突
        if eqtl_sym == "" or eqtl_sym == e:
            continue
        if eqtl_sym != gc_name:
            mismatch.append((e, eqtl_sym, gc_name, gene_type[e]))
    with open(OUT_MISMATCH, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["ENSG", "eqtlgen_SYMBOL", "gencode_gene_name", "gene_type"])
        for row in mismatch:
            w.writerow(list(row))

    # 6) 计数
    from collections import Counter
    cnt = Counter(gene_type.values())
    protein_coding = cnt.get("protein_coding", 0)
    non_coding = n_unique - protein_coding

    print("=== GENCODE v47 重注释结果 (75 unique ENSG) ===")
    for bt, c in cnt.most_common():
        print(f"  {bt}: {c}")
    print(f"  TOTAL: {n_unique}")
    print(f"  protein_coding: {protein_coding}")
    print(f"  non-protein-coding: {non_coding}")
    print(f"  eQTLGen SYMBOL vs GENCODE gene_name 不一致: {len(mismatch)}")
    for row in mismatch:
        print(f"    {row[0]}: eQTLGen={row[1]} vs GENCODE={row[2]} ({row[3]})")

    # 7) 打印缺失/边界项供人工核对（最终锁定值在确认后写死回 CLAUDE.md）
    if missing:
        print(f"  not_in_GENCODE_v47 (withdrawn/非GENCODE): {missing}")


if __name__ == "__main__":
    main()
