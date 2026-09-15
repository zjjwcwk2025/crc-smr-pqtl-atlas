#!/usr/bin/env python3
"""
Compute HEIDI test statistics directly from 1000G VCF for the 75 significant SMR genes.
Bypasses the need for PLINK format by reading VCF directly with cyvcf2.

HEIDI test (Heterogeneity In Dependent Instruments):
  Tests whether the SMR effect at the top eQTL SNP is consistent with
  effects at other SNPs in the cis-region. Rejects if multiple causal
  variants exist (pleiotropy).

  For each SNP i in cis-region (±1Mb of probe):
    d_i = b_GWAS_i - (b_GWAS_top / b_eQTL_top) * b_eQTL_i
    var(d_i) = var(b_GWAS_i) + (b_GWAS_top/b_eQTL_top)² * var(b_eQTL_i)
               - 2*(b_GWAS_top/b_eQTL_top)*cov(b_GWAS_i, b_eQTL_i)
               + terms involving LD

  Simplified: d_i = b_eQTL_i * b_GWAS_top - b_eQTL_top * b_GWAS_i
  Under null of single causal variant, d_i ~ N(0, V)
  V_ij depends on LD between SNPs i and j.

  T_HEIDI = Σ_i Σ_j d_i * V^{-1}_ij * d_j ~ χ² with (m-1) df
  where m = number of SNPs in cis-region

Reference: Zhu et al. 2016, Nature Genetics
"""
import cyvcf2
import numpy as np
import os
import sys
import gzip

# Paths
VCF_DIR = "/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/vcf"
PANEL_FILE = "/ifs1/User/zhouman/project9-v5-crc-atlas/data/1kg_phase3/integrated_call_samples_v3.20130502.ALL.panel"
SMR_RESULTS = "/ifs1/User/zhouman/project9-v5-crc-atlas/results/phase1_bonferroni_significant_annotated.tsv"
GWAS_MA = "/ifs1/User/zhouman/project9-v5-crc-atlas/data/gwas/GCST90255675.ma"
OUTDIR = "/ifs1/User/zhouman/project9-v5-crc-atlas/results"

CIS_WINDOW = 1_000_000  # ±1Mb around probe


def load_eur_samples():
    """Load EUR sample IDs from 1000G panel."""
    eur_ids = set()
    with open(PANEL_FILE) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[2] == 'EUR':
                eur_ids.add(parts[0])
    print(f"Loaded {len(eur_ids)} EUR samples")
    return eur_ids


def load_smr_genes():
    """Load the 75 significant SMR genes with their info."""
    genes = []
    with open(SMR_RESULTS) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 20:
                continue
            genes.append({
                'ensg': parts[0],
                'probe_chr': int(parts[2]),
                'probe_bp': int(parts[4]),
                'top_snp': parts[5],
                'top_snp_chr': int(parts[6]),
                'top_snp_bp': int(parts[7]),
                'a1': parts[8], 'a2': parts[9],
                'freq': float(parts[10]) if parts[10] else 0.5,
                'b_gwas': float(parts[11]),
                'se_gwas': float(parts[12]),
                'p_gwas': float(parts[13]),
                'b_eqtl': float(parts[14]),
                'se_eqtl': float(parts[15]),
                'p_eqtl': float(parts[16]),
                'b_smr': float(parts[17]),
                'se_smr': float(parts[18]),
                'p_smr': float(parts[19]),
                'symbol': parts[22] if len(parts) > 22 else parts[0]
            })
    print(f"Loaded {len(genes)} SMR-significant genes")
    return genes


def compute_ld(vcf_reader, top_snp_id, cis_snps, eur_indices):
    """Compute LD (r²) between top SNP and each cis SNP."""
    n = len(eur_indices)
    r2 = {}

    # Get top SNP genotype
    top_gt = None
    for variant in vcf_reader:
        if variant.ID == top_snp_id:
            gts = np.zeros(n)
            for i, idx in enumerate(eur_indices):
                gt = variant.genotypes[idx]
                if gt[0] == -1 or gt[1] == -1:  # missing
                    gts[i] = np.nan
                else:
                    gts[i] = gt[0] + gt[1]  # 0, 1, 2
            top_gt = gts
            break

    if top_gt is None:
        return r2

    top_valid = ~np.isnan(top_gt)
    top_std = np.nanstd(top_gt)
    if top_std == 0:
        return r2

    # Reset reader to beginning of cis region
    # This is tricky - cyvcf2 doesn't easily reset. We'll use a region query.
    return r2


def compute_heidi_simple(gene, ld_r2_dict, cis_snps_info):
    """
    Simplified HEIDI test using the SMR approach.

    Uses the multi-SNP HEIDI test:
    For each SNP j, compute the deviation:
      d_j = b_eQTL_j * b_SMR - b_GWAS_j  (in standardized units)

    Under null (no pleiotropy), d_j ~ N(0, var(d_j))
    where var(d_j) depends on LD with top SNP.

    Returns: p_HEIDI, n_snps_used
    """
    b_smr = gene['b_smr']
    se_smr = gene['se_smr']
    b_eqtl_top = gene['b_eqtl']
    b_gwas_top = gene['b_gwas']

    # For now, return placeholder — full implementation needs LD computation
    return 1.0, 0


def main():
    eur_ids = load_eur_samples()
    genes = load_smr_genes()

    # Group genes by chromosome for efficient VCF reading
    genes_by_chr = {}
    for g in genes:
        ch = g['probe_chr']
        if ch not in genes_by_chr:
            genes_by_chr[ch] = []
        genes_by_chr[ch].append(g)

    print(f"Genes distributed across chromosomes: {dict((k, len(v)) for k, v in sorted(genes_by_chr.items()))}")

    results = []

    for chr_num, chr_genes in sorted(genes_by_chr.items()):
        vcf_path = os.path.join(VCF_DIR, f"chr{chr_num}.vcf.gz")
        if not os.path.exists(vcf_path):
            print(f"WARNING: VCF not found for chr{chr_num}: {vcf_path}")
            continue

        print(f"\nProcessing chr{chr_num}: {len(chr_genes)} genes...")
        vcf = cyvcf2.VCF(vcf_path)
        eur_indices = []
        for i, sample in enumerate(vcf.samples):
            if sample in eur_ids:
                eur_indices.append(i)
        print(f"  {len(eur_indices)} EUR samples in VCF")

        for gene in chr_genes:
            symbol = gene['symbol']
            probe_chr = gene['probe_chr']
            probe_bp = gene['probe_bp']
            top_snp = gene['top_snp']
            region_start = max(0, probe_bp - CIS_WINDOW)
            region_end = probe_bp + CIS_WINDOW
            region_str = f"{probe_chr}:{region_start}-{region_end}"

            # Query cis-region
            try:
                cis_variants = list(vcf(region_str))
            except Exception as e:
                print(f"  {symbol}: ERROR querying region {region_str}: {e}")
                results.append({
                    'gene': symbol, 'ensg': gene['ensg'], 'chr': probe_chr,
                    'n_cis_snps': 0, 'n_heidi_snps': 0,
                    'p_HEIDI': np.nan, 'status': 'ERROR'
                })
                continue

            n_cis = len(cis_variants)

            # Find top SNP
            top_found = False
            top_b_eqtl = gene['b_eqtl']
            top_b_gwas = gene['b_gwas']

            # Extract genotypes for EUR samples
            n_eur = len(eur_indices)
            genotypes = np.zeros((n_cis, n_eur), dtype=np.int8)
            snp_ids = []
            snp_positions = []

            for vi, variant in enumerate(cis_variants):
                snp_ids.append(variant.ID)
                snp_positions.append(variant.POS)
                for si, ei in enumerate(eur_indices):
                    gt = variant.genotypes[ei]
                    if gt[0] == -1 or gt[1] == -1:
                        genotypes[vi, si] = -1  # missing
                    else:
                        genotypes[vi, si] = gt[0] + gt[1]

                if variant.ID == top_snp:
                    top_idx = vi
                    top_found = True

            if not top_found:
                print(f"  {symbol}: top SNP {top_snp} not found in cis region (n_cis={n_cis})")
                results.append({
                    'gene': symbol, 'ensg': gene['ensg'], 'chr': probe_chr,
                    'n_cis_snps': n_cis, 'n_heidi_snps': 0,
                    'p_HEIDI': np.nan, 'status': 'TOP_SNP_NOT_FOUND'
                })
                continue

            # Filter to SNPs with MAF > 0.01 and low missingness
            top_gt = genotypes[top_idx]
            valid = (top_gt >= 0)
            top_maf = np.mean(top_gt[valid]) / 2
            if top_maf > 0.5:
                top_maf = 1 - top_maf

            # Compute LD (r) between top SNP and all others
            n_used = 0
            for vi in range(n_cis):
                if vi == top_idx:
                    continue
                gt_i = genotypes[vi]
                valid_pair = (top_gt >= 0) & (gt_i >= 0)
                n_valid = np.sum(valid_pair)
                if n_valid < 100:  # too few samples
                    continue

                # Compute r²
                top_std = np.std(top_gt[valid_pair])
                i_std = np.std(gt_i[valid_pair])
                if top_std == 0 or i_std == 0:
                    continue

                corr = np.corrcoef(top_gt[valid_pair], gt_i[valid_pair])[0, 1]

                # Filter by LD: only keep SNPs in moderate LD with top SNP
                # (HEIDI uses all SNPs, but high-LD SNPs dominate the test)
                # For now, count and skip (full HEIDI computation needs GWAS effects)
                n_used += 1

            # HEIDI requires GWAS summary stats for each SNP in the region.
            # These are not available without running GWAS on each SNP.
            # The SMR HEIDI test uses the GWAS summary stats (b, se) from the
            # original GWAS file for each SNP, combined with eQTL effects.
            #
            # Since we have the full GWAS summary stats in GCST90255675.ma,
            # we can match SNPs and compute HEIDI.
            #
            # For now, note the cis SNP counts and proceed.

            print(f"  {symbol}: n_cis={n_cis}, top_found={top_found}, "
                  f"top_maf={top_maf:.3f}, n_pairwise={n_used}")

            results.append({
                'gene': symbol, 'ensg': gene['ensg'], 'chr': probe_chr,
                'top_snp': top_snp, 'n_cis_snps': n_cis,
                'n_ld_computed': n_used, 'top_maf': top_maf,
                'p_HEIDI': np.nan, 'status': 'LD_COMPUTED_NEED_GWAS_MATCH'
            })

        vcf.close()

    # Write results
    outpath = os.path.join(OUTDIR, "phase1_heidi_prelim.tsv")
    with open(outpath, 'w') as f:
        f.write("gene\tensg\tchr\ttop_snp\tn_cis_snps\tn_ld_computed\ttop_maf\tp_HEIDI\tstatus\n")
        for r in results:
            f.write(f"{r['gene']}\t{r['ensg']}\t{r['chr']}\t{r.get('top_snp','')}\t"
                    f"{r['n_cis_snps']}\t{r.get('n_ld_computed','')}\t{r.get('top_maf','')}\t"
                    f"{r['p_HEIDI']}\t{r['status']}\n")

    print(f"\nResults saved to {outpath}")
    print(f"Summary: {len(results)} genes processed")
    status_counts = {}
    for r in results:
        s = r['status']
        status_counts[s] = status_counts.get(s, 0) + 1
    for s, c in status_counts.items():
        print(f"  {s}: {c}")


if __name__ == "__main__":
    main()
