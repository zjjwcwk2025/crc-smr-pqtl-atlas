#!/usr/bin/env python3
"""
Phase 5b: Deep Drug Annotation + Drug Repurposing Score
Data source: Open Targets Platform v4 GraphQL API
DGIdb blocked, OT Genetics DNS fails — skipped gracefully
Integrates: SMR, HEIDI, coloc, SuSiE, TCGA, scRNA for repurposing score
Output: results/phase5b_deep_drug_annotation.csv
"""

import requests, json, time, csv, math, os
from collections import defaultdict

PROJ = "/ifs1/User/zhouman/project9-v5-crc-atlas"
os.chdir(PROJ)

# ── 1. Load all existing results ─────────────────────────────────────────────

print("[1/4] Loading existing results...")

def load_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))

def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))

smr_rows = load_tsv("results/phase1_bonferroni_significant_annotated.tsv")
coloc_rows = load_csv("results/phase3_coloc_top10.csv")
susie_rows = load_csv("results/phase3_susie_finemap.csv")
tcga_rows = load_csv("results/phase6_scrna/TableS_tcga_tumor_vs_normal.csv")
scrna_rows = load_csv("results/phase6_scrna/phase6_celltype_specificity.csv")

coloc_by_ensg = {r['ensg']: r for r in coloc_rows}
susie_by_ensg = {r['ensg']: r for r in susie_rows}
tcga_by_gene = {r['gene']: r for r in tcga_rows}
scrna_by_gene = {r['gene_symbol']: r for r in scrna_rows}

genes = {}
for r in smr_rows:
    sym = r.get('SYMBOL', '').strip()
    ensg = r.get('ENSG_clean', '').strip()
    if not sym or not ensg:
        continue
    p = float(r.get('p_SMR', '1'))
    b = float(r.get('b_SMR', '0'))
    chr_val = r.get('ProbeChr', '')
    pos_val = r.get('Probe_bp', '')
    p_heidi = r.get('p_HEIDI', 'NA')
    try:
        p_heidi = float(p_heidi)
    except:
        p_heidi = None

    if sym not in genes or p < genes[sym]['p_SMR']:
        genes[sym] = {
            'gene': sym, 'ensg': ensg, 'chr': chr_val, 'pos': pos_val,
            'p_SMR': p, 'b_SMR': b, 'p_HEIDI': p_heidi,
        }

for g in genes.values():
    # Defaults
    g.update({
        'PPH4': None, 'PPH3': None, 'n_cs': None, 'top_pip_snp': None, 'tier': None,
        'log2FC': None, 'tcga_pval': None, 'tcga_sig': None, 'tcga_dir': None,
        'top_celltype': None, 'specificity': None, 'expr_fibro': None, 'expr_epithelial': None,
        'heidi_pass': (g['p_HEIDI'] is not None and g['p_HEIDI'] >= 0.01),
    })
    # Coloc
    c = coloc_by_ensg.get(g['ensg'])
    if c:
        g['PPH4'] = float(c.get('PPH4', 0))
        g['PPH3'] = float(c.get('PPH3', 0))
    # SuSiE
    s = susie_by_ensg.get(g['ensg'])
    if s:
        g['n_cs'] = int(s.get('n_cs', 0))
        g['top_pip_snp'] = s.get('top_pip_snp', '')
        g['tier'] = int(s.get('tier', 0))
    # TCGA
    t = tcga_by_gene.get(g['gene'])
    if t:
        g['log2FC'] = float(t.get('log2FC', 0))
        g['tcga_pval'] = float(t.get('p_value', 1))
        g['tcga_sig'] = t.get('significant', '')
        g['tcga_dir'] = t.get('direction', '')
    # scRNA
    sc = scrna_by_gene.get(g['gene'])
    if sc:
        g['top_celltype'] = sc.get('top_celltype', '')
        try:
            g['specificity'] = float(sc.get('specificity', 0) or 0)
        except (ValueError, TypeError):
            g['specificity'] = 0
        try:
            g['expr_fibro'] = float(sc.get('expr_fibro', 0) or 0)
        except (ValueError, TypeError):
            g['expr_fibro'] = 0
        try:
            g['expr_epithelial'] = float(sc.get('expr_epithelial', 0) or 0)
        except (ValueError, TypeError):
            g['expr_epithelial'] = 0

gene_list = sorted(genes.values(), key=lambda x: x['p_SMR'])
n_coloc = sum(1 for g in gene_list if g['PPH4'] is not None and g['PPH4'] > 0.8)
n_tier1 = sum(1 for g in gene_list if g['tier'] and g['tier'] == 1)
print(f"  {len(gene_list)} genes loaded (coloc-passing: {n_coloc}, SuSiE Tier 1: {n_tier1})")

# ── 2. Query Open Targets Platform for ALL genes ─────────────────────────────

print("\n[2/4] Querying Open Targets Platform...")

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

def query_ot(ensg, gene_sym):
    query = (
        'query {'
        f'  target(ensemblId: "{ensg}") {{'
        '    approvedSymbol'
        '    drugAndClinicalCandidates {'
        '      rows {'
        '        drug {'
        '          name'
        '          maximumClinicalStage'
        '          drugType'
        '          mechanismsOfAction { rows { mechanismOfAction } }'
        '          indications { rows { disease { name } } }'
        '        }'
        '      }'
        '    }'
        '  }'
        '}'
    )
    time.sleep(0.15)
    try:
        r = requests.post(OT_URL, json={'query': query}, timeout=30)
        if r.status_code != 200:
            return None
        data = r.json()
    except Exception as e:
        print(f"  {gene_sym}: error - {e}")
        return None

    target = data.get('data', {}).get('target')
    if not target:
        return None

    drugs = target.get('drugAndClinicalCandidates', {}).get('rows') or []
    if not drugs:
        return {'n_drugs': 0, 'n_drugs_unique': 0, 'n_approved': 0, 'n_phase3': 0,
                'n_phase2': 0, 'n_phase1': 0, 'top_drug': '', 'top_stage': '', 'top_type': '',
                'top_moa': '', 'top_indication': '', 'all_drug_names': '',
                'has_crc_drug': False, 'crc_drug_text': ''}

    # Parse drug info
    drug_info = []
    for d in drugs:
        dr = d.get('drug') or {}
        name = dr.get('name') or dr.get('id', '?')
        stage = str(dr.get('maximumClinicalStage') or '0')

        # Stage mapping for ranking
        stage_map = {'APPROVAL': 4, 'PHASE_3': 3, 'PHASE_2': 2, 'PHASE_1': 1}
        stage_num = stage_map.get(stage, 0)

        dtype = dr.get('drugType', '?')

        moas = dr.get('mechanismsOfAction', {}).get('rows', [])
        moa_str = '; '.join([m.get('mechanismOfAction', '') for m in moas[:2] if m.get('mechanismOfAction')])

        inds = dr.get('indications', {}).get('rows', [])
        top_ind = inds[0]['disease']['name'] if inds else ''

        drug_info.append({
            'name': name, 'stage': stage, 'stage_num': stage_num,
            'type': dtype, 'moa': moa_str, 'indication': top_ind,
        })

    # Sort by stage desc
    drug_info.sort(key=lambda x: x['stage_num'], reverse=True)

    n_approved = sum(1 for d in drug_info if d['stage'] == 'APPROVAL')
    n_phase3 = sum(1 for d in drug_info if d['stage'] == 'PHASE_3')
    n_phase2 = sum(1 for d in drug_info if d['stage'] == 'PHASE_2')
    n_phase1 = sum(1 for d in drug_info if d['stage'] == 'PHASE_1')

    top = drug_info[0]
    all_names = ' | '.join(d['name'] for d in drug_info)

    # CRC drugs
    crc_kw = ['colorectal', 'colon', 'rectal', 'bowel', 'crc']
    crc_drugs = []
    for d in drug_info:
        ind_lower = d['indication'].lower()
        if any(kw in ind_lower for kw in crc_kw):
            crc_drugs.append(d['name'])

    return {
        'n_drugs': len(drug_info), 'n_drugs_unique': len(set(d['name'] for d in drug_info)),
        'n_approved': n_approved, 'n_phase3': n_phase3, 'n_phase2': n_phase2, 'n_phase1': n_phase1,
        'top_drug': top['name'], 'top_stage': top['stage'],
        'top_type': top['type'], 'top_moa': top['moa'], 'top_indication': top['indication'],
        'all_drug_names': all_names,
        'has_crc_drug': len(crc_drugs) > 0,
        'crc_drug_text': '; '.join(crc_drugs) if crc_drugs else '',
    }

n_total = len(gene_list)
for i, g in enumerate(gene_list):
    r = query_ot(g['ensg'], g['gene'])
    if r and r.get('n_drugs', 0) > 0:
        g.update(r)
        print(f"  [{i+1}/{n_total}] {g['gene']}: {g['n_drugs']} drugs ({g['n_approved']} approved)")
    else:
        g.update({
            'n_drugs': 0, 'n_drugs_unique': 0,
            'n_approved': 0, 'n_phase3': 0, 'n_phase2': 0, 'n_phase1': 0,
            'top_drug': '', 'top_stage': '', 'top_type': '',
            'top_moa': '', 'top_indication': '', 'all_drug_names': '',
            'has_crc_drug': False, 'crc_drug_text': '',
        })

n_druggable = sum(1 for g in gene_list if g['n_drugs'] > 0)
n_approved  = sum(1 for g in gene_list if g['n_approved'] > 0)
n_crc       = sum(1 for g in gene_list if g['has_crc_drug'])
print(f"\n  Druggable: {n_druggable}/{n_total} ({100*n_druggable/n_total:.0f}%)")
print(f"  Approved drugs: {n_approved} genes")
print(f"  CRC-specific: {n_crc} genes")

# ── 3. Drug Repurposing Potential Score ──────────────────────────────────────

print("\n[3/4] Computing drug repurposing scores...")

max_neg_log_p = -math.log10(min(g['p_SMR'] for g in gene_list))

# Data-driven normalization: max values from actual data, not hardcoded constants
max_log2FC = max(abs(g.get('log2FC', 0) or 0) for g in gene_list)
if max_log2FC == 0:
    max_log2FC = 1.0  # prevent division by zero
max_spec = max(g.get('specificity', 0) or 0 for g in gene_list)
if max_spec == 0:
    max_spec = 1.0

for g in gene_list:
    score_smr       = min(-math.log10(g['p_SMR']) / max_neg_log_p, 1)
    score_coloc     = g['PPH4'] if g['PPH4'] else 0
    score_susie     = 1.0 if (g.get('tier') and g['tier'] == 1) else 0
    score_heidi     = 1.0 if g['heidi_pass'] else 0.3

    lfc = g.get('log2FC')
    score_tcga      = min(abs(lfc) / max_log2FC, 1) if lfc else 0

    spec = g.get('specificity')
    score_celltype  = min(spec / max_spec, 1) if spec else 0

    n_app = g.get('n_approved', 0) or 0
    n_p3  = g.get('n_phase3', 0) or 0
    n_d   = g.get('n_drugs', 0) or 0
    # Drug evidence scoring rubric (definitional, not data-derived):
    # Tiered by clinical trial phase — approved > Phase 3 > other > none
    if n_app > 0:
        score_drugs = 1.0
    elif n_p3 > 0:
        score_drugs = 0.7
    elif n_d > 0:
        score_drugs = 0.4
    else:
        score_drugs = 0

    # Novelty bonus: strong genetics with no existing drugs
    if score_drugs == 0 and (score_coloc > 0.8 or score_susie == 1):
        score_novelty = 0.5
    else:
        score_novelty = 0

    # Composite weights: defined rubric balancing evidence strength
    # Genetic evidence (SMR + coloc + SuSiE + HEIDI) weighted most heavily
    score_genetic  = score_smr * 0.20 + score_coloc * 0.35 + score_susie * 0.30 + score_heidi * 0.15
    # Biological validation (TCGA + scRNA + novelty)
    score_biology  = score_tcga * 0.50 + score_celltype * 0.30 + score_novelty * 0.20
    # Final: genetic (40%), biology (25%), drug pipeline (35%)
    repurposing    = score_genetic * 0.40 + score_biology * 0.25 + score_drugs * 0.35

    if repurposing >= 0.70:
        tier = "A (High — multi-evidence + approved drugs)"
    elif repurposing >= 0.50:
        tier = "B (Strong — coloc/SuSiE + pipeline drugs)"
    elif repurposing >= 0.30:
        tier = "C (Moderate — genetic evidence)"
    else:
        tier = "D (Exploratory — limited evidence)"

    g.update({
        'score_smr': score_smr, 'score_coloc': score_coloc,
        'score_susie': score_susie, 'score_heidi': score_heidi,
        'score_tcga': score_tcga, 'score_celltype': score_celltype,
        'score_drugs': score_drugs,
        'score_genetic': score_genetic, 'score_biology': score_biology,
        'repurposing_score': repurposing, 'drug_tier': tier,
    })

tiers = defaultdict(int)
for g in gene_list:
    t = g['drug_tier']
    if t.startswith('A'): tiers['A'] += 1
    elif t.startswith('B'): tiers['B'] += 1
    elif t.startswith('C'): tiers['C'] += 1
    else: tiers['D'] += 1
print(f"  Tier A: {tiers['A']}, B: {tiers['B']}, C: {tiers['C']}, D: {tiers['D']}")

gene_list.sort(key=lambda x: x['repurposing_score'], reverse=True)

# ── 4. Save output ───────────────────────────────────────────────────────────

print("\n[4/4] Saving results...")

fields = [
    'gene', 'ensg', 'chr', 'pos',
    'p_SMR', 'b_SMR', 'p_HEIDI', 'heidi_pass',
    'PPH4', 'PPH3', 'n_cs', 'top_pip_snp', 'tier',
    'log2FC', 'tcga_pval', 'tcga_sig', 'tcga_dir',
    'top_celltype', 'specificity', 'expr_fibro', 'expr_epithelial',
    'n_drugs', 'n_drugs_unique', 'n_approved', 'n_phase3', 'n_phase2', 'n_phase1',
    'top_drug', 'top_stage', 'top_type', 'top_moa', 'top_indication',
    'all_drug_names', 'has_crc_drug', 'crc_drug_text',
    'score_smr', 'score_coloc', 'score_susie', 'score_heidi',
    'score_tcga', 'score_celltype', 'score_drugs',
    'score_genetic', 'score_biology',
    'repurposing_score', 'drug_tier',
]

with open("results/phase5b_deep_drug_annotation.csv", 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(gene_list)

# ── 5. Summary report ────────────────────────────────────────────────────────

print(f"\n{'='*90}")
print("  DEEP DRUG ANNOTATION + REPURPOSING SCORE")
print(f"{'='*90}\n")

print("── Highest Repurposing Potential (Tier A+B) ──\n")
top = [g for g in gene_list if g['repurposing_score'] >= 0.50]
for g in top[:20]:
    col_str = "coloc✓" if (g['PPH4'] and g['PPH4'] > 0.8) else "coloc✗"
    hei_str = "HEIDI✓" if g['heidi_pass'] else "HEIDI✗"
    if g['n_approved'] > 0:
        drug_str = f"APPROVED({g['n_approved']})"
    elif g['n_drugs'] > 0:
        drug_str = f"{g['n_drugs']} drugs"
    else:
        drug_str = "none"
    crc_flag = " [CRC]" if g['has_crc_drug'] else ""

    print(f"  {g['gene']:<12} score={g['repurposing_score']:.3f}  "
          f"SMR={g['p_SMR']:.1e}  |  {col_str} {hei_str}  |  {drug_str}{crc_flag}")
    if g.get('top_drug'):
        print(f"    → {g['top_drug']} ({g['top_stage']}) — {g['top_indication']}")
    if g.get('top_moa'):
        print(f"    MOA: {g['top_moa']}")

print("\n── Drug Pipeline Summary ──\n")
print(f"  Approved drugs:         {sum(1 for g in gene_list if g['n_approved'] > 0)} genes")
print(f"  Phase 3:                {sum(1 for g in gene_list if g['n_phase3'] > 0)} genes")
print(f"  Phase 2:                {sum(1 for g in gene_list if g['n_phase2'] > 0)} genes")
print(f"  Any pipeline drug:      {sum(1 for g in gene_list if g['n_drugs'] > 0)} genes")
print(f"  CRC-specific drugs:     {sum(1 for g in gene_list if g['has_crc_drug'])} genes")

print("\n── Genes with Approved Drugs ──\n")
approved = sorted([g for g in gene_list if g['n_approved'] > 0],
                  key=lambda x: x['n_approved'], reverse=True)
for g in approved:
    crc_str = f" [CRC drug: {g['crc_drug_text']}]" if g['has_crc_drug'] else ""
    print(f"  {g['gene']}: {g['n_approved']} approved drugs, top={g['top_drug']} → {g['top_indication']}{crc_str}")

print("\n── Novel Targets (coloc/SuSiE + no existing drugs) ──\n")
novel = [g for g in gene_list if g['n_drugs'] == 0 and (g.get('PPH4') is not None and g['PPH4'] > 0.8 or (g.get('tier') and g['tier'] == 1))]
for g in novel:
    print(f"  {g['gene']:<12} coloc={g['PPH4']:.3f}  SuSiE={g.get('tier','-')}  SMR={g['p_SMR']:.1e}  celltype={g.get('top_celltype','-')}")

print(f"\nFull results: results/phase5b_deep_drug_annotation.csv ({len(gene_list)} genes)")
print("Done.")
