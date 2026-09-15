#!/usr/bin/env python3
"""Phase 5: Drug-target annotation via Open Targets Platform v4 GraphQL API
Uses requests library (NOT subprocess+curl — DNS fails in subprocess).
Output: results/phase5_drug_annotation.csv"""

import requests, json, time, csv, sys

API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# === 1. Load SMR genes ===
smr_file = "results/phase1_bonferroni_significant_annotated.tsv"
genes = []
with open(smr_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        sym = row.get('SYMBOL', '').strip()
        ensg = row.get('ENSG_clean', '').strip()
        ps = row.get('p_SMR', '')
        bs = row.get('b_SMR', '')
        if sym and ensg and ps:
            try:
                genes.append({'gene': sym, 'ensg': ensg,
                              'p_SMR': float(ps), 'b_SMR': float(bs) if bs else 0})
            except: pass

# Load coloc PPH4
coloc_pph4 = {}
try:
    with open("results/phase3_coloc_top10.csv") as f:
        for row in csv.DictReader(f):
            coloc_pph4[row['ensg']] = float(row.get('PPH4', 0))
except: pass

# Deduplicate
seen = set()
unique = []
for g in genes:
    if g['gene'] not in seen:
        seen.add(g['gene'])
        g['PPH4'] = coloc_pph4.get(g['ensg'])
        unique.append(g)

unique.sort(key=lambda x: x['p_SMR'])
genes_q = unique[:30]  # top 30 for speed
print(f"Gene set: {len(unique)} unique, querying top {len(genes_q)}")

# === 2. Query Open Targets v4 ===
def query_one(ensg, gene):
    query = (
        'query {'
        f'  target(ensemblId: "{ensg}") {{'
        '    id'
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

    try:
        r = requests.post(API_URL, json={'query': query}, timeout=30)
        if r.status_code != 200:
            print(f"  {gene}: HTTP {r.status_code}")
            return None
        data = r.json()
    except Exception as e:
        print(f"  {gene}: {e}")
        return None

    target = data.get('data', {}).get('target', {})
    drugs = target.get('drugAndClinicalCandidates', {}).get('rows', [])

    if not drugs:
        print(f"  {gene}: 0 drugs")
        return {
            'gene': gene, 'ensg': ensg, 'n_drugs': 0,
            'top_drug': '', 'top_stage': '', 'top_type': '',
            'top_moa': '', 'top_indication': '',
            'all_drugs': '', 'all_stages': '',
            'n_approved': 0, 'n_phase3': 0, 'n_phase2': 0, 'n_phase1': 0, 'n_phase0': 0
        }

    # Parse drug info
    drug_info = []
    for d in drugs:
        dr = d['drug']
        name = dr.get('name', '?')
        stage = str(dr.get('maximumClinicalStage') or '0')
        dtype = dr.get('drugType', '?')

        moas = dr.get('mechanismsOfAction', {}).get('rows', [])
        moa_list = [m.get('mechanismOfAction','') for m in moas[:2] if m.get('mechanismOfAction')]
        moa_str = '; '.join(moa_list) if moa_list else ''

        inds = dr.get('indications', {}).get('rows', [])
        top_ind = inds[0]['disease']['name'] if inds else ''

        drug_info.append((name, stage, dtype, moa_str, top_ind))

    # Stage mapping
    stage_map = {'4': 'APPROVAL', '3': 'PHASE3', '2': 'PHASE2', '1': 'PHASE1', '0': 'PRECLINICAL'}
    approved_count = sum(1 for _, s, _, _, _ in drug_info if s == 'APPROVAL')
    p3_count = sum(1 for _, s, _, _, _ in drug_info if s == 'PHASE_3')

    top = drug_info[0]
    print(f"  {gene}: {len(drugs)} drugs, top={top[0]} ({top[1]})")

    return {
        'gene': gene, 'ensg': ensg, 'n_drugs': len(drugs),
        'top_drug': top[0], 'top_stage': top[1],
        'top_type': top[2], 'top_moa': top[3][:120],
        'top_indication': top[4],
        'all_drugs': ', '.join([d[0] for d in drug_info]),
        'all_stages': ', '.join([d[1] for d in drug_info]),
        'n_approved': approved_count,
        'n_phase3': p3_count,
        'n_phase2': sum(1 for _, s, _, _, _ in drug_info if s in ('PHASE_2', 'PHASE2')),
        'n_phase1': sum(1 for _, s, _, _, _ in drug_info if s in ('PHASE_1', 'PHASE1')),
        'n_phase0': sum(1 for _, s, _, _, _ in drug_info if s in ('0', 'PRECLINICAL', 'PHASE_0', 'PHASE0'))
    }

# === 3. Batch query ===
results = []
for i, g in enumerate(genes_q):
    print(f"[{i+1}/{len(genes_q)}]", end=" ")
    r = query_one(g['ensg'], g['gene'])
    if r:
        r['p_SMR'] = g['p_SMR']
        r['b_SMR'] = g['b_SMR']
        r['PPH4'] = g.get('PPH4')
        results.append(r)
    time.sleep(0.3)

# === 4. Save ===
outfile = "results/phase5_drug_annotation.csv"
fields = ['gene','ensg','p_SMR','b_SMR','PPH4','n_drugs','top_drug',
          'top_stage','top_type','top_moa','top_indication',
          'all_drugs','all_stages',
          'n_approved','n_phase3','n_phase2','n_phase1','n_phase0']
with open(outfile, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(results)

# === 5. Summary ===
print(f"\n{'='*80}")
print(f"SUMMARY: {len(results)} genes queried")
druggable = [r for r in results if r['n_drugs'] > 0]
approved = [r for r in results if r['n_approved'] > 0]
print(f"With known drugs: {len(druggable)}/{len(results)}")
print(f"With approved drugs: {len(approved)}")
print(f"With Phase 3+ pipeline: {sum(1 for r in results if r['n_phase3'] > 0)}")

print(f"\n{'='*80}")
print("TOP DRUGGABLE TARGETS:")
print(f"{'Gene':<10} {'p_SMR':<12} {'PPH4':<8} {'Drugs':<6} {'Top Drug':<22} {'Stage':<14} {'Type':<12} {'Indication'}")
print("-"*110)
for r in sorted(druggable, key=lambda x: x['p_SMR']):
    pph4_str = f"{r['PPH4']:.2f}" if r['PPH4'] is not None else "N/A"
    print(f"{r['gene']:<10} {r['p_SMR']:<12.2e} {pph4_str:<8} {r['n_drugs']:<6} "
          f"{r['top_drug'][:20]:<22} {r['top_stage'][:12]:<14} {r['top_type'][:10]:<12} "
          f"{r['top_indication'][:50]}")

# Drug classification
print(f"\n{'='*80}")
print("DRUG SCORE CLASSIFICATION:")
for r in sorted(results, key=lambda x: x['p_SMR']):
    if r['n_approved'] > 0:
        score = "APPROVED (repurposing)"
    elif r['n_phase3'] > 0:
        score = "Phase 3 pipeline"
    elif r['n_drugs'] > 0:
        score = "Early pipeline"
    else:
        score = "Novel biology"

    is_crc = any(kw in (r.get('top_indication') or '').lower()
                 for kw in ['colorect','colon','rectal','bowel','crc'])
    crc_flag = " [CRC]" if is_crc else ""
    if r['n_drugs'] > 0:
        print(f"  {r['gene']:<10} {score:<25} {r['n_drugs']} drugs{crc_flag}")

print(f"\nSaved: {outfile}")
