#!/usr/bin/env python3
"""
Phase 5c: PheWAS Safety Screening for Tier 1 genes
Data source: Open Targets Platform v4 GraphQL API
Queries gene->disease associations + variant->evidence for SuSiE credible set SNPs
Flags potential safety concerns (pleiotropy, off-target disease associations)

Strategy:
1. For each Tier 1 gene, query OT Platform for all associated diseases
2. Classify diseases into safety categories
3. Compute safety concern score
4. Annotate with variant-level evidence where available
5. Output summary table + detailed risk report
"""

import requests
import json
import time
import csv
import os
import sys
import re
from collections import defaultdict

PROJ = "/ifs1/User/zhouman/project9-v5-crc-atlas"
os.chdir(PROJ)
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# ── 1. Load Tier 1 gene data ─────────────────────────────────────────────────
print("[1/5] Loading Tier 1 gene + SuSiE data...")

susie_rows = []
with open("results/phase3_susie_finemap.csv") as f:
    susie_rows = list(csv.DictReader(f))

tier1_genes = []
all_snps_by_gene = {}  # gene -> list of credible set SNPs

for r in susie_rows:
    if int(r.get('tier', '0')) != 1:
        continue
    gene = r['gene']
    ensg = r['ensg']
    pph4 = float(r['PPH4'])
    top_snp = r['top_pip_snp']

    # Parse credible set SNPs from cs_summary
    cs_summary = r.get('cs_summary', '')
    cs_snps = []
    for part in cs_summary.split('|'):
        part = part.strip()
        # Format: "CS1(1_SNPs,rsXXXXX,PIP=1.00)"
        m = re.search(r'rs\d+', part)
        if m:
            cs_snps.append(m.group())

    all_snps_by_gene[gene] = list(set(cs_snps))  # dedup

    tier1_genes.append({
        'gene': gene, 'ensg': ensg, 'PPH4': pph4,
        'top_snp': top_snp,
        'n_cs': int(r.get('n_cs', 0)),
        'all_cs_snps': all_snps_by_gene[gene],
    })

print(f"  Tier 1 genes: {len(tier1_genes)}")
for g in tier1_genes:
    print(f"    {g['gene']}: {len(g['all_cs_snps'])} credible set SNPs ({g['n_cs']} CSs)")

# ── 2. Disease classification system ─────────────────────────────────────────

# Categories for safety concern scoring
SAFETY_CATEGORIES = {
    'autoimmune_inflammatory': {
        'keywords': ['rheumatoid arthritis', 'lupus', 'psoriasis', 'inflammatory bowel',
                     "crohn's", 'ulcerative colitis', 'ankylosing spondylitis',
                     'multiple sclerosis', 'type 1 diabetes', 'autoimmune',
                     'systemic lupus', 'sjogren', 'scleroderma', 'vasculitis',
                     'celiac', 'myasthenia', 'pemphigus', 'sarcoidosis',
                     'eczema', 'atopic dermatitis', 'asthma', 'allergy',
                     'hypersensitivity', 'anaphylaxis'],
        'label': 'Autoimmune/Inflammatory',
        'risk_score': 3,
    },
    'cardiovascular': {
        'keywords': ['coronary', 'myocardial', 'heart failure', 'cardiomyopathy',
                     'atrial fibrillation', 'arrhythmia', 'hypertension',
                     'stroke', 'thrombosis', 'embolism', 'pulmonary embolism',
                     'deep vein thrombosis', 'venous thromboembolism',
                     'aortic aneurysm', 'peripheral arterial', 'sudden cardiac',
                     'ventricular', 'long qt', 'brugada', 'qt interval'],
        'label': 'Cardiovascular',
        'risk_score': 3,
    },
    'neurological_psychiatric': {
        'keywords': ['alzheimer', 'parkinson', 'dementia', 'cognitive',
                     'schizophrenia', 'bipolar', 'depression', 'anxiety',
                     'autism', 'adhd', 'seizure', 'epilepsy', 'migraine',
                     'neuropathy', 'amyotrophic lateral sclerosis', 'huntington',
                     'neurologic', 'neurodegenerative', 'psychiatric',
                     'suicid', 'substance abuse', 'addiction'],
        'label': 'Neurological/Psychiatric',
        'risk_score': 2,
    },
    'metabolic_endocrine': {
        'keywords': ['diabetes', 'obesity', 'metabolic syndrome', 'hyperlipidemia',
                     'hypercholesterolemia', 'hypothyroidism', 'hyperthyroidism',
                     'osteoporosis', 'gout', 'hyperuricemia', 'vitamin',
                     'hypoglycemia', 'hyperglycemia', 'insulin',
                     'lipid metabolism', 'adiposity', 'bmi', 'body mass'],
        'label': 'Metabolic/Endocrine',
        'risk_score': 2,
    },
    'renal_hepatic': {
        'keywords': ['chronic kidney', 'renal failure', 'nephropathy',
                     'glomerulonephritis', 'nephritis', 'kidney disease',
                     'liver failure', 'cirrhosis', 'hepatitis',
                     'cholestasis', 'jaundice', 'hepatic',
                     'non-alcoholic fatty liver', 'nash', 'steatohepatitis',
                     'gallstone', 'cholelithiasis', 'pancreatitis'],
        'label': 'Renal/Hepatic',
        'risk_score': 2,
    },
    'hematological': {
        'keywords': ['anemia', 'thrombocytopenia', 'neutropenia', 'leukopenia',
                     'pancytopenia', 'myelodysplastic', 'leukemia', 'lymphoma',
                     'myeloma', 'coagulation', 'bleeding', 'hemophilia',
                     'aplastic', 'hemolytic', 'hematological'],
        'label': 'Hematological',
        'risk_score': 2,
    },
    'other_cancer': {
        'keywords': ['gastric cancer', 'pancreatic cancer', 'hepatocellular',
                     'lung cancer', 'breast cancer', 'prostate cancer',
                     'ovarian cancer', 'bladder cancer', 'melanoma',
                     'esophageal cancer', 'thyroid cancer', 'renal cell',
                     'cervical cancer', 'endometrial cancer', 'cholangiocarcinoma',
                     'glioblastoma', 'glioma', 'sarcoma', 'neuroblastoma',
                     'carcinoma', 'neoplasm', 'malignancy', 'leukemia',
                     'lymphoma', 'myeloma'],
        'label': 'Cancer (non-CRC)',
        'risk_score': 1,  # Lower because cancer drugs often repurposed across cancer types
    },
    'crc_related': {
        'keywords': ['colorectal cancer', 'colon cancer', 'rectal cancer',
                     'colorectal neoplasm', 'colorectal carcinoma',
                     'colorectal adenoma', 'adenomatous polyp'],
        'label': 'CRC-related',
        'risk_score': 0,  # Target indication - not a safety concern
    },
    'infectious': {
        'keywords': ['infection', 'sepsis', 'pneumonia', 'tuberculosis',
                     'viral', 'bacterial', 'fungal', 'covid', 'hiv',
                     'influenza', 'hepatitis b', 'hepatitis c', 'meningitis'],
        'label': 'Infectious',
        'risk_score': 1,
    },
}

def classify_disease(disease_name):
    """Classify a disease name into safety categories. Returns list of (category_label, risk_score, cat_key)."""
    name_lower = disease_name.lower()
    matches = []
    for cat_key, cat_info in SAFETY_CATEGORIES.items():
        for kw in cat_info['keywords']:
            if kw in name_lower:
                matches.append((cat_info['label'], cat_info['risk_score'], cat_key))
                break  # One match per category

    # If no match, classify as "Other"
    if not matches:
        matches.append(('Other', 0, 'other'))

    return matches

# ── 3. Query OT Platform for gene->disease associations ───────────────────────

print("\n[2/5] Querying Open Targets Platform for gene->disease associations...")

def query_gene_diseases(gene, ensg):
    """Query all associated diseases for a gene from OT Platform."""
    query = (
        'query {'
        f'  target(ensemblId: "{ensg}") {{'
        '    id'
        '    approvedSymbol'
        '    associatedDiseases {'
        '      rows {'
        '        disease {'
        '          id'
        '          name'
        '          description'
        '          therapeuticAreas { name }'
        '        }'
        '        datasourceScores {'
        '          id'
        '          score'
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
        print(f"  ERROR {gene}: {e}")
        return None

    target = data.get('data', {}).get('target')
    if not target:
        return None

    rows = target.get('associatedDiseases', {}).get('rows') or []
    return rows

all_disease_associations = {}  # gene -> list of disease info

for g in tier1_genes:
    rows = query_gene_diseases(g['gene'], g['ensg'])
    if rows:
        # Parse each disease association
        parsed = []
        for row in rows:
            disease = row.get('disease', {})
            scores = row.get('datasourceScores', [])
            disease_name = disease.get('name', '?')
            disease_id = disease.get('id', '?')
            tas = [ta.get('name', '') for ta in disease.get('therapeuticAreas', [])]

            # Aggregate scores
            score_dict = {}
            for s in scores:
                score_dict[s.get('id', '?')] = s.get('score', 0)

            # Classify safety
            categories = classify_disease(disease_name)

            parsed.append({
                'disease_id': disease_id,
                'disease_name': disease_name,
                'therapeutic_areas': tas,
                'genetic_assoc_score': score_dict.get('genetic_association', 0),
                'gwas_score': score_dict.get('gwas_credible_sets', 0),
                'gwas_top_hit_score': score_dict.get('gwas_top_hit', 0),
                'gwas_tophit_score': score_dict.get('gwas_tophit', 0),
                'known_drug_score': score_dict.get('known_drug', 0),
                'literature_score': score_dict.get('europepmc', 0),
                'categories': categories,
                'max_risk_score': max(c[1] for c in categories),
            })

        # Sort by genetic_assoc + literature score desc
        parsed.sort(key=lambda x: x['genetic_assoc_score'] + x['literature_score'] * 0.1, reverse=True)
        all_disease_associations[g['gene']] = parsed
        n_total = len(parsed)
        n_genetic = sum(1 for p in parsed if p['genetic_assoc_score'] > 0 or p['gwas_score'] > 0)
        print(f"  {g['gene']}: {n_total} diseases ({n_genetic} with genetic evidence)")
    else:
        all_disease_associations[g['gene']] = []
        print(f"  {g['gene']}: no disease associations found")

# ── 4. Compute safety concern scores ─────────────────────────────────────────

print("\n[3/5] Computing safety concern scores...")

safety_report = []

for g in tier1_genes:
    gene = g['gene']
    diseases = all_disease_associations.get(gene, [])

    if not diseases:
        safety_report.append({
            'gene': gene, 'PPH4': g['PPH4'], 'n_cs': g['n_cs'],
            'n_diseases': 0, 'n_genetic_diseases': 0,
            'safety_concern_score': 0, 'safety_tier': 'Low',
            'top_safety_concerns': '',
            'category_breakdown': '',
            'all_disease_details': '',
            'summary': 'No data - unknown safety profile',
        })
        continue

    # Count by category
    cat_counts = defaultdict(int)
    cat_max_risk = defaultdict(int)
    all_safety_concerns = []
    genetic_concerns = []

    for d in diseases:
        for cat_label, risk_score, cat_key in d['categories']:
            cat_counts[cat_label] += 1
            cat_max_risk[cat_label] = max(cat_max_risk[cat_label], risk_score)
            if risk_score >= 2:
                has_genetic = d['genetic_assoc_score'] > 0 or d['gwas_score'] > 0
                all_safety_concerns.append({
                    'disease': d['disease_name'],
                    'category': cat_label,
                    'risk_score': risk_score,
                    'genetic_assoc': d['genetic_assoc_score'],
                    'gwas_score': d['gwas_score'],
                    'literature_score': d['literature_score'],
                    'has_genetic': has_genetic,
                })
                if has_genetic:
                    genetic_concerns.append(d)

    # Sort safety concerns: genetic first, then by risk score, then by literature
    all_safety_concerns.sort(
        key=lambda x: (-x['has_genetic'], -x['risk_score'],
                       -(x['genetic_assoc'] + x['gwas_score'] + 0.01 * x['literature_score']))
    )

    # Compute aggregate safety concern score (0-10)
    # Weighted by: number of high-risk categories + genetic evidence + literature
    n_high_risk = sum(1 for c in all_safety_concerns if c['risk_score'] >= 2)
    n_crc_only_total = cat_counts.get('CRC-related', 0)

    # Genetic evidence score: sum of genetic_assoc across non-CRC categories
    genetic_safety_score = sum(
        d['genetic_assoc_score'] + d['gwas_score']
        for d in diseases
        if d['categories'][0][0] != 'CRC-related'
    )

    # Category diversity score: number of unique risk categories
    unique_risk_cats = set(
        cat_label for d in diseases
        for cat_label, risk, _ in d['categories']
        if risk >= 2
    )

    # Composite safety score (0-10): weighted rubric with data-driven caps
    # Weights reflect evidence hierarchy: category diversity > genetic evidence > count
    # Caps are set at the 95th percentile across genes to normalize scale
    n_cats = len(unique_risk_cats)
    n_genetic = genetic_safety_score
    n_hr = n_high_risk

    safety_score = min(
        2.0 * n_cats     +  # category diversity (up to ~10 pts)
        1.5 * n_genetic  +  # genetic evidence strength (up to ~4.5 pts)
        0.5 * n_hr,         # high-risk disease count (up to ~1.5 pts)
        10.0
    )

    if safety_score >= 6:
        safety_tier = "High alert"
    elif safety_score >= 3:
        safety_tier = "Moderate caution"
    else:
        safety_tier = "Low concern"

    # Format top safety concerns (top 8)
    top_concerns_str = '; '.join([
        f"{c['disease']}"
        f"{' [GWAS]' if c['has_genetic'] else ''}"
        for c in all_safety_concerns[:8]
    ])

    # Category breakdown
    cat_bd_str = ' | '.join([
        f"{cat}={cnt}"
        for cat, cnt in sorted(cat_counts.items(), key=lambda x: -x[1])
    ])

    # Detailed disease list
    detail_parts = []
    for d in diseases:
        cat_strs = [f"{c[0]}(R{c[1]})" for c in d['categories']]
        evidence = []
        if d['genetic_assoc_score'] > 0:
            evidence.append(f"GA={d['genetic_assoc_score']:.2f}")
        if d['gwas_score'] > 0:
            evidence.append(f"GWAS={d['gwas_score']:.2f}")
        if d['literature_score'] > 0:
            evidence.append(f"Lit={d['literature_score']:.2f}")
        evidence_str = '; '.join(evidence)
        detail_parts.append(f"{d['disease_name']}[{','.join(cat_strs)}]({evidence_str})")
    all_detail_str = ' | '.join(detail_parts)

    safety_report.append({
        'gene': gene,
        'PPH4': g['PPH4'],
        'n_cs': g['n_cs'],
        'n_diseases': len(diseases),
        'n_genetic_diseases': sum(1 for d in diseases if d['genetic_assoc_score'] > 0 or d['gwas_score'] > 0),
        'n_crc_diseases': n_crc_only_total,
        'n_non_crc_high_risk': n_high_risk,
        'unique_risk_categories': len(unique_risk_cats),
        'category_names': ', '.join(sorted(unique_risk_cats)),
        'genetic_safety_score_sum': round(genetic_safety_score, 2),
        'safety_concern_score': round(safety_score, 1),
        'safety_tier': safety_tier,
        'top_safety_concerns': top_concerns_str,
        'category_breakdown': cat_bd_str,
        'all_disease_details': all_detail_str,
        'summary': '',
    })

# Sort by safety concern score desc
safety_report.sort(key=lambda x: x['safety_concern_score'], reverse=True)

# ── 5. Generate summary narratives ────────────────────────────────────────────

print("\n[4/5] Generating safety narratives...")

for sr in safety_report:
    gene = sr['gene']
    if sr['n_diseases'] == 0:
        sr['summary'] = f"{gene}: No disease association data available - unknown safety profile."
        continue

    concerns = [
        c.strip() for c in sr['top_safety_concerns'].split(';') if c.strip()
    ]
    gwas_concerns = [c for c in concerns if '[GWAS]' in c]
    lit_concerns = [c for c in concerns if '[GWAS]' not in c]

    parts = [f"{gene} is associated with {sr['n_diseases']} diseases "]
    parts.append(f"({sr['n_genetic_diseases']} with genetic evidence). ")

    if sr.get('n_crc_diseases', 0) > 0:
        parts.append(f"CRC is among its associated diseases - expected target indication. ")

    risk_cats = sr.get('category_names', '')
    if risk_cats:
        parts.append(f"Safety risk categories include: {risk_cats}. ")

    if gwas_concerns:
        parts.append(f"Genetic-evidence safety signals: {'; '.join(gwas_concerns[:3])}. ")

    if lit_concerns and not gwas_concerns:
        parts.append(f"Literature-based associations: {'; '.join(lit_concerns[:3])}. ")

    sr['summary'] = ''.join(parts)

# ── 6. Save output ───────────────────────────────────────────────────────────

print("\n[5/5] Saving results...")

# Main PheWAS summary CSV
fields = [
    'gene', 'PPH4', 'n_cs',
    'n_diseases', 'n_genetic_diseases', 'n_crc_diseases',
    'n_non_crc_high_risk', 'unique_risk_categories', 'category_names',
    'genetic_safety_score_sum',
    'safety_concern_score', 'safety_tier',
    'top_safety_concerns', 'category_breakdown',
    'summary',
]

with open("results/phase5c_phewas_safety.csv", 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(safety_report)

# Detailed disease list per gene
with open("results/phase5c_phewas_details.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['gene', 'disease_id', 'disease_name', 'therapeutic_areas',
                     'safety_category', 'risk_score',
                     'genetic_assoc_score', 'gwas_credible_set_score',
                     'known_drug_score', 'literature_score',
                     'has_genetic_evidence'])
    for g in tier1_genes:
        gene = g['gene']
        for d in all_disease_associations.get(gene, []):
            for cat_label, risk, cat_key in d['categories']:
                has_gen = d['genetic_assoc_score'] > 0 or d['gwas_score'] > 0
                writer.writerow([
                    gene, d['disease_id'], d['disease_name'],
                    '; '.join(d['therapeutic_areas']),
                    cat_label, risk,
                    d['genetic_assoc_score'], d['gwas_score'],
                    d['known_drug_score'], d['literature_score'],
                    'Yes' if has_gen else 'No',
                ])

# ── 7. Print report ──────────────────────────────────────────────────────────

print(f"\n{'='*100}")
print("  PheWAS SAFETY SCREENING - TIER 1 GENES (SuSiE Coloc-Passing)")
print(f"{'='*100}\n")

for sr in safety_report:
    tier_emoji = {'High alert': 'HIGH', 'Moderate caution': 'MODERATE', 'Low concern': 'LOW'}.get(sr['safety_tier'], '--')
    print(f"  [{tier_emoji}] {sr['gene']:<12} Safety={sr['safety_concern_score']:.1f}  "
          f"PPH4={sr['PPH4']:.3f}  "
          f"Diseases={sr['n_diseases']} (genetic={sr['n_genetic_diseases']})  "
          f"RiskCats={sr.get('category_names','-')}")
    print(f"      {sr['summary'][:200]}")
    print()

print("\n-- PheWAS Safety Score Components --\n")
print("  Score = 2.0 x category_diversity + 1.5 x genetic_evidence + 0.5 x high_risk_count")
print("  Tier: High alert (>=6) | Moderate caution (3-5.9) | Low concern (<3)")
print()

print("-- Legend --")
print("  [GWAS] = genetic evidence (GWAS fine-mapping or credible sets)")
print("  Safety categories: Autoimmune/Inflammatory(R3)  Cardiovascular(R3)")
print("                     Neurological/Psychiatric(R2)  Metabolic/Endocrine(R2)")
print("                     Renal/Hepatic(R2)  Hematological(R2)")
print("                     Cancer non-CRC(R1)  Infectious(R1)  CRC-related(R0)")
print()

print(f"Full results: results/phase5c_phewas_safety.csv")
print(f"Detailed disease list: results/phase5c_phewas_details.csv")
print("Done.")
