#!/usr/bin/env python3
"""Make results/figures_rev3 the ONE source of truth for panel files.

The manuscript tex and the supplementary tex both put results/figures_rev3/
first on \graphicspath, but that folder still held the 2026-09-16 vintage of
the main-figure panels while the submission composites were built from the
newer 2026-09-18/19 re-draws that live only in work_bi/ai_user/{rev5,rev6}_safe.
Result: the compiled manuscript PDF showed a different version of Fig. 2-5
than the figure files that get uploaded.

This script copies the current panels into rev3 (originals backed up) and
copies the new compact Fig 1c back into rev6_safe so that folder can no longer
disagree either.
"""
import os, shutil, datetime
import fitz

MM = 72.0 / 25.4
HOME = os.path.expanduser('~')
R3 = os.path.join(HOME, 'project9-v5-crc-atlas/results/figures_rev3')
AI = os.path.join(HOME, 'work_bi/ai_user')
BK = os.path.join(HOME, 'work_bi/rev3_backup_20260926')
os.makedirs(BK, exist_ok=True)

SRC = {
    'Fig1a_conceptual_overview': 'rev6_safe',
    'Fig1b_study_design': 'rev6_safe',
    'Fig2a_coverage_funnel': 'rev6_safe',
    'Fig2b_eqtl_pqtl_classification': 'rev6_safe',
    'Fig2c_eqtl_pqtl_scatter': 'rev6_safe',
    'Fig3a_ccm2_convergence': 'rev5_safe',
    'Fig3b_gtex_colon_scatter': 'rev5_safe',
    'Fig3c_gtex_colon_forest': 'rev5_safe',
    'Fig4a_phewas_safety': 'rev5_safe',
    'Fig4b_drug_repurposing': 'rev5_safe',
    'Fig5a_competitor_overlap': 'rev5_safe',
    'Fig5b_gene_overlap_venn': 'rev5_safe',
}


def dim(p):
    d = fitz.open(p); r = d[0].rect; d.close()
    return r.width / MM, r.height / MM


for name, sub in sorted(SRC.items()):
    src = os.path.join(AI, sub, name + '.pdf')
    dst = os.path.join(R3, name + '.pdf')
    if not os.path.exists(src):
        print('MISSING SOURCE', src)
        continue
    if os.path.exists(dst):
        shutil.copy2(dst, os.path.join(BK, name + '.pdf'))
    shutil.copy2(src, dst)
    w, h = dim(dst)
    print('%-34s <- %-10s %.1f x %.1f mm' % (name, sub, w, h))

# Fig 1c: push the new compact panel back into rev6_safe as well
c_src = os.path.join(R3, 'Fig1c_coloc_pph4.pdf')
c_dst = os.path.join(AI, 'rev6_safe/Fig1c_coloc_pph4.pdf')
if os.path.exists(c_dst):
    shutil.copy2(c_dst, os.path.join(BK, 'Fig1c_coloc_pph4.rev6_old.pdf'))
shutil.copy2(c_src, c_dst)
w, h = dim(c_dst)
print('%-34s -> %-10s %.1f x %.1f mm' % ('Fig1c_coloc_pph4', 'rev6_safe', w, h))
print('SYNCPANELSDONE')