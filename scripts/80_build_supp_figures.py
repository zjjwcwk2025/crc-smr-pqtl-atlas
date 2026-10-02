#!/usr/bin/env python3
"""Build the Supplementary Figure S1-S22 files from the upstream panels.

Supersedes the supplementary part of 77_build_composites.py, which still knows
the pre-2026-10-02 numbering (S1-S19).  Panel order follows
manuscript/crc_smr_supplementary.tex exactly, so the stand-alone files cannot
disagree with the compiled supplement.

Rules (JTM): panel placed at NATIVE size (never scaled); canvas width <= 170 mm,
height <= 225 mm; single-panel figures copied through untouched.
Outputs PDF + 300 dpi LZW RGB TIF per figure.
"""
import os, shutil
import fitz
from PIL import Image

MM = 72.0 / 25.4
FONT = os.path.expanduser("~/work_bi/arialbd.ttf")
BAND = 3.6
GV = 2.2
GH = 5.0
LETTER_PT = 8.0
W_LIMIT = 170.0
H_LIMIT = 225.0

P = os.path.expanduser("~/project9-v5-crc-atlas/results/figures_rev3")
OUT = os.path.expanduser("~/work_bi/supp_figures_20261002")
os.makedirs(OUT, exist_ok=True)

SPEC = [
    ("FigureS1",  [["FigS1a_manhattan"], ["FigS1b_qq_plot"]], None),
    ("FigureS2",  [["FigS_add_coloc_all_pairs"]], None),
    ("FigureS3",  [["FigS2_tier1_regional_plots"]], None),
    ("FigureS4",  [["FigS3_susie_sensitivity"]], None),
    ("FigureS5",  [["FigS4_smr_or_forest"]], None),
    ("FigureS6",  [["FigS5_mr_sensitivity"]], None),
    ("FigureS7",  [["FigS6_ukbppp_pqtl_coloc"]], None),
    ("FigureS8",  [["FigS7_decode_pqtl_coloc"]], None),
    ("FigureS9",  [["FigS8_pqtl_mr_decode_forest"]], None),
    ("FigureS10", [["FigS9_pqtl_mr_ukbppp_forest"]], None),
    ("FigureS11", [["FigS10_pqtl_mr_scatter"]], None),
    ("FigureS12", [["FigS11_pqtl_mr_loo"]], None),
    ("FigureS13", [["FigS12_pqtl_mr_methods_forest"]], None),
    ("FigureS14", [["FigS13_pqtl_mr_funnel_pleiotropy"]], None),
    ("FigureS15", [["FigS14_celltype_dotplot"]], None),
    ("FigureS16", [["FigS15_spatial_tier1_genes"]], None),
    ("FigureS17", [["FigS16a_microenvironment_elbow_plot", "FigS16b_me_zone_barplot"],
                   ["FigS16c_me_tier1_heatmap"]], "ABC"),
    ("FigureS18", [["FigS_add_gtex_scatter"]], None),
    ("FigureS19", [["FigS17_finngen_replication_forest"]], None),
    ("FigureS20", [["FigS18_idg_druggability"]], None),
    ("FigureS21", [["FigS19_phewas_safety_scores"]], None),
    ("FigureS22", [["FigS_add_drug_repurposing"]], None),
]


def panel(p):
    return os.path.join(P, p + ".pdf")


def size(p):
    d = fitz.open(panel(p)); r = d[0].rect; d.close()
    return r.width / MM, r.height / MM


def tif(pdf, out, dpi=300):
    d = fitz.open(pdf)
    pix = d[0].get_pixmap(dpi=dpi, colorspace=fitz.csRGB)
    im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    im.save(out, compression="tiff_lzw", dpi=(dpi, dpi))
    d.close()
    return pix.width, pix.height


missing, warn = [], []
for name, rows, letters in SPEC:
    n_panels = sum(len(r) for r in rows)
    for r in rows:
        for p in r:
            if not os.path.isfile(panel(p)):
                missing.append(p)
    if missing:
        continue
    out_pdf = os.path.join(OUT, name + ".pdf")
    if n_panels == 1:
        shutil.copy2(panel(rows[0][0]), out_pdf)
        w, h = size(rows[0][0])
    else:
        pw = [[size(p) for p in row] for row in rows]
        row_w = [sum(x for x, y in r) + GH * (len(r) - 1) for r in pw]
        row_h = [max(y for x, y in r) for r in pw]
        W = max(row_w)
        H = sum(row_h) + BAND * len(rows) + GV * (len(rows) - 1)
        doc = fitz.open()
        page = doc.new_page(width=W * MM, height=H * MM)
        y, k = 0.0, 0
        for i, row in enumerate(rows):
            x = (W - row_w[i]) / 2.0
            top = y + BAND
            for j, p in enumerate(row):
                px, py = pw[i][j]
                rect = fitz.Rect(x * MM, top * MM, (x + px) * MM, (top + py) * MM)
                src = fitz.open(panel(p))
                page.show_pdf_page(rect, src, 0)
                src.close()
                if letters:
                    page.insert_text(fitz.Point(x * MM + 0.5, (top - 1.0) * MM),
                                     letters[k], fontsize=LETTER_PT,
                                     fontname="ArialBd", fontfile=FONT, color=(0, 0, 0))
                    k += 1
                x += px + GH
            y = top + row_h[i] + GV
        doc.save(out_pdf, garbage=4, deflate=True)
        doc.close()
        w, h = W, H
    flag = "ok  "
    if w > W_LIMIT + 0.5 or h > H_LIMIT:
        flag = "WARN"
        warn.append(name)
    tw, th = tif(out_pdf, os.path.join(OUT, name + ".tif"))
    print("%s %-11s %6.1f x %6.1f mm  panels=%d  tif %dx%d  %.2f MB"
          % (flag, name, w, h, n_panels, tw, th,
             os.path.getsize(os.path.join(OUT, name + ".tif")) / 1e6))

if missing:
    print("MISSING PANELS:", sorted(set(missing)))
print("WARNED:", warn if warn else "none")
print("SUPPDONE")
