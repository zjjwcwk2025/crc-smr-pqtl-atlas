#!/usr/bin/env python3
"""Build submission composites for the CRC atlas (JTM).

Source of truth: results/figures_rev3/ (36 panels).
Rules (2026-10-01 rewrite):
  * every panel placed at NATIVE size - never scaled
  * canvas width <= 170 mm, height <= 225 mm
  * panel letters: Arial Bold 8 pt, UPPERCASE, WITH parentheses, e.g. (A),
    at the top-left of the panel band
  * outputs per figure: PDF (vector master) + TIF 300 dpi RGB + PNG 300 dpi
  * writes a manifest: canvas mm, per-panel placement, SHA256
"""
import os, json, hashlib
import fitz
from PIL import Image

MM = 72.0 / 25.4
FONT = os.path.expanduser("~/work_bi/arialbd.ttf")
BAND = 3.6
GV = 2.2
GH = 3.5
LETTER_PT = 8.0
W_LIMIT = 170.0
H_LIMIT = 225.0

P = os.path.expanduser("~/project9-v5-crc-atlas/results/figures_rev8_cns")
OUT = os.path.expanduser("~/work_bi/fig_composites_20261002b")
os.makedirs(OUT, exist_ok=True)

SPEC = [
    ("Figure1", [["Fig1a_conceptual_overview"], ["Fig1b_study_design"]], "AB"),
    ("Figure2", [["Fig2a_coverage_funnel"],
                 ["Fig2b_eqtl_pqtl_classification", "Fig2c_eqtl_pqtl_scatter"]], "ABC"),
    ("Figure3", [["Fig3a_ccm2_convergence"], ["Fig3c_gtex_colon_forest"]], "AB"),
    ("Figure4", [["Fig4a_idg_druggability"], ["Fig4b_phewas_safety"]], "AB"),
    ("Figure5", [["Fig5a_competitor_overlap"], ["Fig5b_gene_overlap_upset"]], "AB"),
]


def panel(p):
    return os.path.join(P, p + ".pdf")


def size(p):
    d = fitz.open(panel(p)); r = d[0].rect; d.close()
    return r.width / MM, r.height / MM


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


manifest = []
warn = []
for name, rows, letters in SPEC:
    out_pdf = os.path.join(OUT, name + ".pdf")
    pw = [[size(p) for p in row] for row in rows]
    row_w = [sum(x for x, y in r) + GH * (len(r) - 1) for r in pw]
    row_h = [max(y for x, y in r) for r in pw]
    W = max(row_w)
    H = sum(row_h) + BAND * len(rows) + GV * (len(rows) - 1)
    doc = fitz.open()
    page = doc.new_page(width=W * MM, height=H * MM)
    placements = []
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
                page.insert_text(fitz.Point((x + 0.5) * MM, (top - 1.0) * MM),
                                 "(" + letters[k] + ")", fontsize=LETTER_PT,
                                 fontname="ArialBd", fontfile=FONT, color=(0, 0, 0))
                k += 1
            placements.append({"panel": p, "x_mm": round(x, 2), "y_mm": round(top, 2),
                               "w_mm": round(px, 2), "h_mm": round(py, 2),
                               "letter": ("(" + letters[k - 1] + ")") if letters else None})
            x += px + GH
        y = top + row_h[i] + GV
    doc.save(out_pdf, garbage=4, deflate=True)
    doc.close()

    d = fitz.open(out_pdf)
    pix = d[0].get_pixmap(dpi=300, colorspace=fitz.csRGB)
    im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    out_tif = os.path.join(OUT, name + ".tif")
    im.save(out_tif, compression="tiff_lzw", dpi=(300, 300))
    out_png = os.path.join(OUT, name + ".png")
    d[0].get_pixmap(dpi=300, colorspace=fitz.csRGB).save(out_png)
    d.close()

    flag = "ok  "
    if W > W_LIMIT + 0.5 or H > H_LIMIT:
        flag = "WARN"
        warn.append(name)
    manifest.append({
        "figure": name, "canvas_mm": [round(W, 2), round(H, 2)],
        "aspect_h_over_w": round(H / W, 4), "n_panels": sum(len(r) for r in rows),
        "placements": placements,
        "sha256": {"pdf": sha256(out_pdf), "tif": sha256(out_tif), "png": sha256(out_png)},
        "bytes": {k2: os.path.getsize(os.path.join(OUT, name + ext))
                  for k2, ext in (("pdf", ".pdf"), ("tif", ".tif"), ("png", ".png"))},
    })
    print("%s %-8s %6.1f x %6.1f mm  panels=%d  ratio=%.3f" %
          (flag, name, W, H, sum(len(r) for r in rows), H / W))

with open(os.path.join(OUT, "manifest.json"), "w", encoding="utf-8") as fh:
    json.dump({"canvas_limit_mm": [W_LIMIT, H_LIMIT], "figures": manifest,
               "sha256": {"Figure1.pdf": manifest[0]["sha256"]["pdf"]}},
              fh, ensure_ascii=False, indent=2)
print("WARNED:", warn if warn else "none")
print("COMPOSITESDONE")
