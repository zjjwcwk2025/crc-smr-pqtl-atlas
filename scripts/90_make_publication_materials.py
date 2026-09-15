#!/usr/bin/env python3
"""Generate publication-grade figure PNGs (300 DPI) and copy table .tex files.

Deliverables (under project root):
  publication/figures/        -> 30 authoritative vector PDFs (current naming)
  publication/figures_png/    -> 300-DPI PNG for each figure
  publication/tables/         -> 10 authoritative table .tex files
  publication/tables_pdf/     -> each table compiled to standalone PDF
"""
import os
import shutil
import subprocess
import sys

ROOT = "/ifs1/User/zhouman/project9-v5-crc-atlas"
FIG_SRC = os.path.join(ROOT, "results", "figures")
FIG_SCRN_SRC = os.path.join(ROOT, "results", "phase6_scrna")
TAB_SRC = os.path.join(ROOT, "manuscript", "tables")
FIG_DST = os.path.join(ROOT, "publication", "figures")
PNG_DST = os.path.join(ROOT, "publication", "figures_png")
TAB_DST = os.path.join(ROOT, "publication", "tables")
TABPDF_DST = os.path.join(ROOT, "publication", "tables_pdf")

for d in (FIG_DST, PNG_DST, TAB_DST, TABPDF_DST):
    os.makedirs(d, exist_ok=True)

# --- 1. Collect 30 authoritative figure PDFs ------------------------------
figures = []
for f in sorted(os.listdir(FIG_SRC)):
    if f.endswith(".pdf"):
        figures.append(os.path.join(FIG_SRC, f))
# two supplementary figures live in phase6_scrna
for f in ("FigS14_microenvironment_elbow_plot.pdf",
          "FigS17_spatial_chr2_crosscorrelogram.pdf"):
    figures.append(os.path.join(FIG_SCRN_SRC, f))

figures = sorted(set(figures))

# --- 2. Refresh publication/figures (clear stale, copy authoritative) ------
for old in os.listdir(FIG_DST):
    if old.endswith(".pdf"):
        os.remove(os.path.join(FIG_DST, old))

copied = 0
for src in figures:
    dst = os.path.join(FIG_DST, os.path.basename(src))
    shutil.copy2(src, dst)
    copied += 1

# --- 3. Generate 300-DPI PNG ----------------------------------------------
png_ok = 0
png_fail = []
for src in figures:
    base = os.path.splitext(os.path.basename(src))[0]
    out = os.path.join(PNG_DST, base + ".png")
    r = subprocess.run(
        ["pdftoppm", "-singlefile", "-png", "-r", "300", src, out[:-4]],
        capture_output=True, text=True)
    if r.returncode == 0 and os.path.exists(out):
        png_ok += 1
    else:
        png_fail.append((os.path.basename(src), r.stderr.strip()[:200]))

# --- 4. Refresh tables (.tex) ---------------------------------------------
tab_copied = 0
tables = sorted(f for f in os.listdir(TAB_SRC) if f.endswith(".tex"))
for t in tables:
    shutil.copy2(os.path.join(TAB_SRC, t), os.path.join(TAB_DST, t))
    tab_copied += 1

# --- 5. Compile each table to standalone PDF ------------------------------
PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{longtable}
\usepackage[margin=1in]{geometry}
\pagestyle{empty}
\begin{document}
\INPUT
\end{document}
"""

tabpdf_ok = 0
tabpdf_fail = []
for t in tables:
    base = os.path.splitext(t)[0]
    wrapper = os.path.join(ROOT, "publication", f"_tmp_{base}.tex")
    with open(wrapper, "w") as fh:
        fh.write(PREAMBLE.replace("\\INPUT", "\\input{%s}" %
                                  os.path.join(TAB_DST, t)))
    r = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode",
         "-output-directory", ROOT + "/publication",
         wrapper],
        capture_output=True, text=True)
    # pdflatex writes the pdf next to the wrapper (same base name)
    produced = os.path.join(ROOT, "publication", f"_tmp_{base}.pdf")
    target = os.path.join(TABPDF_DST, base + ".pdf")
    if os.path.exists(produced):
        shutil.move(produced, target)
        tabpdf_ok += 1
    else:
        tail = r.stderr.strip().splitlines()[-5:] if r.stderr else []
        tabpdf_fail.append((t, "\n".join(tail)[:300]))
    # cleanup aux/log/tex
    for ext in (".aux", ".log", ".tex"):
        p = os.path.join(ROOT, "publication", f"_tmp_{base}{ext}")
        if os.path.exists(p):
            os.remove(p)

print("=== PUBLICATION MATERIALS SUMMARY ===")
print(f"figures PDF copied : {copied}/30")
print(f"figures PNG (300dpi): {png_ok}/30")
if png_fail:
    print("  PNG FAIL:", png_fail)
print(f"tables .tex copied : {tab_copied}/10")
print(f"tables PDF compiled: {tabpdf_ok}/10")
if tabpdf_fail:
    print("  TABLE PDF FAIL:", tabpdf_fail)
