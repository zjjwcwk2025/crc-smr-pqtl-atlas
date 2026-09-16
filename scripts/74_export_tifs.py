#!/usr/bin/env python3
"""Export publication-grade TIFs.

Two sets are produced under submission/figures/:
  tif_figures/  one TIF per manuscript figure (Figure1-5, FigureS1-S19)
  tif_panels/   one TIF per individual panel (35 panels, incl. Fig1a)

Target 600 dpi; if a file exceeds MAXMB the dpi is stepped down (500/400/300)
so the file stays inside the usual journal upload limit.
"""
import fitz, os, glob
from PIL import Image

BASE = "/ifs1/User/zhouman/project9-v5-crc-atlas"
MAXMB = 10.0

def export(pdf, out, dpi_list=(600, 500, 400, 300)):
    d = fitz.open(pdf)
    pg = d[0]
    chosen = None
    for dpi in dpi_list:
        pix = pg.get_pixmap(dpi=dpi, colorspace=fitz.csRGB)
        im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        im.save(out, compression="tiff_lzw", dpi=(dpi, dpi))
        mb = os.path.getsize(out) / 1e6
        chosen = (dpi, pix.width, pix.height, mb)
        if mb <= MAXMB:
            break
    d.close()
    return chosen

for sub, srcs in [("tif_figures", sorted(glob.glob(BASE + "/submission/figures/Figure*.pdf"))),
                  ("tif_panels",  sorted(glob.glob(BASE + "/results/figures_rev3/*.pdf")) +
                                  [BASE + "/results/figures_rev3/Fig1a_conceptual_overview.pdf"])]:
    outdir = os.path.join(BASE, "submission/figures", sub)
    os.makedirs(outdir, exist_ok=True)
    print("== %s (%d files) ==" % (sub, len(srcs)))
    for p in srcs:
        name = os.path.basename(p)[:-4]
        out = os.path.join(outdir, name + ".tif")
        dpi, w, h, mb = export(p, out)
        flag = "" if mb <= MAXMB else "  <-- still over %.0f MB" % MAXMB
        print("   %-46s %4d dpi  %5dx%-5d  %5.2f MB%s" % (name, dpi, w, h, mb, flag))
