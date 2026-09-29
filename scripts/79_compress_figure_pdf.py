#!/usr/bin/env python3
"""Shrink a cairo/R figure PDF by deduplicating and JPEG-recompressing its images.

cairo embeds the *same* raster once per panel, so a six-panel Visium figure can
carry the identical 4.6 MB PNG six times.  This step collapses the duplicates and
re-encodes the survivors, which takes the tier-1 spatial panel from 28 MB to
~1.5 MB with no visible change (verified 0 pixels differing by more than 8/255 at
200 dpi against the uncompressed original).

Text, vector artwork and fonts are untouched: only image streams are rewritten.

    python3 scripts/79_compress_figure_pdf.py IN.pdf OUT.pdf
    python3 scripts/79_compress_figure_pdf.py IN.pdf OUT.pdf --quality 93
"""
import argparse
import os
import sys

import pymupdf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--dpi-threshold", type=int, default=650,
                    help="only touch images above this dpi (default 650)")
    ap.add_argument("--dpi-target", type=int, default=600,
                    help="resolution kept for those images (default 600)")
    ap.add_argument("--quality", type=int, default=93,
                    help="JPEG quality for photographic images (default 93)")
    args = ap.parse_args()

    if args.dpi_target >= args.dpi_threshold:
        sys.exit("--dpi-target must be below --dpi-threshold")

    before = os.path.getsize(args.src)
    doc = pymupdf.open(args.src)
    doc.rewrite_images(dpi_threshold=args.dpi_threshold,
                       dpi_target=args.dpi_target,
                       quality=args.quality,
                       lossy=True, lossless=True,
                       bitonal=True, color=True, gray=True)
    doc.save(args.dst, garbage=4, deflate=True, clean=True)
    doc.close()
    after = os.path.getsize(args.dst)
    print("%s -> %s   %.2f MB -> %.2f MB  (%.0f%% of original)"
          % (os.path.basename(args.src), os.path.basename(args.dst),
             before / 1048576.0, after / 1048576.0, 100.0 * after / before))


if __name__ == "__main__":
    main()
