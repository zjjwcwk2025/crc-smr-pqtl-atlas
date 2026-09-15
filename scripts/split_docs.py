#!/usr/bin/env python3
"""rev6: split the manuscript into main text + standalone supplementary file (xr cross-refs)."""
import re, sys

TEX  = "manuscript/crc_smr_atlas_rev3.tex"
SUPP = "manuscript/crc_smr_supplementary.tex"

src = open(TEX, encoding="utf-8").read().replace("\r\n", "\n")
marker = "%------------------------------------------------------------------\n\\setcounter{figure}{0}\n"
if src.count(marker) != 1:
    print("marker count =", src.count(marker)); sys.exit(1)
i = src.index(marker)
main_body, supp_block = src[:i], src[i:]

supp_block = supp_block.rstrip()
tail = "\\end{document}"
if not supp_block.endswith(tail):
    print("unexpected end of file"); sys.exit(1)
supp_block = supp_block[:-len(tail)].rstrip() + "\n"

# ---- main document: add xr, drop the supplementary block --------------------
if "\\usepackage{xr}" not in main_body:
    main_body = main_body.replace(
        "\\usepackage{subcaption}",
        "\\usepackage{subcaption}\n"
        "\\usepackage{xr}\n"
        "\\externaldocument{crc_smr_supplementary}", 1)
main_out = main_body.rstrip() + "\n\n\\end{document}\n"
open(TEX, "w", encoding="utf-8").write(main_out)

# ---- supplementary document ------------------------------------------------
preamble = r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[margin=2.5cm]{geometry}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{hyperref}
\usepackage{natbib}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{subcaption}
\captionsetup[subfigure]{labelformat=parens,labelsep=space}
\renewcommand{\thesubfigure}{\alph{subfigure}}
\graphicspath{{results/figures_rev3/}{results/figures/}{results/phase6_scrna/}}

\title{Supplementary Material\\[0.4em]
\large Limited plasma proteomic coverage constrains genetic drug-target prioritization in colorectal cancer}
\author{Jiajie Zhou}
\date{}

\begin{document}
\maketitle

"""
bib = ""
if "\\cite" in supp_block:
    bib = "\n\\bibliographystyle{plainnat}\n\\bibliography{manuscript/crc_smr_atlas}\n"
    print("note: supplementary block cites references; bibliography included")
supp_out = preamble + supp_block.rstrip() + "\n" + bib + "\n\\end{document}\n"
open(SUPP, "w", encoding="utf-8").write(supp_out)

print("main tex  : %d bytes" % len(main_out))
print("supp tex  : %d bytes" % len(supp_out))
print("main figures (main doc)  :", len(re.findall(r"\\label\{fig:(?!S)\w+\}", main_out)))
print("supp figures in supp doc :", len(re.findall(r"\\label\{fig:S\d+\}", supp_out)))