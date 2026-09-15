#!/usr/bin/env python3
"""rev5: trim supplementary figures 27 -> 19 and repair in-text references."""
import re, sys

TEX = "manuscript/crc_smr_atlas_rev3.tex"
src = open(TEX, encoding="utf-8").read().replace("\r\n", "\n")
fails = []

GUARD = r"(?:(?!\\begin\{figure\}).)*?"

def env_re(label):
    return re.compile(r"\\begin\{figure\}\[htbp\]" + GUARD + r"\\label\{fig:" + label + r"\}" + GUARD + r"\\end\{figure\}\n+", re.S)

def drop(label):
    global src
    pat = env_re(label)
    c = len(pat.findall(src))
    if c != 1:
        fails.append("DROP %s: found %d" % (label, c)); return
    src = pat.sub("", src, count=1)
    print("ok  drop figure %s" % label)

def swap(label, new):
    global src
    pat = env_re(label)
    c = len(pat.findall(src))
    if c != 1:
        fails.append("SWAP %s: found %d" % (label, c)); return
    src = pat.sub(lambda m: new, src, count=1)
    print("ok  merge into figure %s" % label)

def rep(pattern, new, label):
    global src
    pat = re.compile(pattern, re.S)
    c = len(pat.findall(src))
    if c != 1:
        fails.append("REP %s: found %d :: %s" % (label, c, pattern[:70])); return
    src = pat.sub(lambda m: new, src, count=1)
    print("ok  %s" % label)

# ---- 1. drop the five figures -------------------------------------------------
for lab in ["S10", "S11", "S12", "S13", "S17"]:
    drop(lab)

# ---- 2. merge S8 + S9 ---------------------------------------------------------
merge89 = r"""\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
\centering
\includegraphics[width=\linewidth]{FigS8_spatial_TMBIM1.pdf}
\caption{}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
\centering
\includegraphics[width=\linewidth]{FigS9_spatial_meCAF_TMBIM1.pdf}
\caption{}
\end{subfigure}
\caption{\textbf{Spatial expression of TMBIM1 and colocalization with myofibroblastic CAF markers in CRC tissue.}
(a) Visium spatial expression of TMBIM1 across $n = 3{,}493$ spots; TMBIM1 was detected in 57\% of
spots with fibroblast-enriched localization (colour scale: SCT-normalized expression).
(b) Spatial colocalization of the myofibroblastic CAF (meCAF) marker-gene signature score with
TMBIM1 expression (Spearman $\rho = 0.091$, $p = 6.5 \times 10^{-8}$, $n = 3{,}493$ spots;
colour scale: SCT-normalized meCAF signature score).}
\label{fig:S8}
\end{figure}
"""
drop("S9")
swap("S8", merge89)

# ---- 3. merge S14 + S15 + S16 -------------------------------------------------
merge1416 = r"""\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.46\textwidth}
\centering
\includegraphics[width=\linewidth]{FigS14_microenvironment_elbow_plot.pdf}
\caption{}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.46\textwidth}
\centering
\includegraphics[width=\linewidth]{FigS15_me_zone_barplot.pdf}
\caption{}
\end{subfigure}
\\[0.6em]
\begin{subfigure}[b]{0.92\textwidth}
\centering
\includegraphics[width=\linewidth]{FigS16_me_tier1_heatmap.pdf}
\caption{}
\end{subfigure}
\caption{\textbf{Spatial microenvironment zoning of CRC Visium data.}
(a) K-means within-cluster sum of squares for $k = 2$ to $k = 8$; the elbow at $k = 6$ defined
the final partition. (b) Composition of the six K-means-defined microenvironment (ME) zones
across 3{,}493 Visium spots; ME6 Stromal matrix was the largest zone (927 spots, 26.5\%),
followed by ME3 Mast cell zone (746, 21.4\%), ME5 T cell territory (613, 17.5\%), ME1
Myeloid-enriched (580, 16.6\%), ME2 Epithelial-dominant (518, 14.8\%) and ME4 B cell niches
(109, 3.1\%). (c) Expression of the 34 coloc-passing genes detected in Visium data across the
six ME zones; BMP2, LAMC1 and PNKD were preferentially expressed in the stromal ME6 zone and
GPBAR1 was detected in a single epithelial ME2 spot (0.03\% of 3{,}493 spots).}
\label{fig:S14}
\end{figure}
"""
drop("S15")
drop("S16")
swap("S14", merge1416)

# ---- 4. repair in-text references --------------------------------------------
rep(r"Supplementary Fig\.~\\ref\{fig:S8\}--S9\)",
    "Supplementary Fig.~\\\\ref{fig:S8})", "text: merged S8")
rep(r"Supplementary Fig\.~\\ref\{fig:S14\}--S16;",
    "Supplementary Fig.~\\\\ref{fig:S14};", "text: merged S14")
rep(r"and pairwise spatial\s+cross-correlation across the chr2:219Mb locus revealed only weak,\s+distance-decaying interdependencies \(Supplementary Fig\.~\\ref\{fig:S17\}; Supplementary Table S9\)\.",
    "Pairwise spatial cross-correlation statistics for the chr2:219Mb locus are provided in\nSupplementary Table S9.",
    "text: drop S17 clause")
rep(r"[,;]\s*Supplementary Figs\.~\\ref\{fig:S10\}--\\ref\{fig:S13\}\)",
    ")", "text: drop S10-S13 ref")

if fails:
    print("\n=== FAILURES ===")
    for f in fails:
        print(" -", f)
    sys.exit(1)

open(TEX, "w", encoding="utf-8").write(src)
print("\nwrote %s" % TEX)
print("supplementary figure environments remaining: %d" % len(re.findall(r"\\label\{fig:S\d+\}", src)))