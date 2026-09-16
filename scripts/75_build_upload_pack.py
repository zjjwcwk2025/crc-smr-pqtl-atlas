#!/usr/bin/env python3
"""Build a flat, self-contained upload pack for journal submission.

Journals that compile LaTeX on their side (Springer Nature / BMC titles such as the
Journal of Translational Medicine) place every uploaded file in a single directory.
The repository layout uses `manuscript/...` paths, which break once the tree is
flattened.  This script rewrites those paths, copies the figures and tables next to
the .tex files, and drops both the .bib and the already compiled .bbl in place.

Output
------
submission/upload_pack/main_text/      main text + 3 tables + 12 figure panels + .bib/.bbl
submission/upload_pack/supplementary/  supplement + 9 tables + 22 figure panels

Both directories are self-contained and compile with two plain pdflatex runs.
Needs the figure sources and the compiled .bbl to be present; if anything is missing the
script prints the list and exits non-zero instead of writing a half pack.

    python3 scripts/75_build_upload_pack.py
"""
import os
import re
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

FIGURE_DIRS = ['results/figures_rev3', 'results/figures', 'results/phase6_scrna',
               'submission/figures']
OUT = 'submission/upload_pack'


def locate(name):
    for d in FIGURE_DIRS:
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return None


def build(tex, sub):
    dst = os.path.join(OUT, sub)
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.makedirs(dst)

    job = os.path.splitext(os.path.basename(tex))[0]
    src = open(tex, encoding='utf-8').read()
    flat = re.sub(r'\\input\{manuscript/tables/', r'\\input{', src)
    flat = re.sub(r'\\bibliography\{manuscript/', r'\\bibliography{', flat)
    flat = re.sub(r'\\graphicspath\{[^\n]*\}', r'\\graphicspath{{./}{figures/}{images/}}', flat)
    open(os.path.join(dst, os.path.basename(tex)), 'w', encoding='utf-8').write(flat)

    copied, missing = 0, []

    def put(source, name=None):
        nonlocal copied
        if os.path.isfile(source):
            shutil.copy2(source, os.path.join(dst, name or os.path.basename(source)))
            copied += 1
            return True
        return False

    for m in re.findall(r'\\input\{([^}]*)\}', flat):
        n = os.path.basename(m)
        if not put(os.path.join('manuscript', 'tables', n)):
            missing.append('table ' + n)

    for m in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}', flat):
        n = os.path.basename(m)
        if os.path.isfile(os.path.join(dst, n)):
            continue
        p = locate(n)
        if p:
            put(p, n)
        else:
            missing.append('figure ' + n)

    for b in re.findall(r'\\bibliography\{([^}]*)\}', flat):
        base = os.path.basename(b)
        if not put(os.path.join('manuscript', base + '.bib')):
            put(base + '.bib')
        # BibTeX names its output after the job, not after the .bib file.
        if not put(job + '.bbl'):
            missing.append('bbl ' + job + '.bbl')

    print('%-14s %2d files copied, %d missing%s'
          % (sub, copied, len(missing), ' ' + ', '.join(missing) if missing else ''))
    return missing


def main():
    missing = build('manuscript/crc_smr_atlas_rev3.tex', 'main_text')
    missing += build('manuscript/crc_smr_supplementary.tex', 'supplementary')
    if missing:
        sys.exit('refusing to finish: %d input file(s) not found. Build the figures and '
                 'the .bbl first.' % len(missing))
    print('upload pack written to ' + OUT)


if __name__ == '__main__':
    main()
