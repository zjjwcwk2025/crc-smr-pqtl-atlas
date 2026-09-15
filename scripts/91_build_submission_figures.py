import re, os, subprocess, sys, shutil
ONLY = [x for x in os.environ.get('ONLY','').split(',') if x]
import fitz
from PIL import Image, ImageOps

BASE = '/ifs1/User/zhouman/project9-v5-crc-atlas'
OUT  = f'{BASE}/submission/figures'
SRC  = f'{OUT}/src'
os.makedirs(SRC, exist_ok=True)
os.chdir(BASE)

def readdoc(p): return open(p, encoding='utf-8').read()

def preamble(path):
    s = readdoc(path)
    i = s.index('\\begin{document}')
    pre = s[:i]
    pre = re.sub(r'\\externaldocument\{[^}]*\}\n?', '', pre)
    pre = pre.replace('\\usepackage[margin=2.5cm]{geometry}',
                      '\\usepackage[paperwidth=6.69in,paperheight=64in,margin=0.5cm]{geometry}')
    pre = pre.replace('\\usepackage{graphicx}', '\\usepackage{graphicx}\n\\usepackage{float}')
    return pre

def spans(s, cmd=r'\\caption'):
    out = []
    for m in re.finditer(cmd + r'(?![a-zA-Z])(?:\s*\[[^\]]*\])?\s*\{', s):
        b = m.end()-1; d=0; j=b; esc=False
        while j < len(s):
            c = s[j]
            if esc: esc=False
            elif c=='\\': esc=True
            elif c=='{': d+=1
            elif c=='}':
                d-=1
                if d==0: break
            j+=1
        out.append((m.start(), j+1))
    return out

def blank_subcaptions(blk):
    """inside subfigure envs: \caption{...} -> \caption{}"""
    res = []
    for m in re.finditer(r'\\begin\{subfigure\}.*?\\end\{subfigure\}', blk, re.S):
        res.append((m.start(), m.end()))
    out = []
    pos = 0
    for a, b in res:
        out.append(blk[pos:a])
        sub = blk[a:b]
        for s0, s1 in reversed(spans(sub)):
            sub = sub[:s0] + '@@SUBCAP@@' + sub[s1:]
        out.append(sub); pos = b
    out.append(blk[pos:])
    return ''.join(out)

def drop_outer_captions(blk):
    for s0, s1 in reversed(spans(blk)):
        blk = blk[:s0] + blk[s1:]
    return blk

def clean_block(blk):
    blk = blank_subcaptions(blk)
    blk = drop_outer_captions(blk)
    blk = blk.replace('@@SUBCAP@@', '\\caption{}')
    pass  # keep \ContinuedFloat: it stops \caption from resetting the subfigure counter
    blk = re.sub('\\\\begin\\{figure\\}\\[[^]]*\\]', '\\\\begin{figure}[H]', blk)
    blk = blk.replace('\\end{figure}', '\\end{figure}')
    return blk

def standalone(pre, body):
    return (pre + '\\begin{document}\n\\thispagestyle{empty}\n\\pagestyle{empty}\n'
            + body + '\n\\end{document}\n')

def build(name, pre, blocks):
    body = '\n\n\\vspace{1.2em}\n\n'.join(clean_block(b) for b in blocks)
    tex = standalone(pre, body)
    src = f'{SRC}/{name}.tex'
    open(src, 'w', encoding='utf-8').write(tex)
    for _ in range(1):
        r = subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',
                            f'-output-directory={SRC}', src],
                           capture_output=True, text=True)
    if r.returncode != 0:
        print(f'  !! pdflatex FAILED for {name}')
        print('\n'.join(r.stdout.split('\n')[-25:]))
        return None
    pdf = f'{SRC}/{name}.pdf'
    return pdf if os.path.exists(pdf) else None

def crop(pdf, out):
    d = fitz.open(pdf)
    page = d[0]
    pix = page.get_pixmap(dpi=100)
    img = Image.frombytes('RGB', (pix.width, pix.height), pix.samples).convert('L')
    bb = ImageOps.invert(img).getbbox()
    if not bb:
        d.close(); return None
    sc = 72.0/100.0
    pad = 8
    r = fitz.Rect(max(0,(bb[0]-pad)*sc), max(0,(bb[1]-pad)*sc),
                  min(page.rect.x1,(bb[2]+pad)*sc), min(page.rect.y1,(bb[3]+pad)*sc))
    nd = fitz.open()
    np_ = nd.new_page(width=r.width, height=r.height)
    np_.show_pdf_page(np_.rect, d, 0, clip=r)
    nd.save(out, deflate=True, garbage=3)
    nd.close(); d.close()
    return (round(r.width/72,2), round(r.height/72,2))

def tiff(pdf, out, dpi=600):
    d = fitz.open(pdf)
    pix = d[0].get_pixmap(dpi=dpi, colorspace=fitz.csRGB)
    im = Image.frombytes('RGB', (pix.width, pix.height), pix.samples)
    im.save(out, compression='tiff_lzw', dpi=(dpi, dpi))
    d.close()
    return f'{pix.width}x{pix.height}'

# ---------- MAIN ----------
pre_main = preamble('manuscript/crc_smr_atlas_rev3.tex')
s = re.sub(r'(?<!\\)%.*', '', readdoc('manuscript/crc_smr_atlas_rev3.tex'))
envs = [m.group(0) for m in re.finditer(r'\\begin\{figure\}.*?\\end\{figure\}', s, re.S)]
groups, cur = [], None
for e in envs:
    if '\\ContinuedFloat' in e:
        cur.append(e)
    else:
        cur = [e]; groups.append(cur)
print(f'MAIN: {len(groups)} figure groups')

report = []
for i, g in enumerate(groups, 1):
    name = f'Figure{i}'
    if ONLY and name not in ONLY: continue
    pdf = build(name, pre_main, g)
    if not pdf: continue
    dst = f'{OUT}/{name}.pdf'
    dim = crop(pdf, dst)
    px = tiff(dst, f'{OUT}/{name}.tif')
    sz = os.path.getsize(f'{OUT}/{name}.tif')/1e6
    report.append((name, dim, px, round(sz,1), [re.search(r'[^/]+\.pdf', x).group(0) for gg in g for x in re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', gg)]))
    print(f'  {name}: {dim} in, tif {px}, {sz:.1f} MB')

# ---------- SUPP ----------
pre_supp = preamble('manuscript/crc_smr_supplementary.tex')
t = re.sub(r'(?<!\\)%.*', '', readdoc('manuscript/crc_smr_supplementary.tex'))
senvs = [m.group(0) for m in re.finditer(r'\\begin\{figure\}.*?\\end\{figure\}', t, re.S)]
print(f'SUPP: {len(senvs)} figures')
for i, e in enumerate(senvs, 1):
    name = f'FigureS{i}'
    if ONLY and name not in ONLY: continue
    pdf = build(name, pre_supp, [e])
    if not pdf: continue
    dst = f'{OUT}/{name}.pdf'
    dim = crop(pdf, dst)
    px = tiff(dst, f'{OUT}/{name}.tif')
    sz = os.path.getsize(f'{OUT}/{name}.tif')/1e6
    report.append((name, dim, px, round(sz,1), re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', e)))
    print(f'  {name}: {dim} in, tif {px}, {sz:.1f} MB')

print()
print('=== REPORT ===')
for r in report:
    print(f'  {r[0]:11s} {str(r[1]):16s} {r[2]:>12s} {r[3]:>7.1f} MB  panels={r[4]}')
open(f'{OUT}/_build_report.txt','w').write('\n'.join(str(r) for r in report))
print('DONE-FIGBUILD')
