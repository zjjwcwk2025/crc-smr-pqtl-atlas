import io, os, re, subprocess, sys
import fitz
from PIL import Image, ImageOps
B = "/ifs1/User/zhouman/project9-v5-crc-atlas"
src = io.open(B + "/manuscript/crc_smr_atlas_rev3.tex", encoding="utf-8").read()

pre = src[:src.index("\\begin{document}")]
pre = re.sub(r"\\externaldocument\{[^}]*\}\n?", "", pre)
pre = pre.replace("\\usepackage[margin=2.5cm]{geometry}",
                  "\\usepackage[paperwidth=6.69in,paperheight=64in,margin=0.5cm]{geometry}")
pre = pre.replace("\\usepackage{graphicx}", "\\usepackage{graphicx}\n\\usepackage{float}")
pre = pre.replace("\\graphicspath{{results/figures_rev3/}{results/figures/}{results/phase6_scrna/}}",
                  "\\graphicspath{{results/figures_rev3/}{results/figures/}{results/phase6_scrna/}}")

i = src.index("\\begin{subfigure}{\\textwidth}", src.index("\\begin{figure}[htbp]"))
j = src.index("\\end{subfigure}", i) + len("\\end{subfigure}")
block = src[i:j]
# drop \caption{...} (balanced)
out, k = [], 0
while k < len(block):
    m = re.compile(r"\\caption(?:\[[^\]]*\])?\{").search(block, k)
    if not m: out.append(block[k:]); break
    out.append(block[k:m.start()]); d = 1; p = m.end()
    while p < len(block) and d:
        if block[p] == "{": d += 1
        elif block[p] == "}": d -= 1
        p += 1
    k = p
block = "".join(out)

tex = pre + "\\begin{document}\n\\thispagestyle{empty}\\pagestyle{empty}\n\\begin{figure}[H]\n" + block + "\n\\end{figure}\n\\end{document}\n"
os.makedirs(B + "/submission/figures/src", exist_ok=True)
tp = B + "/submission/figures/src/Fig1a.tex"
io.open(tp, "w", encoding="utf-8").write(tex)
r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                    "-output-directory=" + B + "/submission/figures/src", tp],
                   capture_output=True, text=True)
if r.returncode != 0:
    print("pdflatex FAILED"); print("\n".join(r.stdout.split("\n")[-20:])); sys.exit(1)
print("pdflatex ok")

pdf = B + "/submission/figures/src/Fig1a.pdf"
d = fitz.open(pdf); pg = d[0]
pix = pg.get_pixmap(dpi=100)
img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples).convert("L")
bb = ImageOps.invert(img).getbbox()
sc = 72.0/100.0; pad = 8
rect = fitz.Rect(max(0,(bb[0]-pad)*sc), max(0,(bb[1]-pad)*sc),
                 min(pg.rect.x1,(bb[2]+pad)*sc), min(pg.rect.y1,(bb[3]+pad)*sc))
nd = fitz.open(); np_ = nd.new_page(width=rect.width, height=rect.height)
np_.show_pdf_page(np_.rect, d, 0, clip=rect)
out = B + "/results/figures_rev3/Fig1a_conceptual_overview.pdf"
nd.save(out, deflate=True, garbage=3)

# report + min font
d2 = fitz.open(out); p2 = d2[0]
sizes = [round(s["size"],2) for b in p2.get_text("dict")["blocks"] for l in b.get("lines",[]) for s in l.get("spans",[]) if s["text"].strip()]
print("Fig1a_conceptual_overview.pdf  %.2f x %.2f in   min font %.2f pt" % (p2.rect.width/72, p2.rect.height/72, min(sizes) if sizes else -1))
d2.close(); nd.close(); d.close()
