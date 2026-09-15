import re, sys
P = "manuscript/crc_smr_supplementary.tex"
s = open(P, encoding="utf-8").read()
old = "Expression of the 34 coloc-passing genes detected in Visium data across the\nsix ME zones;"
new = "Expression of the 15 Tier 1 genes detected in Visium data across the six ME\nzones;"
if s.count(old) != 1:
    print("caption anchor not found (%d)" % s.count(old)); sys.exit(1)
s = s.replace(old, new, 1)
open(P, "w", encoding="utf-8").write(s)
print("fixed caption gene count (34 -> 15)")