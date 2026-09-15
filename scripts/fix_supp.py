import re, sys
P = "manuscript/crc_smr_supplementary.tex"
s = open(P, encoding="utf-8").read()
n = 0
s, k = re.subn(r"Table~\\ref\{tab:tier1\}", "Table~1", s); n += k
s, k = re.subn(r"\(Table~\\ref\{tab:mr_sensitivity\}\)", "(Table~3)", s); n += k
if n != 2:
    print("expected 2 replacements, got", n); sys.exit(1)
open(P, "w", encoding="utf-8").write(s)
print("fixed", n, "references; remaining tab: refs:", len(re.findall(r"\\ref\{tab:", s)))