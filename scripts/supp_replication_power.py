import csv, math, io, os

P = '/ifs1/User/zhouman/project9-v5-crc-atlas'
src = os.path.join(P, 'results', 'phase4_finngen_replication.csv')

def Phi(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))

rows = list(csv.DictReader(io.open(src, encoding='utf-8')))
out = []
for r in rows:
    gene = r['gene']
    try:
        b = float(r['b_SMR_discovery'])
        se = float(r['ivw_se'])
    except (ValueError, KeyError):
        continue
    if se <= 0:
        continue
    z = abs(b) / se
    pw = Phi(z - 1.96) + Phi(-z - 1.96)
    out.append((gene, abs(b), se, z, pw, r.get('replication', '')))

print('gene\tabs_b_disc\tSE_finngen\tz\tpower\treplication')
for g, b, se, z, pw, rep in out:
    print('%s\t%.4f\t%.4f\t%.3f\t%.3f\t%s' % (g, b, se, z, pw, rep))

n = len(out)
if n:
    pws = [x[4] for x in out]
    print()
    print('n genes          : %d' % n)
    print('mean power       : %.3f' % (sum(pws) / n))
    print('min / max power  : %.3f / %.3f' % (min(pws), max(pws)))
    print('median power     : %.3f' % sorted(pws)[n // 2])
    print('observed repl    : %d/%d full' % (sum(1 for x in out if x[5].strip().lower() == 'full'), n))