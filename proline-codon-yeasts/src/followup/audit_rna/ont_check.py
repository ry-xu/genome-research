"""(e) RNA RDDs vs same-strain (PYCC 7102) ONT WGS. pileup parse copied from src/rdd.py with BQ as an argument.
usage: ont_check.py ref.fa rdd.tsv ont.bam rna.bam genes.tsv"""
import sys, csv, subprocess
from collections import Counter
fa, rddtsv, ont, rna, genes = sys.argv[1:6]
def pileup(bam, Q):
    cmd = ['samtools', 'mpileup', '-B', '-Q', str(Q), '-q', '20', '-d', '1000000', '-f', fa, bam]
    out = {}
    for line in subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines():
        c, p, ref, d, bases = line.split('\t')[:5]
        cnt, i = Counter(), 0
        while i < len(bases):
            b = bases[i]
            if b == '^': i += 2; continue
            if b in '+-':
                j = i + 1
                while j < len(bases) and bases[j].isdigit(): j += 1
                i = j + int(bases[i + 1:j]); continue
            if b in '.,': cnt[ref.upper()] += 1
            elif b.upper() in 'ACGT': cnt[b.upper()] += 1
            i += 1
        out[int(p)] = cnt
    return out
g = ''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper()
po = pileup(ont, 7); pr = pileup(rna, 20)
COMP = str.maketrans('ACGT', 'TGCA')
n_conf = n_tot = 0
print('pos\tgene\tcp\ttype\tRNA alt/depth\tONT counts (genome + strand)\tONT frac of RNA alt\tverdict')
for r in csv.DictReader(open(rddtsv), delimiter='\t'):
    p = int(r['pos']); ref = g[p - 1]
    alt = r['type'].split('>')[1].replace('U', 'T')
    if r['strand'] == '-1': alt = alt.translate(COMP)
    c = po.get(p, Counter()); d = sum(c.values())
    f = c[alt] / d if d else float('nan')
    v = 'DNA variant (strain)' if d >= 10 and f >= 0.5 else ('DNA mixed' if d >= 10 and f >= 0.05 else ('no DNA support' if d >= 10 else 'ONT depth<10'))
    n_tot += 1; n_conf += v.startswith('DNA variant')
    print(f"{p}\t{r['gene']}\t{r['codon_pos']}\t{r['type']}\t{r['alt_n']}/{r['depth']}\t{dict(c)}\t{f:.2f}\t{v}")
print(f'RDDs confirmed as PYCC 7102 DNA variants: {n_conf} of {n_tot}')
# converse: coding positions where ONT (depth>=20) has a non-ref allele >= 0.8 and RNA depth >= 20
cod = {int(x['pos1']) for x in csv.DictReader(open(genes), delimiter='\t')}
dv = [p for p in cod if sum(po.get(p, Counter()).values()) >= 20 and max([v for b, v in po[p].items() if b != g[p-1]] or [0]) / sum(po[p].values()) >= 0.8]
both = [p for p in dv if sum(pr.get(p, Counter()).values()) >= 20]
print(f'coding positions where PYCC 7102 ONT differs from BK074543 (alt>=0.8, depth>=20): {len(dv)}; of these with RNA depth>=20: {len(both)}')
# ONT depth across coding positions
d20 = sum(1 for p in cod if sum(po.get(p, Counter()).values()) >= 20)
print(f'coding positions with ONT depth>=20 (BQ>=7): {d20} of {len(cod)}')
