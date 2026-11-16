# Probe how cpos/strand map onto the FASTA (exploration only)
import csv, sys
G = "data/rna_chicken/genomes/"
SITES = "results/rna_chicken/sites.tsv"
def readfa(p):
    d = {}; name = None; buf = []
    for line in open(p):
        if line.startswith('>'):
            if name: d[name] = ''.join(buf).upper()
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: d[name] = ''.join(buf).upper()
    return d
comp = str.maketrans('ACGTN', 'TGCAN')
rc = lambda s: s.translate(comp)[::-1]
rows = [r for r in csv.DictReader(open(SITES), delimiter='\t') if r['genome'] == 'GCA_947856825.1']
fa = readfa(G + 'GCA_947856825.1.fna')
ok = {}
for r in rows[:2000]:
    c = int(r['cpos']); s = fa[r['contig']]
    for name, f in [('0b+', lambda: s[c:c+3]), ('1b+', lambda: s[c-1:c+2]),
                    ('0b-', lambda: rc(s[c-2:c+1])), ('1b-', lambda: rc(s[c-3:c]))]:
        if r['strand'] == name[-1]:
            ok.setdefault(name, [0, 0])
            ok[name][0] += f() == r['codon']; ok[name][1] += 1
print(ok)
