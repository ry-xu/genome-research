# Probe: reproduce old_frame_codons / plus1_frame_codons of site_class.tsv (exploration only)
import csv, sys, collections
G = "data/rna_chicken/genomes/"
R = "results/rna_chicken/"
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
STOP = {'TAA', 'TAG', 'TGA'}
def run(seq, start, cap=60):
    n = 0; i = start
    while n < cap:
        cod = seq[i:i+3]
        if len(cod) < 3: return n, 'end'
        if cod in STOP: return n, 'stop'
        n += 1; i += 3
    return n, 'cap'
cls = {r['site_id']: r for r in csv.DictReader(open(R + 'site_class.tsv'), delimiter='\t')}
gen = sys.argv[1]
fa = readfa(G + gen + '.fna')
rcfa = {}
tab = collections.Counter(); bad = []
for r in csv.DictReader(open(R + 'sites.tsv'), delimiter='\t'):
    if r['genome'] != gen or r['site_id'] not in cls: continue
    c = int(r['cpos']); s = fa[r['contig']]; L = len(s)
    if r['strand'] == '+':
        seq = s; p = c
    else:
        if r['contig'] not in rcfa: rcfa[r['contig']] = rc(s)
        seq = rcfa[r['contig']]; p = L - 1 - c
    assert seq[p:p+3] == r['codon']
    for off_old in (3,):
        for off_p1 in (1, 4):
            o, ow = run(seq, p + off_old); q, qw = run(seq, p + off_p1)
            k = cls[r['site_id']]
            match = (o == int(k['old_frame_codons']) and q == int(k['plus1_frame_codons']))
            tab[(off_old, off_p1, match)] += 1
            if not match and off_p1 == 1: bad.append((r['site_id'], o, ow, q, qw, k['old_frame_codons'], k['plus1_frame_codons'], k['class']))
print(tab)
print(bad[:20])
