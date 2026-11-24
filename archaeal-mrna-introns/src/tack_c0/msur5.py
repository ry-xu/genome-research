"""Look for the c0 site in extra Methanosuratincola genomes: upstream flank CACGGTAT (exact 8-mer preceded by a
fuzzy context match) followed within 60 nt by TTAGC; prints what lies between. usage: msur5.py GENOME.fna.gz ..."""
import gzip, sys, re
rc = lambda t: t.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
UP = 'CCCAGCTTCACGGTAT'
for f in sys.argv[1:]:
    G = {}; n = None
    for l in gzip.open(f, 'rt'):
        if l[0] == '>': n = l[1:].split()[0]; G[n] = []
        else: G[n].append(l.strip().upper())
    hits = []
    for c, s in G.items():
        s = ''.join(s)
        for strand, t in (('+', s), ('-', rc(s))):
            for m in re.finditer('CACGGTAT', t):
                ctx = t[m.start() - 8:m.end()]
                mm = sum(a != b for a, b in zip(ctx, UP)) if len(ctx) == 16 else 99
                if mm > 3: continue
                tail = t[m.end():m.end() + 70]
                k = tail.find('TTAGC')
                hits.append((c, strand, m.end(), mm, tail[:k] if k >= 0 else None, tail[:40]))
    print(f.split('/')[-1], sum(len(''.join(v)) for v in G.values()), 'bp;', len(hits), 'site(s)')
    for h in hits:
        ins = h[4]
        print(f'  {h[0]} {h[1]} end {h[2]} flank mismatches {h[3]}: insert {ins!r} ({len(ins) if ins is not None else "-"} nt) next40 {h[5]}')
