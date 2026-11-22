"""T5 reference (RULES_T5.md): Aquicultoria 16S loci (+-150 nt) with classes, census 16S decoys, Deneb's 23S ref.

usage: uv run python src/t5_ref.py out/t5
Writes ref16.fa, class16.tsv (locus, class, 16S start, 16S end within the locus, 0-based half-open), ref.fa (23S + 16S).
"""
import json, os, subprocess, sys, tempfile
import t4_census as C

CG = ['GCA_009993245.1', 'GCF_001871795.1', 'GCF_002779175.1', 'GCF_002784845.1', 'GCF_002791875.1',
      'GCA_002783345.1', 'GCA_002782845.1', 'GCA_002789655.1']
REF23 = 'data/external/expr_ref.fa'
FL = 150
IDX = {l.split('\t')[0]: l.split('\t')[1] for l in open(C.ATLAS + 'index/genomes.tsv')}


def hits16(seqs, tmp):
    f = os.path.join(tmp, 'g.fa')
    with open(f, 'w') as o:
        for k, v in seqs.items(): o.write(f'>{k}\n{v}\n')
    tbl = os.path.join(tmp, 'g.tbl')
    subprocess.run([C.CMSEARCH, '--cpu', '1', '--hmmonly', '--noali', '-E', '1e-10', '--tblout', tbl, C.CM, f],
                   check=True, stdout=subprocess.DEVNULL)
    hs = []
    for l in open(tbl):
        if l[0] == '#': continue
        x = l.split(); mf, mt, a, b, st = int(x[5]), int(x[6]), int(x[7]), int(x[8]), x[9]
        if mt - mf + 1 >= 1000:
            hs.append((x[0], min(a, b), max(a, b), st, mf, mt, float(x[14])))
    keep = []
    for h in sorted(hs, key=lambda h: -h[6]):
        if not any(k[0] == h[0] and k[3] == h[3] and min(k[2], h[2]) - max(k[1], h[1]) > 0.5 * min(k[2] - k[1], h[2] - h[1]) for k in keep):
            keep.append(h)
    return keep


def locus(seqs, h):
    c, lo, hi, st = h[:4]; s = seqs[c]
    a, b = max(0, lo - 1 - FL), min(len(s), hi + FL)
    seq = s[a:b]; s0, s1 = lo - 1 - a, hi - a
    if st == '-':
        seq = C.rc(seq); s0, s1 = len(seq) - s1, len(seq) - s0
    return seq, s0, s1


def refine(seqs, h, tmp, wide=400):
    """T5-A1: 16S span of a target locus from a full-CM cmsearch of the hit +-400 nt; locus = span +-150."""
    c, lo, hi, st = h[:4]; s = seqs[c]
    a, b = max(0, lo - 1 - wide), min(len(s), hi + wide)
    w = s[a:b] if st == '+' else C.rc(s[a:b])
    f = os.path.join(tmp, 'w.fa'); open(f, 'w').write(f'>w\n{w}\n'); t = os.path.join(tmp, 'w.tbl')
    subprocess.run([C.CMSEARCH, '--cpu', '1', '--noali', '-E', '1e-5', '--toponly', '--tblout', t, C.CM, f],
                   check=True, stdout=subprocess.DEVNULL)
    x = max((l.split() for l in open(t) if l[0] != '#'), key=lambda x: float(x[14]))
    p0, p1 = int(x[7]) - 1, int(x[8])  # 0-based half-open in the oriented window
    q0, q1 = max(0, p0 - FL), min(len(w), p1 + FL)
    return w[q0:q1], p0 - q0, p1 - q0, (int(x[5]), int(x[6]))


def ident(q, ref_fa, tmp):
    """best (matches/alignment columns) of q vs ref over alignments >= 1,000 columns."""
    qf = os.path.join(tmp, 'q.fa'); open(qf, 'w').write(f'>q\n{q}\n')
    out = subprocess.run(['minimap2', '-c', '-x', 'asm20', ref_fa, qf], capture_output=True, text=True, check=True).stdout
    best = {}
    for l in out.splitlines():
        f = l.split('\t')
        if int(f[10]) >= 1000:
            best[f[5]] = max(best.get(f[5], 0), int(f[9]) / int(f[10]))
    return best


def main(od):
    os.makedirs(od, exist_ok=True); tmp = tempfile.mkdtemp(dir=od)
    aq = open('out/aq65.txt').read().split()
    loci = []
    for acc in aq:
        seqs = C.read_fasta(C.ATLAS + IDX[acc])
        for h in hits16(seqs, tmp):
            seq, s0, s1 = locus(seqs, h)
            if acc in CG:
                seq, s0, s1, m = refine(seqs, h, tmp)
                print('refined', acc, h[:6], 'full CM model', m, '16S length', s1 - s0, flush=True)
            loci.append((f'{acc}|{h[0]}|{h[1]}|{h[2]}|{h[3]}|cm{h[4]}-{h[5]}', seq, s0, s1, acc, h))
    sec = [l for l in loci if l[4] == 'GCA_009993245.1']
    print('A. secundus GCA_009993245.1 16S hits:', [(l[0], l[3] - l[2]) for l in sec], flush=True)
    div = [l for l in sec if l[5][0] == 'JAACTV010000069.1']; nrm = [l for l in sec if l[5][0] != 'JAACTV010000069.1']
    assert len(div) == 1 and len(nrm) == 1, (div, nrm)
    rf = os.path.join(od, 'sec16.fa')
    with open(rf, 'w') as o:
        for name, l in (('div', div[0]), ('norm', nrm[0])): o.write(f'>{name}\n{l[1][l[2]:l[3]]}\n')
    rows = []
    for name, seq, s0, s1, acc, h in loci:
        cls = 'decoy'
        if acc in CG:
            b = ident(seq[s0:s1], rf, tmp)
            cls = 'div' if b.get('div', 0) >= 0.97 else 'norm' if b.get('norm', 0) >= 0.97 else 'other'
            print(name, cls, {k: round(v, 4) for k, v in b.items()}, flush=True)
        rows.append((name, cls, s0, s1, seq))
    # census decoys: first copy classed 'same' per multi-copy species, no flank
    seen = set(r[4] for r in rows)
    n_dec = 0
    with open(os.path.join(od, 'ref16.fa'), 'w') as o, open(os.path.join(od, 'class16.tsv'), 'w') as t:
        for name, cls, s0, s1, seq in rows:
            o.write(f'>{name}\n{seq}\n'); t.write(f'{name}\t{cls}\t{s0}\t{s1}\n')
        for l in open('out/t4/census_v2.jsonl'):
            r = json.loads(l)
            if r['acc'] in aq or not r.get('copies'): continue
            c = next((c for c in r['copies'] if c['cls'] == 'same'), None)
            if c is None: continue
            s = C.read_fasta(C.ATLAS + IDX[r['acc']])[c['contig']][c['start'] - 1:c['end']]
            if c['strand'] == '-': s = C.rc(s)
            nm = f"dec|{r['acc']}|{c['contig']}|{c['start']}|{c['end']}"
            o.write(f'>{nm}\n{s}\n'); t.write(f'{nm}\tdecoy\t0\t{len(s)}\n'); n_dec += 1
    with open(os.path.join(od, 'ref.fa'), 'w') as o:
        for p in (REF23, os.path.join(od, 'ref16.fa')):
            for l in open(p): o.write(l)
    print('loci', len(rows), 'census decoys', n_dec, 'classes', {k: sum(r[1] == k for r in rows) for k in ('div', 'norm', 'other', 'decoy')})


if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(__file__))
    main(sys.argv[1])
