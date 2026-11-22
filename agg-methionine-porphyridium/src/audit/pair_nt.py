"""Independent re-check of codon_pair.py (audit, 2026-11-11).

Different route from codon_pair.py: each shared CDS pair is aligned at the NUCLEOTIDE level (MAFFT --globalpair, and,
as a second aligner, a codon-unaware Biopython global nucleotide alignment), and a codon of A is paired with a codon of
B only when all three of its bases align, gap-free, to the three bases of one in-frame codon of B. A third route aligns
proteins with MAFFT after masking AGG as X in both isolates (so the aligner cannot favour AGG<->Arg or AGG<->Met).

Also checks CDS extraction: own table-4 translation == the record's /translation qualifier.

usage: python -I pair_nt.py A.gb B.gb outdir
"""
import sys, os, subprocess, tempfile
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

ARG = {'AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'}


def cds(path):
    rec = SeqIO.read(path, 'genbank')
    out, dup = {}, []
    for f in rec.features:
        if f.type != 'CDS':
            continue
        g = f.qualifiers.get('gene', [None])[0]
        if g is None:
            continue
        g = g.lower()
        nt = str(f.extract(rec.seq)).upper()
        assert f.qualifiers.get('codon_start', ['1'])[0] == '1'
        prot = str(Seq(nt).translate(table=4))
        tq = f.qualifiers.get('translation', [''])[0]
        ok = prot.rstrip('*') == tq or prot[1:].rstrip('*') == tq[1:]  # alt starts translate as M in qualifier
        if g in out:
            dup.append(g)
        out[g] = (nt, ok)
    return out, dup


def mafft_pair(s1, s2, kind):
    d = tempfile.mkdtemp()
    p = os.path.join(d, 'in.fa')
    with open(p, 'w') as fh:
        fh.write(f'>a\n{s1}\n>b\n{s2}\n')
    args = ['mafft', '--globalpair', '--maxiterate', '1000', '--quiet']
    args += ['--nuc'] if kind == 'nt' else ['--amino']
    r = subprocess.run(args + [p], capture_output=True, text=True, check=True)
    seqs, cur = {}, None
    for line in r.stdout.splitlines():
        if line.startswith('>'):
            cur = line[1:].strip(); seqs[cur] = ''
        else:
            seqs[cur] += line.strip().upper()
    return seqs['a'], seqs['b']


def bio_nt_pair(s1, s2):
    al = PairwiseAligner()
    al.mode = 'global'
    al.match_score, al.mismatch_score = 2, -1
    al.open_gap_score, al.extend_gap_score = -6, -1
    a = al.align(s1, s2)[0]
    return a[0], a[1]


def codon_pairs_from_nt(ga, gb):
    """gapped nt strings -> list of (i, j) codon index pairs where A codon i maps base-for-base to B codon j"""
    ia = ib = 0
    amap = {}  # A base pos -> B base pos (None if gap)
    for x, y in zip(ga, gb):
        if x != '-' and y != '-':
            amap[ia] = ib
        elif x != '-':
            amap[ia] = None
        if x != '-':
            ia += 1
        if y != '-':
            ib += 1
    out = []
    for i in range(ia // 3):
        b = [amap.get(3 * i + k) for k in range(3)]
        if None in b:
            continue
        if b[1] == b[0] + 1 and b[2] == b[0] + 2 and b[0] % 3 == 0:
            out.append((i, b[0] // 3))
    return out


def codon_pairs_from_prot(ga, gb):
    ia = ib = 0
    out = []
    for x, y in zip(ga, gb):
        if x != '-' and y != '-':
            out.append((ia, ib))
        if x != '-':
            ia += 1
        if y != '-':
            ib += 1
    return out


def masked_prot(nt):
    cs = [nt[k:k + 3] for k in range(0, len(nt) - len(nt) % 3, 3)]
    return ''.join('X' if c == 'AGG' else str(Seq(c).translate(table=4)).replace('*', 'X') for c in cs)


def tally(pairs_by_gene, a, b, label, fh, sites_fh=None):
    tot = defaultdict(Counter)
    same = n = 0
    for g, pairs in sorted(pairs_by_gene.items()):
        na, nb = a[g][0], b[g][0]
        for i, j in pairs:
            ca, cb = na[3 * i:3 * i + 3], nb[3 * j:3 * j + 3]
            if len(ca) < 3 or len(cb) < 3:
                continue
            n += 1
            same += ca == cb
            tot[ca][cb] += 1
            if sites_fh is not None and ca == 'AGG':
                sites_fh.write(f'{label}\t{g}\t{i + 1}\t{ca}\t{cb}\n')
    fh.write(f'\n## {label}\n')
    fh.write(f'codon-paired positions {n}, identical codons {same} ({same / n:.3f})\n')
    agg = tot['AGG']
    fh.write(f'B at A=AGG ({sum(agg.values())}): {dict(agg.most_common())}\n')
    fh.write(f'   -> AGG {agg["AGG"]}, ATG {agg["ATG"]}, other Arg {sum(v for k, v in agg.items() if k in ARG and k != "AGG")}, '
             f'else {sum(v for k, v in agg.items() if k not in ARG and k != "ATG")}\n')
    oth = Counter()
    for c in ('AGA', 'CGA', 'CGC', 'CGG', 'CGT'):
        oth.update(tot[c])
    nsite = sum(oth.values())
    nsame = sum(tot[c][c] for c in ('AGA', 'CGA', 'CGC', 'CGG', 'CGT'))
    darg = sum(v for k, v in oth.items() if k in ARG) - nsame
    fh.write(f'B at A in AGA/CGN ({nsite}): same codon {nsame}, different Arg codon {darg} (of which AGG {oth["AGG"]}), '
             f'ATG {oth["ATG"]}, else {nsite - nsame - darg - oth["ATG"]}: {dict(oth.most_common())}\n')
    atg = tot['ATG']
    fh.write(f'B at A=ATG ({sum(atg.values())}): {dict(atg.most_common())}\n')
    return tot


def main(A, B, outdir):
    a, da = cds(A)
    b, db = cds(B)
    shared = sorted(set(a) & set(b))
    fh = open(os.path.join(outdir, 'pair_nt.txt'), 'w')
    sf = open(os.path.join(outdir, 'pair_nt_AGG_sites.tsv'), 'w')
    sf.write('method\tgene\tcodon_index_A\tA\tB\n')
    fh.write(f'A={A}\nB={B}\nshared genes ({len(shared)}): {shared}\n')
    fh.write(f'duplicate gene names: A {da}, B {db}\n')
    fh.write(f'translation check (own table-4 translation == /translation): '
             f'A {sum(v[1] for v in a.values())}/{len(a)}, B {sum(v[1] for v in b.values())}/{len(b)}; '
             f'failures A {[g for g, v in a.items() if not v[1]]}, B {[g for g, v in b.items() if not v[1]]}\n')
    res = {}
    for label, fn in (('mafft_nt', lambda x, y: codon_pairs_from_nt(*mafft_pair(x, y, 'nt'))),
                      ('biopython_nt', lambda x, y: codon_pairs_from_nt(*bio_nt_pair(x, y))),
                      ('mafft_prot_AGGmasked', lambda x, y: codon_pairs_from_prot(*mafft_pair(masked_prot(x), masked_prot(y), 'aa')))):
        pairs = {g: fn(a[g][0], b[g][0]) for g in shared}
        res[label] = tally(pairs, a, b, label, fh, sf)
    fh.close(); sf.close()


if __name__ == '__main__':
    main(*sys.argv[1:])
