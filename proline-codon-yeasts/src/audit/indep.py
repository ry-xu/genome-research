"""Independent check: six-frame table-4 ORFs per genome, local-aligned (Biopython, BLOSUM62) to Entelexis
sorbosivorans (PQ436977.1) reference proteins; each genome's codons projected onto reference coordinates.
At reference positions where BOTH relatives carry a CCN (Pro) codon, tally the query's codon; same for TTR (Leu)."""
import sys, glob
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.Align import PairwiseAligner, substitution_matrices

recdir, genes = sys.argv[1], sys.argv[2].split(',')
queries = sys.argv[3].split(','); rels = sys.argv[4].split(',')
STOPS = {'TAA', 'TAG'}
T4 = {c: str(Seq(c).translate(table=4)) for c in (a + b + d for a in 'TCAG' for b in 'TCAG' for d in 'TCAG')}

ref = SeqIO.read(f'{recdir}/PQ436977.1.gb', 'genbank')
refp = {}
for f in ref.features:
    if f.type == 'CDS' and f.qualifiers.get('gene', [''])[0] in genes:
        refp[f.qualifiers['gene'][0]] = f.qualifiers['translation'][0]

al = PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                     open_gap_score=-11, extend_gap_score=-1)

def orfs(seq):
    L = len(seq)
    out = []
    for strand, s in ((1, seq), (-1, reverse_complement(seq))):
        sx = s + s[:6000]
        for fr in range(3):
            start = fr
            for i in range(fr, len(sx) - 2, 3):
                if sx[i:i + 3] in STOPS:
                    if start < L and i - start >= 90:
                        out.append(sx[start:i])
                    start = i + 3
    return out

def project(acc):
    seq = str(SeqIO.read(f'{recdir}/{acc}.gb', 'genbank').seq).upper()
    os_ = orfs(seq)
    prots = [''.join(T4.get(o[3*i:3*i+3], 'X') for i in range(len(o)//3)) for o in os_]
    res = {}
    for g, rp in refp.items():
        best = {}  # ref pos -> (score, codon, aa)
        for o, p in zip(os_, prots):
            a = al.align(rp, p)[0]
            if a.score < 60:
                continue
            for (r0, r1), (q0, q1) in zip(*a.aligned):
                for k in range(r1 - r0):
                    r, q = r0 + k, q0 + k
                    if a.score > best.get(r, (-1,))[0]:
                        best[r] = (a.score, o[3*q:3*q+3], p[q])
        res[g] = {r: v[1] for r, v in best.items()}
    return res

proj = {a: project(a) for a in queries + rels}
for g in refp:
    print(f'## {g} (ref len {len(refp[g])}); covered positions:', {a: len(proj[a][g]) for a in proj})
grand = defaultdict(Counter); grandL = defaultdict(Counter)
for q in queries:
    for g in refp:
        P_ = Counter(); L_ = Counter()
        for r in proj[q][g]:
            rc = [proj[x][g].get(r) for x in rels]
            if all(c and c.startswith('CC') for c in rc):
                P_[proj[q][g][r]] += 1
            if all(c in ('TTA', 'TTG') for c in rc):
                L_[proj[q][g][r]] += 1
        grand[q].update(P_); grandL[q].update(L_)
        print(q, g, 'opposite relatives\' CCN:', dict(P_.most_common()), '| opposite relatives\' TTR:', dict(L_.most_common(6)))
print()
for q in queries:
    P_, L_ = grand[q], grandL[q]
    ctn = sum(v for k, v in P_.items() if k.startswith('CT')); ccn = sum(v for k, v in P_.items() if k.startswith('CC'))
    ctnL = sum(v for k, v in L_.items() if k.startswith('CT')); ttrL = sum(v for k, v in L_.items() if k in ('TTA', 'TTG'))
    print(f'{q} TOTAL at relatives-CCN positions: n={sum(P_.values())} CTN={ctn} CCN={ccn} other={sum(P_.values())-ctn-ccn} {dict(P_.most_common(8))}')
    print(f'{q} TOTAL at relatives-TTR positions: n={sum(L_.values())} TTR={ttrL} CTN={ctnL} {dict(L_.most_common(6))}')
