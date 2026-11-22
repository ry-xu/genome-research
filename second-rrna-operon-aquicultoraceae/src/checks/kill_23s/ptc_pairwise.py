"""Independent (structure-free) check of PTC base calls: Biopython global pairwise alignments
E.coli 23S -> each locus directly, and E.coli -> same-genome normal copy -> divergent copy (composed).
Prints +-6 nt windows around PTC positions (E. coli numbering) for each locus, and flags whether the
direct and composed mappings agree."""
import sys
from Bio import Align
def rd(p):
    return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper().replace("U", "T")
E = rd("s_ecoli23.fa")
seqs = {k: rd(f"s_{k}.fa") for k in ["dx_norm23", "dx_div23", "kmv_div23", "cad_div23"]}
al = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
al.end_gap_score = 0
def mapping(a, b):
    """dict: 0-based pos in a -> 0-based pos in b (aligned columns only)"""
    aln = al.align(a, b)[0]
    m = {}
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1): m[s1 + k] = s2 + k
    return m, aln.score
maps = {}
for k, s in seqs.items():
    maps[("E", k)], sc = mapping(E, s)
maps[("n", "dx_div23")], _ = mapping(seqs["dx_norm23"], seqs["dx_div23"])
maps[("n", "kmv_div23")], _ = mapping(seqs["dx_norm23"], seqs["kmv_div23"])
maps[("n", "cad_div23")], _ = mapping(seqs["dx_norm23"], seqs["cad_div23"])
pos = [2057, 2058, 2059, 2061, 2062, 2251, 2252, 2253, 2447, 2449, 2450, 2451, 2452, 2453, 2503, 2504, 2505, 2506, 2553, 2583, 2584, 2585, 2602]
W = 6
def win(s, m, p):
    q = m.get(p - 1)
    if q is None: return "-" * (2 * W + 1), "-"
    return s[max(0, q - W):q].lower() + s[q] + s[q + 1:q + 1 + W].lower(), s[q]
print("ecoli_pos\tE.coli\t" + "\t".join(seqs) + "\tcomposed_div(dx,kmv,cad)\tagree")
for p in pos:
    row = [str(p), win(E, {i: i for i in range(len(E))}, p)[0]]
    direct = {}
    for k, s in seqs.items():
        w, b = win(s, maps[("E", k)], p); row.append(w); direct[k] = b
    comp = []
    for k in ["dx_div23", "kmv_div23", "cad_div23"]:
        q = maps[("E", "dx_norm23")].get(p - 1)
        r = maps[("n", k)].get(q) if q is not None else None
        comp.append(seqs[k][r] if r is not None else "-")
    agree = all(comp[i] == direct[k] for i, k in enumerate(["dx_div23", "kmv_div23", "cad_div23"]))
    row += [",".join(comp), str(agree)]
    print("\t".join(row))
