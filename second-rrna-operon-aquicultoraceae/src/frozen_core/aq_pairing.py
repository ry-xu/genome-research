"""Helix-pairing test for divergent vs normal LSU copies: is the divergent rRNA still folded (compensatory changes)
or drifting like a pseudogene?

usage: python -I src/aq_pairing.py STO COLUMNS.tsv LOCI_TABLE.tsv
  LOCI_TABLE: name<TAB>class (div|norm); names are lsu_loci.sto row names.
  Uses census base-pair partners (columns.tsv 'partner', 0-based cols). For each locus and each pair of columns with
  both residues present: canonical (WC or GU) or not; and whether the pair differs from the census consensus pair at
  one side only (single) or both sides (double). A double change that stays canonical is a compensatory change.
  Null: the same counts with each column paired to a random column at the same distance class (shuffled partners).
"""
import random
import sys

CAN = {"AU", "UA", "GC", "CG", "GU", "UG"}
sto, colf, tab = sys.argv[1:4]
cons, part = {}, {}
for l in open(colf):
    f = l.rstrip("\n").split("\t")
    if f[0] == "col":
        hdr = f
        continue
    d = dict(zip(hdr, f))
    cons[int(d["col"])] = d["cons"]
    if d["partner"] not in ("", "NA", "nan", "-1"):
        part[int(d["col"])] = int(float(d["partner"]))
pairs = sorted((i, j) for i, j in part.items() if i < j)
rows = {}
for l in open(sto):
    if l.startswith("#") or not l.strip() or l.startswith("//"):
        continue
    n, s = l.split()
    rows[n] = rows.get(n, "") + s.upper().replace("T", "U")
cls = dict(l.rstrip("\n").split("\t") for l in open(tab))
rng = random.Random(1)
cols = sorted(cons)
shuf = [(i, rng.choice(cols)) for i, j in pairs]


def score(r, prs):
    n = can = single = single_can = double = double_can = 0
    for i, j in prs:
        a, b = r[i], r[j]
        if a not in "ACGU" or b not in "ACGU":
            continue
        n += 1
        c = a + b in CAN
        can += c
        da, db = a != cons[i], b != cons[j]
        if da and db:
            double += 1
            double_can += c
        elif da or db:
            single += 1
            single_can += c
    return n, can, single, single_can, double, double_can


print("name\tclass\tn_pairs\tf_canonical\tsingle\tsingle_can\tdouble\tdouble_can\tnull_f_canonical")
for nm, c in cls.items():
    n, can, s, sc, d, dc = score(rows[nm], pairs)
    nn, ncan, *_ = score(rows[nm], shuf)
    print(f"{nm}\t{c}\t{n}\t{can / n:.3f}\t{s}\t{sc}\t{d}\t{dc}\t{ncan / nn:.3f}")
