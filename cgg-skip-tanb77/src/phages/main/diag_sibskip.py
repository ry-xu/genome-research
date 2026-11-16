"""Post-freeze diagnostic: the skip signature (in last 16 codons, +1 frame from 2nd base open >= 50 codons) for CGG
and other codons, per set, phages in GC bins 24-30. usage: python -I diag_sibskip.py OUT.tsv"""
import csv, sys
from collections import defaultdict
import pyrodigal
sys.path.insert(0, "src")
from codon_stats import read_fasta, revcomp, open_len, STOPS
CODS = ["CGG", "CGA", "CGT", "CGC", "AGG", "AGA", "GGG", "CCG", "GGA", "AAG"]
sets = {r["id"]: r["set"] for r in csv.DictReader(open("data/groups/run_ids.tsv"), delimiter="\t")
        if r["set"] in ("T1", "T0", "Tf", "R") and 24 <= int(r["gc_bin"]) <= 30}
gf = pyrodigal.GeneFinder(meta=True); stops = STOPS[11]
tot = defaultdict(lambda: defaultdict(int)); sk = defaultdict(lambda: defaultdict(int))
for name, seq in read_fasta(open("data/phages.fna"), set(sets)):
    if len(seq) < 5000:
        continue
    s0 = sets[name]; rc = revcomp(seq)
    for g in gf.find_genes(seq.encode()):
        if g.partial_begin or g.partial_end:
            continue
        s = seq[g.begin - 1:g.end] if g.strand == 1 else revcomp(seq[g.begin - 1:g.end])
        if len(s) < 300 or len(s) % 3:
            continue
        cod = [s[i:i + 3] for i in range(0, len(s) - 3, 3)]
        for i, c in enumerate(cod):
            if c not in CODS:
                continue
            tot[s0][c] += 1
            if len(cod) - i <= 16:
                st = (g.begin - 1 + 3 * i + 1, seq) if g.strand == 1 else ((len(seq) - g.end) + 3 * i + 1, rc)
                sk[s0][c] += open_len(st[1], st[0], stops) >= 50
w = open(sys.argv[1], "w"); w.write("set\tcodon\tn\tskiplike\tshare\n")
for s0 in ("T1", "T0", "Tf", "R"):
    print(s0, "  ".join(f"{c}:{sk[s0][c]}/{tot[s0][c]}={sk[s0][c] / max(1, tot[s0][c]):.4f}" for c in CODS))
    for c in CODS:
        w.write(f"{s0}\t{c}\t{tot[s0][c]}\t{sk[s0][c]}\t{sk[s0][c] / max(1, tot[s0][c]):.5f}\n")
