"""ARAGORN on the short tRNAscan-SE Undet calls (< 120 bp) of the pig TANB77 representatives, with 150 nt flanks.
usage: python -I undet_short.py <fna_dir> <out...>"""
import sys
from pathlib import Path

import pyaragorn

f = pyaragorn.RNAFinder(11, tmrna=False, linear=True)
for out in sys.argv[2:]:
    g = Path(out).name.split(".")[0]
    seqs, name = {}, None
    for line in open(Path(sys.argv[1]) / f"{g}.fna"):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    for line in open(out):
        x = line.split()
        if len(x) < 9 or not x[1].isdigit() or x[4] != "Undet":
            continue
        a, b = sorted((int(x[2]), int(x[3])))
        if b - a + 1 >= 120:
            continue
        s = seqs[x[0]][max(0, a - 151):b + 150]
        hits = [(h.amino_acids, h.anticodon, h.anticodon_length, h.energy) for h in f.find_rna(s)]
        print(x[0], a, b, b - a + 1, hits or "no ARAGORN call", sep="\t")
