"""Anticodon-loop length of every tRNAscan-SE call: the run of '.' holding the anticodon's first base, bounded by '>' and '<'.
usage: python -I loops.py <ss files...>; prints genome, name, type, anticodon, loop length, score; then a summary."""
import re
import sys
from collections import Counter

rows = []
for f in sys.argv[1:]:
    txt = open(f).read().split("\n\n")
    for b in txt:
        m = re.search(r"^(\S+) \((\d+)-(\d+)\).*\nType: (\S+)\s+Anticodon: (\S+) at (\d+)-(\d+).*Score: ([\d.]+)", b, re.M)
        if not m:
            continue
        seq = re.search(r"^Seq: (\S+)", b, re.M)[1]
        st = re.search(r"^Str: (\S+)", b, re.M)[1]
        a = int(m[6]) - 1
        if "intron" in b.lower():
            pass
        i = a
        while i > 0 and st[i - 1] == ".":
            i -= 1
        j = a
        while j < len(st) and st[j] == ".":
            j += 1
        rows.append((m[1].split("_1.")[0] if "_" in m[1] else m[1], m[1], m[4], m[5], j - i, float(m[8]),
                     seq[i:j], st[i - 1:j + 1]))
for r in rows:
    if r[4] != 7 or r[2] in ("Gly", "Undet", "Sup") or r[3] in ("CCG", "NNN"):
        print(*r, sep="\t")
print("loop lengths:", Counter(r[4] for r in rows))
print("Gly anticodons:", Counter(r[3] for r in rows if r[2] == "Gly"))
print("genomes:", len({r[0] for r in rows}), "tRNAs:", len(rows))
