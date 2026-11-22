"""Map E. coli 23S positions to contig coordinates in cmalign'd Aquicultoria LSU loci, via the census column table.

usage: python -I src/aq_ptc_coords.py STO FNA COLUMNS.tsv LOCUS[,LOCUS...] POS[,POS...]
  STO = cmalign --matchonly Pfam output (match column i = census col i, 0-based, as in columns.tsv);
  locus names are the lsu_loci.fna headers (ACC|CONTIG|START|END|STRAND|hmm...). Coordinates are found by
  locating the ungapped +-10 match-column context in the locus sequence; only '+'-strand loci are handled.
"""
import re
import sys

sto, fna, colf, loci, pos = sys.argv[1:6]
cols = {}
for l in open(colf):
    f = l.rstrip("\n").split("\t")
    if f[0] == "col":
        hdr = f
        continue
    d = dict(zip(hdr, f))
    cols[int(d["ecoli_pos"])] = (int(d["col"]), d["cons"], float(d["f_cons"]))
rows = {}
for l in open(sto):
    if l.startswith("#") or not l.strip() or l.startswith("//"):
        continue
    n, s = l.split()
    rows[n] = rows.get(n, "") + s
seqs, name = {}, None
for l in open(fna):
    if l.startswith(">"):
        name = l[1:].split()[0]
        seqs[name] = ""
    else:
        seqs[name] += l.strip().upper()
ec = rows["NR_103073.1"]
print("locus\tecoli_pos\tecoli\tcons\tf_cons\tlocus_base\tcontig\tcoord")
for nm in loci.split(","):
    r, loc = rows[nm], seqs[nm]
    acc, contig, start, end, strand = nm.split("|")[:5]
    assert strand == "+"
    for p in map(int, pos.split(",")):
        i, cons, fc = cols[p]
        b = r[i].upper()
        clean = lambda x: x.replace("-", "").replace(".", "").upper().replace("U", "T")
        L, R = clean(r[max(0, i - 10):i]), clean(r[i + 1:i + 11])
        hits = [m.start() for m in re.finditer("(?=" + L + b.replace("U", "T") + R + ")", loc)] if b not in "-." else []
        coord = int(start) + hits[0] + len(L) if len(hits) == 1 else None
        print(f"{acc}\t{p}\t{ec[i].upper()}\t{cons}\t{fc:.4f}\t{b}\t{contig}\t{coord}")
