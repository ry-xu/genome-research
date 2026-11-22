"""Count RNA reads on divergent vs normal 23S loci.

usage: python -I src/expr_count.py LOCI.paf LOCI_CLASS.tsv [MIN_ID]
  LOCI.paf        minimap2 -c PAF lines whose target is a locus (acc|contig|from|to|strand|hmmF-hmmT)
  LOCI_CLASS.tsv  locus name, class (div/norm), ...
Read identity = matches / alignment length (PAF cols 10/11); reads below MIN_ID (default 0.98) are dropped.
Prints per locus and per class: reads, reads per kb of locus, and the E. coli-model span they fall in
(locus offsets are converted to model positions linearly between hmmF and hmmT, so spans are approximate).
"""
import collections
import sys


def main(paf, cls_path, min_id=0.98):
    cls = {l.split("\t")[0]: l.split("\t")[1] for l in open(cls_path)}
    n = collections.Counter()
    bins = collections.defaultdict(collections.Counter)
    tlen = {}
    for line in open(paf):
        f = line.split("\t")
        t, tl, ts, te, m, alen = f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
        if t not in cls or alen == 0 or m / alen < min_id or alen < 60:  # decoys can also contain "|"
            continue
        n[t] += 1
        tlen[t] = tl
        parts = t.split("|")
        if len(parts) > 5:  # loci found by blastn instead of nhmmer carry no model span; they are counted but not binned
            hf, ht = map(int, parts[5][3:].split("-"))
            mid = hf + (ts + te) / 2 / tl * (ht - hf)
            bins[cls.get(t, "?")][int(mid // 500) * 500] += 1
    tot = collections.Counter()
    for t in sorted(n, key=lambda x: (cls.get(x, "?"), -n[x])):
        c = cls.get(t, "?")
        tot[c] += n[t]
        print(f"{c}\t{t}\t{n[t]}\t{n[t] / tlen[t] * 1000:.1f}/kb")
    print("class totals:", dict(tot))
    for c, b in bins.items():
        print(c, "reads by model position (500-nt bins):", dict(sorted(b.items())))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 0.98)
