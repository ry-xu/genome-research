"""Model-position distribution of CG hits (id>=.98, alen>=60), duplicates, region-matched counts.
usage: positions.py CLASS.tsv PAF..."""
import sys, collections
CG = ["PFKS", "PFIW", "MNXI", "JAACTV", "PFTP", "PFNG", "PEXG", "PFFU"]
cls = {l.split("\t")[0]: l.split("\t")[1] for l in open(sys.argv[1])}
for paf in sys.argv[2:]:
    hits = []
    for l in open(paf):
        f = l.split("\t"); t = f[5]
        if t not in cls or not any(k in t.split("|")[1] for k in CG): continue
        m, al = int(f[9]), int(f[10])
        if al < 60 or m / al < 0.98: continue
        hf, ht = map(int, t.split("|")[5][3:].split("-")); tl = int(f[6])
        ts, te = int(f[7]), int(f[8])
        a = hf + ts / tl * (ht - hf); b = hf + te / tl * (ht - hf)
        hits.append((cls[t], f[0], t, ts, te, a, b, int(f[2]) < 151, f[4]))
    print("==", paf)
    for c in ("div", "norm"):
        h = [x for x in hits if x[0] == c]
        if not h: print(" ", c, 0); continue
        bins = collections.Counter(int((x[5] + x[6]) / 2 // 250) * 250 for x in h)
        inwin = sum(1 for x in h if x[5] >= 558 and x[6] <= 2122)
        inwin8 = sum(1 for x in h if x[5] >= 8 and x[6] <= 2122)
        spots_win = len({x[1] for x in h if x[5] >= 558 and x[6] <= 2122})
        uniq = len({(x[2].split("|")[0], x[3], x[8]) for x in h})
        # duplicates across loci are not comparable (different genomes), so also key by model start rounded
        uniq_model = len({(round(x[5]), x[8], x[7]) for x in h})
        r1 = collections.Counter((x[7], x[8]) for x in h)
        print(f"  {c}: aln={len(h)} spots={len({x[1] for x in h})} fully_in_model[558,2122]={inwin} (spots {spots_win}) fully_in[8,2122]={inwin8} distinct(model_start,strand,mate)={uniq_model}")
        print(f"     (isR1, strand) counts: {dict(r1)}")
        print(f"     model-pos 250-bins: {dict(sorted(bins.items()))}")
