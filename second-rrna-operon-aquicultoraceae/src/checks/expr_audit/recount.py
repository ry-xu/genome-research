"""Independent recount of PAF hits on CG loci. usage: recount.py CLASS.tsv PAF..."""
import sys, collections
CG = ["PFKS", "PFIW", "MNXI", "JAACTV", "PFTP", "PFNG", "PEXG", "PFFU"]
cls = {}
for l in open(sys.argv[1]):
    f = l.rstrip("\n").split("\t"); cls[f[0]] = f[1]
def iscg(t): return any(k in t.split("|")[1] for k in CG)
for paf in sys.argv[2:]:
    lines = [l.rstrip("\n").split("\t") for l in open(paf)]
    out = collections.Counter(); spots = collections.defaultdict(set); idh = collections.defaultdict(list)
    lens = collections.Counter(); strands = collections.Counter(); notcls = 0
    for f in lines:
        q, ql, qs, qe, st, t = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5]
        m, al = int(f[9]), int(f[10])
        if t not in cls: notcls += 1; continue
        if not iscg(t): continue
        c = cls[t]; ident = m / al
        idh[c].append(ident)
        for thr, minlen in [(0.98, 0), (0.98, 60), (0.95, 0), (0.90, 0), (0.0, 0)]:
            if ident >= thr and al >= minlen:
                out[(c, thr, minlen, "aln")] += 1
                spots[(c, thr, minlen)].add(q)
        if ident >= 0.98:
            lens[(c, al // 20 * 20)] += 1; strands[(c, st)] += 1
    print("==", paf, "lines", len(lines), "targets not in class table", notcls)
    for thr, minlen in [(0.98, 0), (0.98, 60), (0.95, 0), (0.90, 0), (0.0, 0)]:
        print(f" id>={thr} alen>={minlen}: " + "  ".join(f"{c}: aln={out[(c,thr,minlen,'aln')]} spots={len(spots[(c,thr,minlen)])}" for c in ("div", "norm")))
    both = spots[("div", 0.98, 0)] & spots[("norm", 0.98, 0)]
    print(" spots hitting both classes (id>=.98):", len(both))
    print(" strand (id>=.98):", dict(strands))
    print(" aligned length bins (id>=.98):", dict(sorted(lens.items())))
    for c in ("div", "norm"):
        v = sorted(idh[c])
        if v: print(f" {c} identity: n={len(v)} min={v[0]:.3f} median={v[len(v)//2]:.3f} n<0.98={sum(x<0.98 for x in v)}")
