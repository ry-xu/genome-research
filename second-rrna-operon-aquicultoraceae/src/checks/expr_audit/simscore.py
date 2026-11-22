"""Score simulated reads mapped to ref.fa. usage: simscore.py CLASS.tsv SIM.fq SIM.paf
Recovery = fraction of simulated reads with >=1 alignment that passes the pipeline filter (id>=0.98, alen>=60) on a CG locus of the SAME class.
Also counts reads captured by the other class or by decoys/other loci."""
import sys, collections
CG = ["PFKS", "PFIW", "MNXI", "JAACTV", "PFTP", "PFNG", "PEXG", "PFFU"]
cls = {l.split("\t")[0]: l.split("\t")[1] for l in open(sys.argv[1])}
iscg = lambda t: "|" in t and t in cls and any(k in t.split("|")[1] for k in CG)
src = {}
for i, l in enumerate(open(sys.argv[2])):
    if i % 4 == 0:
        q = l[1:].strip(); pp = q.split("|"); loc = "|".join(pp[1:7]); L, e = pp[8], pp[9]
        src[q] = (cls[loc], int(L), float(e))
res = collections.defaultdict(lambda: collections.Counter())
best = {}
for l in open(sys.argv[3]):
    f = l.split("\t"); q, t = f[0], f[5]; m, al = int(f[9]), int(f[10])
    sc, L, e = src[q]
    passf = al >= 60 and m / al >= 0.98
    if iscg(t) and passf: tag = "same" if cls[t] == sc else "OTHER_CLASS"
    elif passf: tag = "decoy_or_nonCG_pass"
    else: tag = "fail_filter"
    rank = {"same": 0, "OTHER_CLASS": 1, "decoy_or_nonCG_pass": 2, "fail_filter": 3}[tag]
    if q not in best or rank < best[q]: best[q] = rank
names = ["same", "OTHER_CLASS", "decoy_or_nonCG_pass", "fail_filter", "unmapped"]
for q, (sc, L, e) in src.items():
    res[(sc, L, e)][names[best.get(q, 4)]] += 1
print("class\tfraglen\terr\tn\t" + "\t".join(names) + "\trecovery")
for key in sorted(res):
    c = res[key]; n = sum(c.values())
    print("\t".join(map(str, key)) + f"\t{n}\t" + "\t".join(str(c[x]) for x in names) + f"\t{c['same'] / n:.3f}")
