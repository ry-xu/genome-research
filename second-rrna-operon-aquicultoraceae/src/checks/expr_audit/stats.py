"""Binomial/upper-bound statistics on the claimed counts, by alignment and by spot, with window matching and dedup.
usage: stats.py CLASS.tsv DNA.paf RNA.paf..."""
import sys, collections, random, math
CG = ["PFKS", "PFIW", "MNXI", "JAACTV", "PFTP", "PFNG", "PEXG", "PFFU"]
cls = {l.split("\t")[0]: l.split("\t")[1] for l in open(sys.argv[1])}
def load(paf, win=None):
    aln = collections.Counter(); spots = collections.defaultdict(set); frags = collections.defaultdict(set)
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    for l in open(paf):
        f = l.split("\t"); t = f[5]
        if t not in cls or not any(k in t.split("|")[1] for k in CG): continue
        m, al = int(f[9]), int(f[10])
        if al < 60 or m / al < 0.98: continue
        hf, ht = map(int, t.split("|")[5][3:].split("-")); tl = int(f[6])
        a = hf + int(f[7]) / tl * (ht - hf); b = hf + int(f[8]) / tl * (ht - hf)
        if win and not (a >= win[0] and b <= win[1]): continue
        c = cls[t]; aln[c] += 1; spots[c].add(f[0]); per[c][f[0]].append((round(a), round(b)))
    for c in per:
        for q, iv in per[c].items():
            frags[c].add((min(x[0] for x in iv), max(x[1] for x in iv)))  # fragment = span of its mates in model coords
    return {c: (aln[c], len(spots[c]), len(frags[c])) for c in ("div", "norm")}
def ub_r(n_rna_norm, d, nm, B=200000, seed=1):
    """95% upper bound on per-copy RNA ratio r=div/norm given 0 div of n RNA, DNA div d / norm nm.
    analytic: Clopper-Pearson q_u with p point estimate; MC: q~Beta(1,n+1)(Jeffreys-ish), p~Beta(d+1,nm+1)."""
    p = d / (d + nm); n = n_rna_norm
    qu = 1 - 0.05 ** (1 / n)
    r_an = qu / (1 - qu) * (1 - p) / p
    rnd = random.Random(seed); rs = []
    for _ in range(B):
        q = rnd.betavariate(0.5, n + 0.5); pp = rnd.betavariate(d + 0.5, nm + 0.5)
        rs.append(q / (1 - q) * (1 - pp) / pp)
    rs.sort()
    return p, n * p / (1 - p), n * p, (1 - p) ** n, r_an, rs[int(0.95 * B)]
dna = sys.argv[2]; rnas = sys.argv[3:]
for label, win in [("all positions", None), ("model window [8,2122] (both classes covered)", (8, 2122)), ("model window [558,2122] (all 16 loci covered)", (558, 2122))]:
    D = load(dna, win); R = [load(r, win) for r in rnas]
    print(f"### {label}")
    print("  DNA div/norm (aln, spots, dedup-fragments):", D["div"], D["norm"])
    for r, x in zip(rnas, R): print("  RNA", r.split("/")[-1], "div", x["div"], "norm", x["norm"])
    for k, unit in enumerate(("alignments", "spots", "dedup fragments")):
        d, nm = D["div"][k], D["norm"][k]
        nr = sum(x["norm"][k] for x in R); dr = sum(x["div"][k] for x in R)
        p, e_given_norm, e_binom_total, p0, r_an, r_mc = ub_r(nr, d, nm)
        print(f"  [{unit}] DNA share p={p:.4f}; RNA div={dr} norm={nr}; expected div given norm count = {e_given_norm:.1f}; "
              f"n*p = {e_binom_total:.1f}; P(0|n,p) = {p0:.3g}; 95% UB on per-copy ratio r: analytic {r_an:.4f}, MC with DNA uncertainty {r_mc:.4f}")
        per = [(x["norm"][k], (1 - p) ** x["norm"][k]) for x in R]
        print(f"     per-sample (norm n, P(0)): {[(a, float(f'{b:.3g}')) for a, b in per]}")
        # overdispersion sensitivity: design effect 1+(m-1)rho within a sample, m = sample size
        for rho in (0.01, 0.05, 0.2):
            neff = sum(x["norm"][k] / (1 + (x["norm"][k] - 1) * rho) for x in R)
            qu = 1 - 0.05 ** (1 / neff)
            print(f"     rho={rho}: n_eff={neff:.1f}, 95% UB r={qu / (1 - qu) * (1 - p) / p:.4f}")
