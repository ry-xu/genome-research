"""Pass 3: follow-ups to analyse.py (context-matched background throughout, family decoy, taxonomy)."""
import sys, collections, random, itertools, math
sys.argv = ["x"]
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    exec(open("results/skip_context_chicken/recheck/analyse.py").read().split("# ---------------- break (a)")[0])
P.clear()
fam_of = {}
for r in csv.DictReader(open(SEL), delimiter="\t"):
    fam_of[r["acc"]] = r["family"]

pr("== Follow-ups (context-matched background: codon before in-frame CGA/CGC/CGT, per genome) ==")
o, e = obs_exp(SK, "cgn"); o2, e2 = obs_exp(IF, "cgn"); o3, e3 = obs_exp([r for r in A if r["given_class"] == "unclear"], "cgn")
pr(f"CTT: skip {o['CTT']:.0f}/{e['CTT']:.1f} = {ratio(o, e, 'CTT'):.3f}; in-frame {o2['CTT']:.0f}/{e2['CTT']:.1f} = {ratio(o2, e2, 'CTT'):.3f}; "
   f"unclear {o3['CTT']:.0f}/{e3['CTT']:.1f} = {ratio(o3, e3, 'CTT'):.3f}")
pr("five enriched (skip / in-frame):", {c: (round(ratio(o, e, c), 2), round(ratio(o2, e2, c), 2)) for c in ENRICHED})
pr("four depleted (skip / in-frame):", {c: (round(ratio(o, e, c), 2), round(ratio(o2, e2, c), 2)) for c in DEPLETED})
o4, e4 = obs_exp([r for r in A if r["given_class"] == "unclear"])
pr(f"(GGN bg) unclear class CTT {o4['CTT']:.0f}/{e4['CTT']:.1f} = {ratio(o4, e4, 'CTT'):.3f}")
pg = collections.defaultdict(lambda: [0, 0.0])
for r in SK:
    pg[r["genome"]][0] += r["up"] == "CTT"
    pg[r["genome"]][1] += BG["cgn"][r["genome"]]["CTT"]
lo, hi, _ = boot([pg[g] for g in sorted(pg)])
pr(f"bootstrap over genomes (1,000, seed 1), cgn background: {lo:.3f} - {hi:.3f}")
ex = sorted(((v[0] - v[1], g) for g, v in pg.items()), reverse=True)
for k in (5, 10, 20):
    top = {g for _, g in ex[:k]}
    oo, ee = obs_exp([r for r in SK if r["genome"] not in top], "cgn")
    pr(f"drop top {k} genomes by CTT excess (cgn bg): {ratio(oo, ee, 'CTT'):.3f}")
# all CAG-269 out
oo, ee = obs_exp([r for r in SK if genus.get(r["genome"]) != "CAG-269"], "cgn")
oo2, ee2 = obs_exp([r for r in SK if genus.get(r["genome"]) != "CAG-269"], "ggn")
pr(f"genus CAG-269 removed: CTT {oo['CTT']:.0f}/{ee['CTT']:.1f} = {ratio(oo, ee, 'CTT'):.3f} (cgn), {ratio(oo2, ee2, 'CTT'):.3f} (ggn)")
# GTDB family
for fm in sorted({fam_of.get(r["genome"]) for r in SK}):
    for nm, S in (("skip", SK), ("inframe", IF)):
        S2 = [r for r in S if fam_of.get(r["genome"]) == fm]
        if not S2:
            continue
        oo, ee = obs_exp(S2, "cgn")
        oo2, ee2 = obs_exp(S2, "ggn")
        lam = ee["CTT"]
        k = int(oo["CTT"])
        p_low = sum(math.exp(-lam + i * math.log(lam) - math.lgamma(i + 1)) for i in range(k + 1))
        pr(f"family {fm} {nm}: genomes {len({r['genome'] for r in S2})} n={len(S2)} CTT {k}/{lam:.1f} = {k / lam:.2f} (cgn; Poisson P(<=obs) {p_low:.3g}); ggn {ratio(oo2, ee2, 'CTT'):.2f}")

# family-link decoy: reversed peptides keep composition, break position
pep = {r["site_id"]: r["pep_up"] for r in A if len(r["pep_up"]) == 20}
ids = list(pep)
rng = random.Random(5)
def ident(a, b): return sum(1 for x, y in zip(a, b) if x == y)
real = collections.Counter(); dec = collections.Counter()
for _ in range(200000):
    a, b = rng.choice(ids), rng.choice(ids)
    if a == b: continue
    real[ident(pep[a], pep[b])] += 1
    dec[ident(pep[a], pep[b][::-1])] += 1
for t in (6, 8, 10, 12):
    pr(f"pairs with >= {t}/20 identical: random real pairs {sum(v for k, v in real.items() if k >= t)}, "
       f"reversed-peptide decoys {sum(v for k, v in dec.items() if k >= t)} (of 200,000)")

def families(th):
    idx = collections.defaultdict(list)
    for sid, s in pep.items():
        for i in range(18):
            w = s[i:i + 3]
            if "*" in w or "X" in w: continue
            idx[(i, w)].append(sid)
    par = {s: s for s in pep}
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for key, lst in idx.items():
        for a, b in itertools.combinations(lst, 2):
            if ident(pep[a], pep[b]) >= th:
                ra, rb = find(a), find(b)
                if ra != rb: par[ra] = rb
    return {s: find(s) for s in pep}

for th in (8, 10, 12):
    fam = families(th)
    out = []
    for nm, S in (("skip", SK), ("inframe", IF)):
        S2 = [r for r in S if r["site_id"] in fam]
        size = collections.Counter(fam[r["site_id"]] for r in S2)
        w = {r["site_id"]: 1.0 / size[fam[r["site_id"]]] for r in S2}
        og, eg = obs_exp(S2, "ggn", weight=w)
        oc, ec = obs_exp(S2, "cgn", weight=w)
        out.append(f"{nm}: {len(size)} families, largest {max(size.values())}, CTT ggn {ratio(og, eg, 'CTT'):.2f} cgn {ratio(oc, ec, 'CTT'):.2f}")
    pr(f"families at >= {th}/20: " + " | ".join(out))
    if th == 10:
        # one site per family, family-level bootstrap of the CTT ratio (skip-like, ggn bg)
        S2 = [r for r in SK if r["site_id"] in fam]
        byf = collections.defaultdict(list)
        for r in S2: byf[fam[r["site_id"]]].append(r)
        units = []
        for f, rs in byf.items():
            n = len(rs)
            units.append([sum(r["up"] == "CTT" for r in rs) / n, sum(BG["ggn"][r["genome"]]["CTT"] for r in rs) / n])
        lo, hi, _ = boot(units)
        pr(f"  family-weighted bootstrap over {len(units)} skip-like families (>=10/20, ggn bg, seed 1): {lo:.2f} - {hi:.2f}")
# exact duplicate upstream peptides among CTT skip-like sites
c = collections.Counter(r["pep_up"] for r in SK if r["up"] == "CTT" and r["pep_up"])
pr("CTT skip-like sites sharing an identical 20-aa upstream peptide: groups", sum(1 for v in c.values() if v > 1),
   "sites in them", sum(v for v in c.values() if v > 1))
open(OUT + "analysis2.txt", "w").write("\n".join(P) + "\n")
# codons ranked by skip-like ratio under the context-matched background (obs >= 20)
o, e = obs_exp(SK, "cgn"); o2, e2 = obs_exp(IF, "cgn")
rk = sorted([c for c in CODONS if o[c] >= 20], key=lambda c: -ratio(o, e, c))
pr("top 8 skip-like codons, cgn bg, obs >= 20 (obs/exp=ratio | in-frame ratio): " +
   "; ".join(f"{c} {o[c]:.0f}/{e[c]:.1f}={ratio(o, e, c):.2f} | {ratio(o2, e2, c):.2f}" for c in rk[:8]))
pr(f"CTT share of sites: skip {sum(r['up'] == 'CTT' for r in SK)}/{len(SK)}, in-frame {sum(r['up'] == 'CTT' for r in IF)}/{len(IF)}")
open(OUT + "analysis2.txt", "w").write("\n".join(P) + "\n")
