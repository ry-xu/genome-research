"""Pass 2: rule numbers (P1-P3, decision, bootstrap) and break attempts, from out/ tables of recheck.py."""
import csv, json, random, collections, itertools, sys

BASE = ""
OUT = BASE + "results/skip_context_chicken/recheck/out/"
SEL = BASE + "results/rna_chicken/genomes_selected.tsv"
SITES = BASE + "results/rna_chicken/sites.tsv"
B = "TCAG"
CODONS = [a + b + c for a in B for b in B for c in B]
STOP = {"TAA", "TAG", "TGA"}
ENRICHED = ["CTT", "TTG", "TCA", "CCA", "GAC"]
DEPLETED = ["TTT", "AAA", "GAA", "TAT"]
LEU = ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"]

chk = json.load(open(OUT + "checks.json"))
excl = set(chk["excluded"])
bgraw = json.load(open(OUT + "background.json"))
usage = json.load(open(OUT + "cds_usage.json"))
genus = {}
for r in csv.DictReader(open(SEL), delimiter="\t"):
    genus[r["acc"]] = r["genus"]
crun = {}
for r in csv.DictReader(open(SITES), delimiter="\t"):
    if r["set"] == "A":
        crun[r["site_id"]] = int(r["run_end"]) - int(r["run_start"])


def freqs(name):
    out = {}
    for g, d in bgraw.items():
        c = d[name]
        t = sum(c.values())
        out[g] = {k: c.get(k, 0) / t for k in CODONS} if t else None
    return out


BG = {k: freqs(k) for k in ("ggn", "cgg", "cgn", "all")}

rows = list(csv.DictReader(open(OUT + "sites_up.tsv"), delimiter="\t"))
A = [r for r in rows if r["set"] == "A" and r["genome"] not in excl]
Aall = [r for r in rows if r["set"] == "A"]
C = [r for r in rows if r["set"] == "C" and r["genome"] not in excl]
RND = [r for r in csv.DictReader(open(OUT + "random_up.tsv"), delimiter="\t") if r["genome"] not in excl]

P = []


def pr(*a):
    s = " ".join(str(x) for x in a)
    P.append(s)
    print(s)


def obs_exp(sites, bgname="ggn", codon_key="up", weight=None):
    """obs[c], exp[c]; exp = sum over genomes of freq_g(c) * n sites in g (per-site sum, same thing)."""
    obs = collections.Counter()
    exp = collections.Counter()
    f = BG[bgname]
    for r in sites:
        w = 1.0 if weight is None else weight[r["site_id"]]
        obs[r[codon_key]] += w
        fg = f[r["genome"]]
        for c in CODONS:
            exp[c] += w * fg[c]
    return obs, exp


def ratio(o, e, c):
    return o[c] / e[c] if e[c] > 0 else float("nan")


# ---------------- rule numbers ----------------
gens_all = sorted({r["genome"] for r in Aall})
gens = sorted({r["genome"] for r in A})
pr("== Inputs ==")
pr("Set A genomes:", len(gens_all), "| excluded (on cgg_test_T list):", len(excl), sorted(excl), "| used:", len(gens))
pr("site codons reproduced:", chk["checks"].get("codon_ok"), "bad:", chk["checks"].get("codon_BAD", 0),
   "| Set A class + both run lengths reproduced:", chk["checks"].get("A_class_match"),
   "mismatch:", chk["checks"].get("A_class_MISMATCH", 0), "| internal stops in CDS calls:", chk["checks"]["internal_stops"])
cls_n = collections.Counter(r["given_class"] for r in A)
pr("Set A sites used by class:", dict(cls_n))
bad_up = [r for r in A if len(r["up"]) != 3 or any(ch not in "ACGT" for ch in r["up"])]
pr("Set A sites with no clean upstream codon:", len(bad_up), "| upstream codon is a stop:", sum(r["up"] in STOP for r in A))
pr("Set A sites in-frame inside their called CDS:", collections.Counter((r["given_class"], r["in_cds"]) for r in A))
pr("upstream codon inside the called CDS (codon_idx>=1):",
   collections.Counter((r["given_class"], int(r["in_cds"] == "1" and int(r["codon_idx"]) >= 1)) for r in A))

SK = [r for r in A if r["given_class"] == "skip"]
IF = [r for r in A if r["given_class"] == "inframe"]
os_, es_ = obs_exp(SK)
oi, ei = obs_exp(IF)
pr("\n== All codons, observed / expected (background: codon before in-frame GGN, per genome) ==")
pr("codon aa  skip_obs skip_exp skip_ratio | inframe_obs inframe_exp inframe_ratio")
table = []
for c in sorted(CODONS, key=lambda c: -ratio(os_, es_, c) if es_[c] else 0):
    table.append((c, os_[c], es_[c], ratio(os_, es_, c), oi[c], ei[c], ratio(oi, ei, c)))
    pr(f"{c} {'':2} {os_[c]:6.0f} {es_[c]:8.1f} {ratio(os_, es_, c):6.2f} | {oi[c]:6.0f} {ei[c]:8.1f} {ratio(oi, ei, c):6.2f}")
with open(OUT + "ratios_recheck.tsv", "w") as fh:
    fh.write("codon\tskip_obs\tskip_exp\tskip_ratio\tinframe_obs\tinframe_exp\tinframe_ratio\n")
    for t in table:
        fh.write("\t".join([t[0], f"{t[1]:.0f}", f"{t[2]:.3f}", f"{t[3]:.4f}", f"{t[4]:.0f}", f"{t[5]:.3f}", f"{t[6]:.4f}"]) + "\n")

r_ctt_s = ratio(os_, es_, "CTT")
r_ctt_i = ratio(oi, ei, "CTT")
p1 = r_ctt_s >= 2
p2 = r_ctt_s > r_ctt_i
p3n = sum(ratio(os_, es_, c) > 1 for c in ENRICHED)
p3 = p3n >= 4
pr("\n== Predictions ==")
pr(f"P1 CTT skip-like: {os_['CTT']:.0f} obs / {es_['CTT']:.2f} exp = {r_ctt_s:.3f} (>= 2?) -> {'HOLDS' if p1 else 'FAILS'}")
pr(f"P2 CTT in-frame-like: {oi['CTT']:.0f} / {ei['CTT']:.2f} = {r_ctt_i:.3f}; skip {r_ctt_s:.3f} > in-frame? -> {'HOLDS' if p2 else 'FAILS'}")
pr("P3 skip-like ratios:", {c: round(ratio(os_, es_, c), 3) for c in ENRICHED}, f"-> {p3n}/5 > 1 ->", "HOLDS" if p3 else "FAILS")
pr("Depleted four, skip-like:", {c: (int(os_[c]), round(es_[c], 1), round(ratio(os_, es_, c), 3)) for c in DEPLETED})
pr("Depleted four, in-frame-like:", {c: (int(oi[c]), round(ei[c], 1), round(ratio(oi, ei, c), 3)) for c in DEPLETED})
pr("Enriched five, in-frame-like:", {c: round(ratio(oi, ei, c), 3) for c in ENRICHED})
dec = "REPLICATED" if (p1 and p2 and p3) else ("NOT REPLICATED" if not p1 else "neither (P1 holds, P2 or P3 fails; rule silent)")
pr("Decision:", dec)

# bootstrap over genomes
per_g = collections.defaultdict(lambda: [0, 0.0])
for r in SK:
    per_g[r["genome"]][0] += r["up"] == "CTT"
    per_g[r["genome"]][1] += BG["ggn"][r["genome"]]["CTT"]


def boot(units, n=1000, seed=1):
    rng = random.Random(seed)
    vals = []
    for _ in range(n):
        o = e = 0.0
        for _ in range(len(units)):
            u = units[rng.randrange(len(units))]
            o += u[0]
            e += u[1]
        vals.append(o / e)
    vals.sort()
    return vals[int(0.025 * n)], vals[int(0.975 * n) - 1], vals


gsk = sorted(per_g)
lo, hi, _ = boot([per_g[g] for g in gsk])
pr(f"\nBootstrap over the {len(gsk)} genomes with >= 1 skip-like site (1,000, random.Random(1)): CTT ratio 95% "
   f"percentile interval {lo:.3f} - {hi:.3f}")
units_all = [per_g[g] if g in per_g else [0, 0.0] for g in gens]
lo2, hi2, _ = boot(units_all)
pr(f"Same over all {len(gens)} used genomes (zero-site genomes included): {lo2:.3f} - {hi2:.3f}")

# ---------------- break (a): genomes, genera, genes, families ----------------
pr("\n== Break (a): is the CTT excess carried by a few genomes, genera or genes? ==")
pg = []
for g in gsk:
    o, e = per_g[g]
    n = sum(1 for r in SK if r["genome"] == g)
    pg.append((o - e, g, int(o), e, n))
pg.sort(reverse=True)
pr("genomes with >= 1 CTT skip-like site:", sum(1 for x in pg if x[2] > 0), "of", len(pg))
pr("top 10 genomes by CTT excess (obs - exp): genome genus n_skip obs exp ratio")
for x in pg[:10]:
    pr(f"  {x[1]} {genus.get(x[1], '?')} {x[4]} {x[2]} {x[3]:.2f} {x[2] / x[3]:.2f}")
tot_o = sum(x[2] for x in pg)
pr(f"share of CTT skip-like sites in the top 5 genomes by excess: {sum(x[2] for x in pg[:5])}/{tot_o}")
for label, key in (("excess", lambda x: x[0]), ("CTT count", lambda x: x[2]), ("skip-like site count", lambda x: x[4])):
    top = {x[1] for x in sorted(pg, key=key, reverse=True)[:5]}
    sub = [r for r in SK if r["genome"] not in top]
    o, e = obs_exp(sub)
    pr(f"drop top 5 by {label}: CTT {o['CTT']:.0f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f} over {len(sub)} skip-like sites")
for k in (10, 20):
    top = {x[1] for x in pg[:k]}
    sub = [r for r in SK if r["genome"] not in top]
    o, e = obs_exp(sub)
    pr(f"drop top {k} by excess: CTT {o['CTT']:.0f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f}")
loo = []
for g in gsk:
    sub = [r for r in SK if r["genome"] != g]
    o, e = obs_exp(sub)
    loo.append(ratio(o, e, "CTT"))
pr(f"leave-one-genome-out CTT ratio range: {min(loo):.3f} - {max(loo):.3f}")
# per-genome ratio distribution among genomes with expected >= 1
big = [(x[2] / x[3]) for x in pg if x[3] >= 1]
big.sort()
pr(f"per-genome CTT ratio for genomes with exp >= 1 ({len(big)}): median {big[len(big) // 2]:.2f}, "
   f"share > 1: {sum(b > 1 for b in big)}/{len(big)}, share >= 2: {sum(b >= 2 for b in big)}/{len(big)}")
# genus
pgen = collections.defaultdict(lambda: [0, 0.0, 0])
for r in SK:
    gg = genus.get(r["genome"], "?")
    pgen[gg][0] += r["up"] == "CTT"
    pgen[gg][1] += BG["ggn"][r["genome"]]["CTT"]
    pgen[gg][2] += 1
pr("by genus (n_skip, CTT obs, exp, ratio), genera with >= 100 skip-like sites:")
for gg, (o, e, n) in sorted(pgen.items(), key=lambda kv: -kv[1][2]):
    if n >= 100:
        pr(f"  {gg}: {n} {o} {e:.1f} {o / e:.2f}")
rest = [(o, e, n) for gg, (o, e, n) in pgen.items() if n < 100]
pr(f"  other genera pooled: {sum(x[2] for x in rest)} {sum(x[0] for x in rest)} {sum(x[1] for x in rest):.1f} "
   f"{sum(x[0] for x in rest) / sum(x[1] for x in rest):.2f}")
# genes
per_cds = collections.Counter(r["cds"] for r in SK if r["up"] == "CTT")
pr("CTT skip-like sites per CDS:", collections.Counter(per_cds.values()))
one = {}
for r in SK:
    one.setdefault(r["cds"], []).append(r)
rng = random.Random(3)
sub = [rng.choice(v) for v in one.values()]
o, e = obs_exp(sub)
pr(f"one random skip-like site per CDS ({len(sub)} sites): CTT {o['CTT']:.0f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f}")

# families: anchored upstream peptide (codons -21..-2, the -1 codon left out), single linkage at >= 8/20 identical
pep = {r["site_id"]: r["pep_up"] for r in A if len(r["pep_up"]) == 20}
idx = collections.defaultdict(list)
for sid, s in pep.items():
    for i in range(18):
        w = s[i:i + 3]
        if "*" in w or "X" in w:
            continue
        idx[(i, w)].append(sid)
par = {s: s for s in pep}


def find(x):
    while par[x] != x:
        par[x] = par[par[x]]
        x = par[x]
    return x


def ident(a, b):
    return sum(1 for x, y in zip(a, b) if x == y)


seen = set()
for key, lst in idx.items():
    if len(lst) > 2000:
        continue
    for a, b in itertools.combinations(lst, 2):
        if (a, b) in seen:
            continue
        seen.add((a, b))
        if ident(pep[a], pep[b]) >= 8:
            ra, rb = find(a), find(b)
            if ra != rb:
                par[ra] = rb
# decoy calibration: random pairs of sites
ids = list(pep)
rng = random.Random(5)
dec_hits = sum(ident(pep[rng.choice(ids)], pep[rng.choice(ids)]) >= 8 for _ in range(200000))
pr(f"\nFamilies (anchored 20-aa upstream peptide, >= 8/20 identical, single linkage): decoy random pairs "
   f"reaching the threshold: {dec_hits}/200000")
fam = {s: find(s) for s in pep}
for clsname, S in (("skip", SK), ("inframe", IF)):
    S2 = [r for r in S if r["site_id"] in fam]
    size = collections.Counter(fam[r["site_id"]] for r in S2)
    w = {r["site_id"]: 1.0 / size[fam[r["site_id"]]] for r in S2}
    o, e = obs_exp(S2, weight=w)
    ctt_f = collections.Counter(fam[r["site_id"]] for r in S2 if r["up"] == "CTT")
    pr(f"{clsname}: {len(S2)} sites in {len(size)} families (largest {max(size.values())}); CTT sites in "
       f"{len(ctt_f)} families (largest {max(ctt_f.values()) if ctt_f else 0}); family-weighted CTT "
       f"{o['CTT']:.1f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f}")
    if clsname == "skip":
        top = ctt_f.most_common(5)
        pr("  largest CTT families (size, genera):",
           [(n, sorted({genus.get(r['genome']) for r in S2 if fam[r['site_id']] == f})) for f, n in top])

# ---------------- break (b): does the class rule make the bias by itself? ----------------
pr("\n== Break (b): the class rule applied to other sites ==")
for nm, S in (("Set C skip-like (CGA/CGC/CGT, same walks, same rule, my classes)", [r for r in C if r["my_class"] == "skip"]),
              ("Set C in-frame-like", [r for r in C if r["my_class"] == "inframe"]),
              ("Set C all", C)):
    for bgn in ("ggn", "cgn"):
        o, e = obs_exp(S, bgn)
        pr(f"{nm} [bg {bgn}]: n={len(S)} CTT {o['CTT']:.0f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f}; "
           + " ".join(f"{c}={ratio(o, e, c):.2f}" for c in ENRICHED[1:] + DEPLETED))
for cg in ("CGA", "CGC", "CGT"):
    S = [r for r in C if r["my_class"] == "skip" and r["codon"] == cg]
    o, e = obs_exp(S, "cgn")
    pr(f"  Set C skip-like {cg}: n={len(S)} CTT {o['CTT']:.0f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f}")
# random in-frame CDS positions classed by the same rule, background = all in-frame codons of the genome
for nm in ("skip", "inframe", "unclear"):
    S = [r for r in RND if r["my_class"] == nm]
    o, e = obs_exp(S, "all")
    pr(f"random CDS positions {nm}-like: n={len(S)} CTT {o['CTT']:.0f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f}")
S = [r for r in RND if r["my_class"] == "skip" and r["codon"].startswith("GG")]
o, e = obs_exp(S, "ggn")
pr(f"random positions at an in-frame GGN, skip-like (bg ggn): n={len(S)} CTT {o['CTT']:.0f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f}")
# random positions matched to skip-like geometry: within 30 codons of the CDS end
S = [r for r in RND if 0 <= int(r["codons_to_cds_end"]) <= 30]
o, e = obs_exp(S, "all")
pr(f"random positions within 30 codons of the called stop: n={len(S)} CTT {o['CTT']:.0f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f}")

# ---------------- break (c): other artefacts ----------------
pr("\n== Break (c): other backgrounds, synonyms, nucleotides, geometry ==")
for bgn, lab in (("cgg", "codon before in-frame CGG in CDS calls"), ("cgn", "codon before in-frame CGA/CGC/CGT"),
                 ("all", "all in-frame codons")):
    o, e = obs_exp(SK, bgn)
    o2, e2 = obs_exp(IF, bgn)
    pr(f"bg = {lab}: CTT skip {ratio(o, e, 'CTT'):.3f}, in-frame {ratio(o2, e2, 'CTT'):.3f}; P3 count "
       f"{sum(ratio(o, e, c) > 1 for c in ENRICHED)}/5")
# own-gene usage background
for nm, S in (("skip", SK), ("inframe", IF)):
    o = sum(r["up"] == "CTT" for r in S)
    e = 0.0
    n = 0
    for r in S:
        u = usage.get(r["cds"])
        if not u:
            continue
        t = sum(u.values())
        e += u.get("CTT", 0) / t
        n += 1
    pr(f"own-CDS codon usage as background, {nm}: CTT {o}/{e:.2f} = {o / e:.3f} (n={n})")
# synonymous share among Leu
for nm, S in (("skip", SK), ("inframe", IF)):
    o, e = obs_exp(S)
    lo_ = sum(o[c] for c in LEU)
    le_ = sum(e[c] for c in LEU)
    pr(f"{nm}: Leu codons obs/exp {lo_:.0f}/{le_:.1f} = {lo_ / le_:.3f}; CTT share of Leu obs {o['CTT'] / lo_:.3f} "
       f"vs exp {e['CTT'] / le_:.3f}; per Leu codon: " + " ".join(f"{c}={ratio(o, e, c):.2f}" for c in LEU))
# nucleotide by position
for nm, S in (("skip", SK), ("inframe", IF)):
    o, e = obs_exp(S)
    parts = []
    for pos in range(3):
        for nt in "ACGT":
            on = sum(o[c] for c in CODONS if c[pos] == nt)
            en = sum(e[c] for c in CODONS if c[pos] == nt)
            parts.append(f"{pos + 1}{nt}={on / en:.2f}")
    pr(f"{nm} upstream-codon nucleotide ratios: " + " ".join(parts))
    xtt = {c: round(ratio(o, e, c), 2) for c in ("ATT", "CTT", "GTT", "TTT")}
    ctx = {c: round(ratio(o, e, c), 2) for c in ("CTA", "CTC", "CTG", "CTT")}
    pr(f"  NTT: {xtt}  CTN: {ctx}")
# -2 codon
for nm, S in (("skip", SK), ("inframe", IF)):
    S2 = [r for r in S if len(r["up2"]) == 3]
    o, e = obs_exp(S2, "all", codon_key="up2")
    pr(f"{nm}: codon -2 (two before the CGG, bg all codons) CTT {ratio(o, e, 'CTT'):.2f}; top 5: "
       + " ".join(f"{c}={ratio(o, e, c):.2f}" for c in sorted(CODONS, key=lambda c: -ratio(o, e, c) if e[c] > 5 else 0)[:5]))
# geometry
for nm, f in (("inside called CDS, in frame, >= 30 codons from start", lambda r: r["in_cds"] == "1" and int(r["codon_idx"]) >= 30),
              ("inside called CDS, in frame, < 30 codons from start", lambda r: r["in_cds"] == "1" and int(r["codon_idx"]) < 30),
              ("not in frame in its called CDS (past an earlier skip or beyond the call)", lambda r: r["in_cds"] == "0"),
              ("C not in a C run (run length 1)", lambda r: crun[r["site_id"]] == 1),
              ("+1 frame closed upstream within 10 codons", lambda r: int(r["up_plus1_open"]) < 10),
              ("+1 frame open >= 60 codons upstream", lambda r: int(r["up_plus1_open"]) >= 60)):
    S = [r for r in SK if f(r)]
    S2 = [r for r in IF if f(r)]
    o, e = obs_exp(S)
    o2, e2 = obs_exp(S2)
    pr(f"{nm}: skip n={len(S)} CTT {o['CTT']:.0f}/{e['CTT']:.1f} = {ratio(o, e, 'CTT'):.2f} | in-frame n={len(S2)} "
       f"{ratio(o2, e2, 'CTT'):.2f}")
# old-frame stop distance after the CGG (skip-like)
for lo_, hi_ in ((0, 5), (6, 15), (16, 30)):
    S = [r for r in SK if lo_ <= int(r["old"]) <= hi_]
    o, e = obs_exp(S)
    pr(f"skip-like with old-frame stop {lo_}-{hi_} codons after: n={len(S)} CTT {ratio(o, e, 'CTT'):.2f}")
# excluded genomes, for reference
SKx = [r for r in Aall if r["genome"] in excl and r["given_class"] == "skip"]
o, e = obs_exp(SKx)
pr(f"\n(reference) the {len(excl)} excluded genomes, skip-like: n={len(SKx)} CTT {o['CTT']:.0f}/{e['CTT']:.2f} = {ratio(o, e, 'CTT'):.3f}")
open(OUT + "analysis.txt", "w").write("\n".join(P) + "\n")
