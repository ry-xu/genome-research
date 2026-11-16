"""Score code-walls groups (RULES.md). Reads per-phage stats (codon_stats.py) and run_ids.tsv.

usage: python -I score.py RUN_IDS.tsv STATS.tsv[,STATS2.tsv...] MODE [--target T1] [--ref R] [--perm 2000] [--seed 1]
MODE:
  groups   target vs ref: GC-bin-stratified Mantel-Haenszel rate ratio of O/E (and of skip/O), within-bin
           permutation p, phage-genus cluster bootstrap CI, per-phage classes.
  null     null world: 200 fake targets drawn from ref within bins (target's bin sizes), permutation p each.
  recall   P3 at group level: fake targets (target's bin sizes) drawn from ref, their stats taken from --plantstats
           (CGG removed at share f); reference = the remaining ref phages, unplanted. Share with p < 0.01 and RR.
  heldout  N control: k from a random half of ref, avoid-call share in the other half.
Per-phage avoid call: E >= 5, Poisson P(X <= O | k_bin * E) < 0.01, rr < 0.25 (k_bin = pooled ref O/E in the bin).
Host-style: skip >= 3 and s = skip/O above the ref 99th percentile (phages with O >= 1).
"""
import argparse
import csv
import math
import random
from collections import defaultdict

import numpy as np


def pois_cdf(o, mu):
    if mu <= 0:
        return 1.0
    t = math.exp(-mu)
    s = t
    for i in range(1, o + 1):
        t *= mu / i
        s += t
    return min(1.0, s)


def load(run_ids, stats_files, meta=None):
    sets, binof, host = {}, {}, {}
    for r in csv.DictReader(open(run_ids), delimiter="\t"):
        sets[r["id"]] = r["set"]; binof[r["id"]] = int(r["gc_bin"]); host[r["id"]] = r["host"]
    ph = {}
    for f in stats_files.split(","):
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r["id"] in sets:
                ph[r["id"]] = dict(O=int(r["O"]), E=float(r["E"]), skip=int(r["skip_sites"]),
                                   ccg=int(r.get("trna_ccg", 0) or 0), tcg=int(r.get("trna_tcg", 0) or 0),
                                   set=sets[r["id"]], bin=binof[r["id"]],
                                   **{k: float(r[k]) if r.get(k) not in (None, "", "NA") else 0.0
                                      for k in ("O_CGA", "E_CGA", "O_CGT", "E_CGT", "O_CGC", "E_CGC")},
                                   genus=(host[r["id"]].split("g__")[1].split(";")[0] if "g__" in host[r["id"]] else ""))
    return ph


def pclust(path):
    """vOTU genome id -> UHGV vGENUS (from uhgv_taxonomy) (cluster for bootstrap); falls back to the id itself."""
    import subprocess
    out = {}
    p = subprocess.run(["zstd", "-dcq", path], capture_output=True, text=True).stdout.splitlines()
    ix = {k: i for i, k in enumerate(p[0].split("\t"))}
    for line in p[1:]:
        f = line.split("\t")
        g = [x for x in f[ix["uhgv_taxonomy"]].split(";") if x.startswith("vGENUS-")]
        out[f[0]] = g[0] if g else f[0]
    return out


def mh_rr(Ot, Et, Or, Er):
    """Mantel-Haenszel rate ratio over bins; inputs are per-bin arrays (E plays person-time)."""
    tot = Et + Er
    ok = (Et > 0) & (Er > 0)
    num = np.sum((Ot * Er / np.where(tot > 0, tot, 1))[ok])
    den = np.sum((Or * Et / np.where(tot > 0, tot, 1))[ok])
    return num / den if den > 0 else float("nan")


def bin_sums(ids, ph, bins, key_o, key_e):
    O = np.zeros(len(bins)); E = np.zeros(len(bins)); bi = {b: i for i, b in enumerate(bins)}
    for i in ids:
        p = ph[i]
        O[bi[p["bin"]]] += p[key_o]; E[bi[p["bin"]]] += p[key_e]
    return O, E


def group_stat(tids, rids, ph, bins):
    Ot, Et = bin_sums(tids, ph, bins, "O", "E"); Or, Er = bin_sums(rids, ph, bins, "O", "E")
    St, Ut = bin_sums(tids, ph, bins, "skip", "O"); Sr, Ur = bin_sums(rids, ph, bins, "skip", "O")
    return mh_rr(Ot, Et, Or, Er), mh_rr(St, Ut, Sr, Ur)


def perm_p(tids, rids, ph, bins, nperm, rng):
    """Two-sided within-bin label permutation p for log MH RR (O/E) and (skip/O); vectorized over permutations."""
    obs = group_stat(tids, rids, ph, bins)
    nrs = np.random.default_rng(rng.randrange(2**32))
    tset = set(tids)
    byb = defaultdict(list)
    for i in list(tids) + list(rids):
        byb[ph[i]["bin"]].append(i)
    num = [np.zeros(nperm), np.zeros(nperm)]; den = [np.zeros(nperm), np.zeros(nperm)]
    for b, lst in byb.items():
        nt = sum(1 for i in lst if i in tset)
        if nt == 0 or nt == len(lst):
            continue
        O = np.array([ph[i]["O"] for i in lst], float); E = np.array([ph[i]["E"] for i in lst], float)
        S = np.array([ph[i]["skip"] for i in lst], float)
        idx = np.argsort(nrs.random((nperm, len(lst))), axis=1)[:, :nt]
        Ot, Et, St = O[idx].sum(1), E[idx].sum(1), S[idx].sum(1)
        Oa, Ea, Sa = O.sum(), E.sum(), S.sum()
        for j, (xt, yt, xa, ya) in enumerate(((Ot, Et, Oa, Ea), (St, Ot, Sa, Oa))):
            xr, yr = xa - xt, ya - yt
            ok = (yt > 0) & (yr > 0)
            num[j] += np.where(ok, xt * yr / ya if ya > 0 else 0, 0)
            den[j] += np.where(ok, xr * yt / ya if ya > 0 else 0, 0)
    ps = []
    for j in range(2):
        with np.errstate(divide="ignore", invalid="ignore"):
            st = num[j] / den[j]
        if math.isnan(obs[j]):
            ps.append(float("nan")); continue
        if obs[j] <= 0:  # RR = 0: one-sided extreme; p = share of permutations as extreme
            ps.append((np.sum(np.nan_to_num(st, nan=1.0) <= 0) + 1) / (nperm + 1)); continue
        lo = np.abs(np.log(np.where(st > 0, st, 1e-9)))
        ps.append((np.sum(lo >= abs(math.log(obs[j])) - 1e-12) + 1) / (nperm + 1))
    return obs, ps


def cluster_boot(tids, rids, ph, bins, clus, nboot, rng):
    ct = defaultdict(list); cr = defaultdict(list)
    for i in tids:
        ct[clus.get(i, i)].append(i)
    for i in rids:
        cr[clus.get(i, i)].append(i)
    kt, kr = list(ct), list(cr)
    res = []
    for _ in range(nboot):
        bt = [i for k in rng.choices(kt, k=len(kt)) for i in ct[k]]
        br = [i for k in rng.choices(kr, k=len(kr)) for i in cr[k]]
        res.append(group_stat(bt, br, ph, bins))
    a = np.array(res, dtype=float)
    return [tuple(np.nanpercentile(a[:, j], [2.5, 97.5])) for j in range(2)]


def classes(ids, ph, kbin, s99):
    c = defaultdict(int)
    for i in ids:
        p = ph[i]; k = kbin.get(p["bin"], float("nan"))
        s = p["skip"] / p["O"] if p["O"] else 0.0
        host = p["skip"] >= 3 and s > s99
        ok = p["E"] >= 5 and k == k and k > 0  # A3: no class from a bin whose reference O/E is 0 or missing
        av = ok and pois_cdf(p["O"], k * p["E"]) < 0.01 and (p["O"] / p["E"]) / k < 0.25
        anc = ok and (p["O"] / p["E"]) / k >= 0.5 and not host
        c["host-style" if host else "avoid" if av else "ancestral" if anc else "unresolved"] += 1
        c["E>=5"] += p["E"] >= 5
        c["tRNA_CCG"] += p["ccg"] > 0
        c["tRNA_TCG"] += p["tcg"] > 0
    return dict(c)


def kfit(rids, ph):
    O = defaultdict(float); E = defaultdict(float)
    for i in rids:
        O[ph[i]["bin"]] += ph[i]["O"]; E[ph[i]["bin"]] += ph[i]["E"]
    return {b: (O[b] / E[b] if E[b] > 0 else float("nan")) for b in E}


def s99_of(rids, ph):
    v = [ph[i]["skip"] / ph[i]["O"] for i in rids if ph[i]["O"] >= 1]
    return float(np.percentile(v, 99)) if v else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_ids"); ap.add_argument("stats"); ap.add_argument("mode")
    ap.add_argument("--target", default="T1"); ap.add_argument("--ref", default="R")
    ap.add_argument("--perm", type=int, default=2000); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--nnull", type=int, default=200)
    ap.add_argument("--plantstats", default=None)
    ap.add_argument("--logo", action="store_true", help="also leave out each target host genus in turn")
    ap.add_argument("--meta", default="data/uhgv/votus_metadata_extended.tsv.zst")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    ph = load(a.run_ids, a.stats)
    refs = set(a.ref.split(","))
    rids_all = [i for i, p in ph.items() if p["set"] in refs]
    if a.mode == "heldout":
        r = rids_all[:]; rng.shuffle(r); h = len(r) // 2
        kb = kfit(r[:h], ph); s99 = s99_of(r[:h], ph)
        c = classes(r[h:], ph, kb, s99)
        n5 = c["E>=5"]
        print(f"N heldout: n={len(r) - h} E>=5={n5} avoid={c.get('avoid', 0)} share_of_E5={c.get('avoid', 0) / max(1, n5):.4f}"
              f" host-style={c.get('host-style', 0)} classes={c}")
        return
    if a.mode in ("null", "recall"):
        if a.mode == "recall":
            pl = load(a.run_ids, a.plantstats)
            for i in list(pl):
                ph["PL_" + i] = pl[i]
        # fake targets: the target set's per-bin sizes, drawn from ref; the rest of ref is the reference
        tset = set(a.target.split(","))
        tsz = defaultdict(int)  # target bin sizes from run_ids (target stats are not read here)
        for r in csv.DictReader(open(a.run_ids), delimiter="\t"):
            if r["set"] in tset:
                tsz[int(r["gc_bin"])] += 1
        byb = defaultdict(list)
        for i in rids_all:
            byb[ph[i]["bin"]].append(i)
        bins = sorted(set(tsz))
        ps, rrs = [], []
        for k in range(a.nnull):
            ft, fr = [], []
            for b in bins:
                s = byb[b][:]; rng.shuffle(s); n = min(tsz[b], len(s) // 2); ft += s[:n]; fr += s[n:]
            if a.mode == "recall":
                ft = ["PL_" + i for i in ft if "PL_" + i in ph]
            obs, p = perm_p(ft, fr, ph, bins, a.perm, rng)
            ps.append(p); rrs.append(obs[0])
        ps = np.array(ps)
        print(f"{a.mode}: median MH RR O/E = {np.nanmedian(rrs):.3f}")
        for j, nm in enumerate(["O/E", "skip/O"]):
            print(f"null {nm}: n={len(ps)} share p<0.05 = {np.mean(ps[:, j] < 0.05):.3f}; share p<0.01 = "
                  f"{np.mean(ps[:, j] < 0.01):.3f}; median p = {np.median(ps[:, j]):.3f}")
        return
    if a.mode == "groups":
        tset = set(a.target.split(","))
        tids = [i for i, p in ph.items() if p["set"] in tset]
        bins = sorted({ph[i]["bin"] for i in tids})
        rids = [i for i in rids_all if ph[i]["bin"] in bins]
        kb = kfit(rids, ph); s99 = s99_of(rids, ph)
        obs, p = perm_p(tids, rids, ph, bins, a.perm, rng)
        clus = pclust(a.meta)
        ci = cluster_boot(tids, rids, ph, bins, clus, 1000, rng)
        nclt = len({clus.get(i, i) for i in tids})
        print(f"target {a.target}: n={len(tids)} phage-genus clusters={nclt}; ref {a.ref}: n={len(rids)}; bins={bins}")
        print(f"  sum O={sum(ph[i]['O'] for i in tids)} sum E={sum(ph[i]['E'] for i in tids):.1f}; "
              f"ref sum O={sum(ph[i]['O'] for i in rids)} sum E={sum(ph[i]['E'] for i in rids):.1f}")
        print(f"  MH rate ratio O/E (target/ref) = {obs[0]:.3f}  perm p = {p[0]:.4f}  cluster-boot 95% CI = "
              f"{ci[0][0]:.3f}-{ci[0][1]:.3f}")
        print(f"  MH ratio skip/O (target/ref)   = {obs[1]:.3f}  perm p = {p[1]:.4f}  cluster-boot 95% CI = "
              f"{ci[1][0]:.3f}-{ci[1][1]:.3f}")
        for x in ("CGA", "CGT", "CGC"):
            Ot, Et = bin_sums(tids, ph, bins, "O_" + x, "E_" + x); Or, Er = bin_sums(rids, ph, bins, "O_" + x, "E_" + x)
            rx = mh_rr(Ot, Et, Or, Er)
            print(f"  sibling {x}: MH RR (target/ref) = {rx:.3f}; RR_use(CGG)/RR({x}) = {obs[0] / rx:.3f}")
        print(f"  classes target: {classes(tids, ph, kb, s99)}")
        print(f"  classes ref:    {classes(rids, ph, kb, s99)}   (s99={s99:.3f})")
        if a.logo:
            for g in sorted({ph[i]["genus"] for i in tids}):
                sub = [i for i in tids if ph[i]["genus"] != g]
                o2, p2 = perm_p(sub, rids, ph, bins, a.perm, rng)
                print(f"  leave out {g or '(no genus)'} (n={len(tids) - len(sub)}): RR_use={o2[0]:.3f} p={p2[0]:.4f}  "
                      f"RR_skip={o2[1]:.3f} p={p2[1]:.4f}")
        return
    raise SystemExit("mode?")


if __name__ == "__main__":
    main()
