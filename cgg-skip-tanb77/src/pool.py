"""Pooled frame-break test over a clade: summarize.py's two routes on several genomes at once (README Gate, pooled).
usage: python -I src/pool.py <pools.tsv>                       pools.tsv rows: pool <tab> <outdir>/<name>
       python -I src/pool.py --set <set.tsv> <outdir> <rank> [seed,seed,...]
           pools = the set's scanned genomes grouped by a GTDB rank column; each seed adds a shuffled copy of the
           same pools (same sizes, same genomes, pool names shuf<seed>_<i>), the empirical null
-> TSV on stdout, one row per pool with >= 2 members used.
Small and divergent genomes have few aligned loci, so a code in a few percent of genes gives each genome 2-5 motif
proteins, below S_GENES; a clade sharing the code pools them.
A member is used if it has a background table and fewer than MAX_BREAK clean breaks per 1,000 loci (ONT-type indel
load, README 12:56). In --set mode the cap is applied before grouping and shuffling.
Route S: each subject protein with a prog break anywhere in the pool counts once, under the motif read by the most
members (ties: first member, then table order); the commonest motif's proteins and their share of all prog proteins
(same thresholds as one genome).
Route K: one clean window per subject protein across the pool (first member listed wins). A window from genome g holds
4-mer k with g's own background fraction q_g(k), so the expected count is sum n_g q_g(k) and the tail is
Poisson-binomial (exact, per-genome binomials convolved): a genome mix cannot fake enrichment through composition.
kmer_eff = (observed - expected) / windows. Every 4-mer is tested; the pool passes if any 4-mer passes all three
thresholds (reported: the best passing one, else the best lp). kmer_junction: share of the reported 4-mer's windows
that are -1 gene junctions (two adjacent genes aligned to one fused subject; the artifact class the pool audit,
scratch/review_pool, found at small effect in large pools).
"""
import csv
import gzip
import math
import random
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize import K, K_EFF, K_LP, K_N, KMERS, MIN_ALN, S_GENES, S_SHARE, flag  # noqa: E402

MAX_BREAK = 100
LG = [0.0]


def lgamma_table(n):
    while len(LG) <= n + 1:
        LG.append(math.lgamma(len(LG) + 1))   # LG[i] = log(i!)


def binom_pmf(n, q):
    lgamma_table(n)
    k = np.arange(n + 1)
    lg = np.array(LG[:n + 1])
    return np.exp(lg[n] - lg - lg[::-1] + k * math.log(q) + (n - k) * math.log(1 - q))


def tail_lp(obs, parts):
    """-log10 P(X >= obs), X = sum of independent Bin(n, q) over parts [(n, q)]"""
    if obs <= 0:
        return 0.0
    pmf = np.array([1.0])
    for n, q in parts:
        if n:
            pmf = np.convolve(pmf, binom_pmf(n, min(max(q, 1e-9), 1 - 1e-9)))
    t = float(pmf[obs:].sum()) if obs < len(pmf) else 0.0
    return 300.0 if t <= 1e-300 else -math.log10(t)


def junction(r):
    """-1 break whose abandoned-frame stop (within 3 codons) overlaps a start codon in the new frame: two adjacent
    genes aligned to one fused subject (from scratch/review_pool/atga_junction.py)"""
    if r["type"] != "-1" or r["a_stop"] == "nan" or int(r["a_stop"]) > 3:
        return False
    s = r["stop_motif"].replace("-", "")
    return any(s[i:i + 3] in ("ATG", "GTG", "TTG") for i in (2, 5)) or s[5:9] == "ATGA"


def load(path):
    """one genome's tables reduced to what pool() needs"""
    path = Path(path)
    outdir, name = path.parent, path.name
    loci = sum(1 for h in csv.DictReader(gzip.open(outdir / f"{name}.hits.tsv.gz", "rt"), delimiter="\t")
               if int(h["alen"]) >= MIN_ALN)
    fs = list(csv.DictReader(gzip.open(outdir / f"{name}.fs.tsv.gz", "rt"), delimiter="\t"))
    flag(fs)
    bgf = outdir / f"{name}.bg.tsv.gz"
    bg = [l.strip() for l in gzip.open(bgf, "rt")][1:] if bgf.exists() else []
    m = dict(path=str(path), name=name, loci=loci, clean=sum(1 for r in fs if r["clean"]),
             prog=[(r["type"], r["stop_motif"][3:], r["sseqid"]) for r in fs if r["prog"]], bg=bool(bg), win=[])
    if bg:
        L = len(bg[0])
        c = Counter(k for w in bg for k in {w[i:i + K] for i in range(len(w) - K + 1)})
        m["bgfrac"] = {k: (c[k] + 1) / (len(bg) + 2) for k in KMERS}
        m["win"] = [(r["sseqid"], r["win_seq"], junction(r)) for r in fs
                    if r["clean"] and len(r.get("win_seq", "")) == L]
    return m


def usable(m):
    return m["bg"] and 1000 * m["clean"] / max(m["loci"], 1) < MAX_BREAK


def pool(pname, members):
    used = [m for m in members if usable(m)]
    # route S: one motif per subject protein, the one read by the most members
    votes = {}
    for i, m in enumerate(used):
        for t, s, sid in dict.fromkeys(m["prog"]):
            votes.setdefault(sid, {}).setdefault((t, s), set()).add(i)
    mot = Counter(max(v.items(), key=lambda x: len(x[1]))[0] for v in votes.values())   # ties: first inserted
    top, topn = max(mot.items(), key=lambda x: (x[1], x[0]), default=(("", ""), 0))
    nprog = len(votes)
    # route K: one window per subject protein across the pool, each tagged with its member's background
    win, seen = [], set()
    for i, m in enumerate(used):
        for sid, w, j in m["win"]:
            if sid not in seen:
                seen.add(sid)
                win.append((i, w, j))
    per = Counter(i for i, _, _ in win)
    has = defaultdict(int)
    for _, w, _ in win:
        for k in {w[a:a + K] for a in range(len(w) - K + 1)}:
            has[k] += 1
    best, best_pass = ("", 0, 0.0, 0.0, 0.0), None
    for k in KMERS:
        exp = sum(n * used[i]["bgfrac"][k] for i, n in per.items())
        if has[k] <= exp:
            continue
        lp = tail_lp(has[k], [(n, used[i]["bgfrac"][k]) for i, n in per.items()])
        eff = (has[k] - exp) / len(win)
        row = (k, has[k], exp, lp, eff)
        if lp > best[3]:
            best = row
        if lp >= K_LP and has[k] >= K_N and eff >= K_EFF and (best_pass is None or lp > best_pass[3]):
            best_pass = row
    km, kn, kexp, klp, eff = best_pass or best
    kj = [j for _, w, j in win if km and any(w[a:a + K] == km for a in range(len(w) - K + 1))]
    route_s = topn >= S_GENES and topn / max(nprog, 1) >= S_SHARE
    route_k = best_pass is not None
    return dict(pool=pname, members=len(members), used=len(used), loci=sum(m["loci"] for m in used),
                clean=sum(m["clean"] for m in used), prog=sum(len(m["prog"]) for m in used), prog_genes=nprog,
                top_motif=f"{top[0]}:{top[1]}", top_genes=topn, top_share=f"{topn / max(nprog, 1):.2f}",
                windows=len(win), kmer=km, kmer_n=kn, kmer_exp=f"{kexp:.1f}", kmer_lp=f"{klp:.1f}",
                kmer_eff=f"{eff:.2f}", kmer_junction=f"{sum(kj) / max(len(kj), 1):.2f}",
                route_s=int(route_s), route_k=int(route_k), flagged=int(route_s or route_k),
                names=",".join(m["name"] for m in used))


def main():
    groups = defaultdict(list)
    if sys.argv[1] == "--set":
        setf, outdir, rank = sys.argv[2], Path(sys.argv[3]), sys.argv[4]
        seeds = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else []
        rows = [r for r in csv.DictReader(open(setf), delimiter="\t") if (outdir / f"{r['gca']}.fs.tsv.gz").exists()]
        n = Counter(r[rank] for r in rows)
        rows = [r for r in rows if n[r[rank]] >= 2]          # only genomes that could be pooled are loaded
        with ProcessPoolExecutor(4) as ex:
            ms = list(ex.map(load, [outdir / r["gca"] for r in rows], chunksize=8))
        for r, m in zip(rows, ms):
            if usable(m):                                     # the cap comes before grouping and shuffling
                groups[r[rank]].append(m)
        groups = {p: v for p, v in groups.items() if len(v) >= 2}
        real = list(groups.items())
        for s in seeds:
            allm = [m for _, v in real for m in v]
            random.Random(s).shuffle(allm)
            it = iter(allm)
            groups.update({f"shuf{s}_{i}": [next(it) for _ in v] for i, (_, v) in enumerate(real)})
    else:
        for line in open(sys.argv[1]):
            p, path = line.rstrip("\n").split("\t")
            groups[p].append(load(path))
    w = None
    for p, ms in groups.items():
        if len(ms) < 2:
            continue
        row = pool(p, ms)
        if row["used"] < 2:
            continue
        if w is None:
            w = csv.DictWriter(sys.stdout, list(row), delimiter="\t", lineterminator="\n")
            w.writeheader()
        w.writerow(row)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
