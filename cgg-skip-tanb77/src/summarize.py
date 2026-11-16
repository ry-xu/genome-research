"""Per-genome frame-break summary from fs_scan tables.
usage: python -I src/summarize.py <outdir> [name ...]   (default: every <name>.fs.tsv.gz in outdir) -> TSV on stdout

A break is CLEAN when it is not one of a PAIR (opposite breaks within PAIR_NT nt), its alignment has >= MIN_ALN aa,
both 15-column flanks are fully aligned with identity >= FLANK_ID (a working gene on both sides, not a decaying
pseudogene), the longest homopolymer touching its tie window is <= MAX_HP, and no run of > MAX_HP lies within NEAR nt
either side of the window (the nanopore/454 indel class; another agent's audit: run-edge breaks pass a touching-run test).
It is PROGRAMMED-LIKE (prog) when it passes the same pair, length and flank tests and the abandoned frame reads a stop
within A_STOP codons of the tie window (prfB, dnaX, IS elements, Euplotes all stop the old frame at or just after the
shift). prog has no homopolymer test: prfB's site sits next to its GGGGG Shine-Dalgarno run.
Columns: loci (hits >= MIN_ALN aa), breaks, clean, prog, prog_per_kloci, frac_prog (prog / clean), top_motif /
top_genes / top_share (route S: the commonest stop-anchored motif among prog, type + P-site codon + stop, in distinct
subject proteins, and its share of prog proteins), loci_with_prog,
and the motif-agnostic route: kmer / kmer_n / kmer_bg / kmer_lp = the 4-mer (gene strand, homopolymers excluded) most
over-represented in clean break windows (WINSEQ nt from 3 codons before the tie window) against random windows of the
same length in the same genome's aligned loci: count of clean breaks whose window has it, the background fraction of
windows with it, and -log10 of the binomial upper tail; kmer_eff = window fraction minus background fraction.
route_s / route_k / flagged apply the frozen gate (README).
"""
import csv
import gzip
import itertools
import math
import sys
from collections import Counter
from pathlib import Path

K = 4
KMERS = ["".join(x) for x in itertools.product("ACGT", repeat=K) if len(set(x)) > 1]


def binom_lsf(k, n, p):
    """-log10 P(X >= k), X ~ Bin(n, p), summed in log space."""
    if k <= 0:
        return 0.0
    p = min(max(p, 1e-9), 1 - 1e-9)
    terms = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * math.log(p)
             + (n - i) * math.log(1 - p) for i in range(k, n + 1)]
    m = max(terms)
    return -(m + math.log(sum(math.exp(t - m) for t in terms))) / math.log(10)


def kmer_route(wins, bg):
    """best (kmer, count, background fraction, -log10 p) over KMERS"""
    if not wins or not bg:
        return "", 0, 0.0, 0.0
    obs = Counter(k for w in wins for k in {w[i:i + K] for i in range(len(w) - K + 1)})
    bgc = Counter(k for w in bg for k in {w[i:i + K] for i in range(len(w) - K + 1)})
    best = ("", 0, 0.0, 0.0)
    for k in KMERS:
        q = (bgc[k] + 1) / (len(bg) + 2)
        lp = binom_lsf(obs[k], len(wins), q)
        if lp > best[3]:
            best = (k, obs[k], q, lp)
    return best

PAIR_NT = 150     # opposite breaks this close in one alignment = a detour around an in-frame stop or a double indel
                  # (was 30 until the 12:42 audit amendment: detour spacing is flat out to ~150 nt)
FLANK_ID = 0.20
MAX_HP = 4
A_STOP = 3
MIN_ALN = 50
NEAR = 6
# the gate (README, frozen 2026-11-11 12:25 EST)
S_GENES, S_SHARE = 6, 0.25
K_LP, K_N, K_EFF = 6.0, 8, 0.20


def run_over(x, a, b):
    """longest single-base run of x overlapping x[a:b], measured over all of x (not cut at a or b)"""
    best, i = 0, 0
    for _, g in itertools.groupby(x.upper()):
        n = len(list(g))
        if i < b and i + n > a:
            best = max(best, n)
        i += n
    return best


def near_run(r):
    """longest run reaching from NEAR nt before the tie window to NEAR nt past the site codon (win_seq starts 9 nt
    before the window; site_ctx is 12 nt either side of the site)"""
    n = 3 * (int(r["k_hi"]) - int(r["k_lo"])) + 3 + 2 * NEAR
    return max(run_over(r["win_seq"], 9 - NEAR, 9 - NEAR + n), run_over(r["site_ctx"], 12 - NEAR, 15 + NEAR))


def flag(fs):
    """set r["pair"], r["clean"], r["prog"] on every break row; return the number of pairs"""
    by = {}
    for r in fs:
        by.setdefault((r["contig"], r["strand"], r["sseqid"]), []).append(r)
    npair = 0
    for v in by.values():
        v.sort(key=lambda r: int(r["gpos"]))
        for i, a in enumerate(v):
            for b in v[i + 1:]:
                if int(b["gpos"]) - int(a["gpos"]) > PAIR_NT:
                    break
                if a["type"] != b["type"]:
                    a["pair"] = b["pair"] = "1"
                    npair += 1
    for r in fs:
        good = ("pair" not in r and int(r["aln_aa"]) >= MIN_ALN and r["left_n"] == r["right_n"] == "15"
                and float(r["left_id"]) >= FLANK_ID and float(r["right_id"]) >= FLANK_ID)
        r["good"] = good
        r["clean"] = good and int(r["win_hp"]) <= MAX_HP and near_run(r) <= MAX_HP
        # no homopolymer test here: the frame-A stop is the anchor, and real sites sit next to runs (prfB's GGGGG
        # Shine-Dalgarno, A-AAA-AAG slippery sites); per-gene motif counting keeps run indels from concentrating
        r["prog"] = good and r["a_stop"] != "nan" and int(r["a_stop"]) <= A_STOP
    return npair


def summarize(outdir, name):
    hits = list(csv.DictReader(gzip.open(outdir / f"{name}.hits.tsv.gz", "rt"), delimiter="\t"))
    fs = list(csv.DictReader(gzip.open(outdir / f"{name}.fs.tsv.gz", "rt"), delimiter="\t"))
    loci = sum(1 for h in hits if int(h["alen"]) >= MIN_ALN)
    npair = flag(fs)
    clean = [r for r in fs if r["clean"]]
    prog = [r for r in fs if r["prog"]]
    # motif support counted in distinct subject proteins, so IS-element copies (one transposase) count once
    mot = Counter(m for m in {(r["type"], r["stop_motif"][3:], r["sseqid"]) for r in prog})
    mot = Counter((t, m) for t, m, _ in mot)
    top, topn = max(((m, n) for m, n in mot.items()), key=lambda x: (x[1], x[0]), default=(("", ""), 0))
    nprog_genes = len({r["sseqid"] for r in prog})
    lp = len({(r["contig"], r["strand"], r["sseqid"]) for r in prog})
    bgf = outdir / f"{name}.bg.tsv.gz"
    bg = [l.strip() for l in gzip.open(bgf, "rt")][1:] if bgf.exists() else []
    one = {}
    for r in clean:                     # one window per subject protein, for the same reason
        one.setdefault(r["sseqid"], r)
    km, kn, kq, klp = kmer_route([r["win_seq"] for r in one.values() if len(r.get("win_seq", "")) == len(bg[0])] if bg else [],
                                 bg)
    route_s = topn >= S_GENES and topn / max(nprog_genes, 1) >= S_SHARE
    eff = kn / max(len(one), 1) - kq
    route_k = klp >= K_LP and kn >= K_N and eff >= K_EFF
    return dict(name=name, loci=loci, breaks=len(fs), pairs=npair, clean=len(clean), clean_genes=len(one),
                clean_p1=sum(r["type"] == "+1" for r in clean), clean_m1=sum(r["type"] == "-1" for r in clean),
                prog=len(prog), prog_per_kloci=f"{1000 * len(prog) / max(loci, 1):.2f}",
                frac_prog=f"{len(prog) / max(len(clean), 1):.2f}", prog_genes=nprog_genes,
                top_motif=f"{top[0]}:{top[1]}", top_genes=topn, top_share=f"{topn / max(nprog_genes, 1):.2f}",
                loci_with_prog=lp, kmer=km, kmer_n=kn, kmer_bg=f"{kq:.3f}", kmer_lp=f"{klp:.1f}", kmer_eff=f"{eff:.2f}",
                route_s=int(route_s), route_k=int(route_k), flagged=int(route_s or route_k))


def main():
    outdir = Path(sys.argv[1])
    names = sys.argv[2:] or sorted(p.name[:-10] for p in outdir.glob("*.fs.tsv.gz"))
    w = None
    for n in names:
        r = summarize(outdir, n)
        if w is None:
            w = csv.DictWriter(sys.stdout, list(r), delimiter="\t", lineterminator="\n")
            w.writeheader()
        w.writerow(r)


if __name__ == "__main__":
    main()
