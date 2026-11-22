"""Stage 2 for one called (unit, stop codon, meaning): is the stop read as a sense codon on the organism's own genes?

usage: stage2.py UNITS.tsv UNIT CODON MEANING THREADS > report.txt
Re-aligns the unit's bodies exactly as stopsense.py does and keeps per-body data. Each test is also run on a
like-for-like control: the sense codon for the meaning (or, if none is close, any sense codon) whose inner count is
nearest the stop's, so codon-specific effects (AT-rich codons sit in AT-rich genes) are compared, not assumed away.
  dispersion   stops per body vs Poisson with the body's inner length: observed / expected bodies with >= 1, and the
               Pearson dispersion index. A minority source (contaminant, prey) puts stops in few bodies.
  composition  GC3 of bodies with the codon vs without (median, Mann-Whitney p), and the codon-usage distance between
               the two pools vs 1,000 random splits of the same sizes (percentile).
  pident       best-hit identity, bodies with vs without.
  conserved    share of the stop's columns aligned to the meaning where the 10 flanking columns are >= 60% identical,
               vs the same share for the meaning's own codons (ident at conserved columns).
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # runs under python -I
from stopsense import AA, CI, CODE1, CODONS, EDGE, REF, TOK, bodies  # noqa: E402

STD = set("ACDEFGHIKLMNPQRSTVWY")
FLANK, CONS = 5, 0.6


def align(prefix, threads):
    body = {k: (contig, s, b) for k, contig, s, b in bodies(prefix)}
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        fa, hf = os.path.join(d, "q.fa"), os.path.join(d, "h.tsv")
        with open(fa, "w") as f:
            for k, (_, _, b) in body.items():
                if len(b) >= 3 * (2 * EDGE + 1):
                    f.write(f">u|{k}\n{b}\n")
        subprocess.run(["diamond", "blastx", "-q", fa, "-d", REF, "-o", hf, "--fast", "-k", "1", "--max-hsps", "1",
                        "-e", "1e-5", "--threads", str(threads), "--quiet", "--outfmt", "6", "qseqid", "sseqid",
                        "qstart", "qend", "sstart", "send", "slen", "pident", "qframe", "btop"], check=True)
        hits = [l.rstrip("\n").split("\t") for l in open(hf)]
    out = []
    for q, s, qs, qe, ss, se, sl, pid, fr, btop in hits:
        k = int(q.split("|", 1)[1])
        contig, _, b = body[k]
        qs = int(qs)
        if int(fr) != 1 or (qs - 1) % 3:
            continue
        cols, pos = [], qs - 1
        for n, p in TOK.findall(btop):
            if n:
                for _ in range(int(n)):
                    c = b[pos:pos + 3]
                    cols.append((c, CODE1.get(c, "X"), True))
                    pos += 3
            elif p[0] == "-":
                cols.append(None)
            elif p[1] == "-":
                cols.append(None)
                pos += 3
            else:
                cols.append((b[pos:pos + 3], p[1], False))
                pos += 3
        inner = []
        for j in range(EDGE, len(cols) - EDGE):
            x = cols[j]
            if x is None or x[0] not in CI or x[1] not in STD:
                continue
            fl = [cols[i] for i in range(j - FLANK, j + FLANK + 1) if i != j]
            loc = sum(1 for y in fl if y is not None and y[2]) / len(fl)
            inner.append((x[0], x[1], loc))
        cod = [b[i:i + 3] for i in range(0, len(b) - 2, 3)]
        gc3 = [c[2] in "GC" for c in cod if c in CI and AA[CI[c]] != "*"]
        use = np.zeros(64)
        for c in cod:
            if c in CI:
                use[CI[c]] += 1
        out.append(dict(contig=contig, sseqid=s, pid=float(pid), inner=inner, gc3=np.mean(gc3) if gc3 else np.nan,
                        use=use))
    return out


def tests(B, codon, meaning, rng):
    c = np.array([sum(1 for x in b["inner"] if x[0] == codon) for b in B], float)
    L = np.array([len(b["inner"]) for b in B], float)
    S, C = c.sum(), L.sum()
    e = S * L / C
    has = c > 0
    r = dict(n=int(S), bodies=len(B), obs_pos=int(has.sum()), exp_pos=float((1 - np.exp(-e)).sum()),
             disp=float(((c - e) ** 2 / np.where(e > 0, e, 1)).sum() / max(len(B) - 1, 1)))
    g = np.array([b["gc3"] for b in B])
    a, z = g[has & ~np.isnan(g)], g[~has & ~np.isnan(g)]
    r.update(gc3_with=float(np.median(a)) if len(a) else np.nan, gc3_without=float(np.median(z)) if len(z) else np.nan,
             gc3_p=float(mannwhitneyu(a, z).pvalue) if len(a) and len(z) else np.nan)
    U = np.array([b["use"] for b in B])
    sense = [i for i, cc in enumerate(CODONS) if AA[i] != "*"]

    def dist(mask):
        p, q = U[mask][:, sense].sum(0), U[~mask][:, sense].sum(0)
        return float(np.abs(p / p.sum() - q / q.sum()).sum() / 2)

    d0 = dist(has)
    null = [dist(rng.permutation(has)) for _ in range(1000)]
    r.update(use_dist=d0, use_pct=float(np.mean(np.array(null) < d0)))
    pid = np.array([b["pid"] for b in B])
    r.update(pid_with=float(np.median(pid[has])), pid_without=float(np.median(pid[~has])))
    cons = [x for b in B for x in b["inner"] if x[0] == codon and x[2] >= CONS]
    r.update(cons_n=len(cons), cons_share=float(np.mean([x[1] == meaning for x in cons])) if cons else np.nan)
    return r


def main(units_tsv, unit, codon, meaning, threads):
    prefix = next(l.split("\t")[2] for l in open(units_tsv) if l.split("\t")[0] == unit)
    B = align(prefix, threads)
    rng = np.random.default_rng(1)
    counts = {cc: sum(1 for b in B for x in b["inner"] if x[0] == cc) for cc in CODONS}
    n0 = counts[codon]
    own = [cc for cc in CODONS if AA[CI[cc]] == meaning]
    pool = own if any(0.25 <= counts[cc] / max(n0, 1) <= 4 for cc in own) else [cc for cc in CODONS if AA[CI[cc]] != "*"]
    ctrl = min(pool, key=lambda cc: abs(np.log((counts[cc] + 1) / (n0 + 1))))
    own_cons = [x for b in B for x in b["inner"] if AA[CI[x[0]]] == meaning and x[2] >= CONS]
    print(f"unit {unit}  codon {codon}  meaning {meaning}  bodies aligned {len(B)}  control codon {ctrl} "
          f"({AA[CI[ctrl]]})")
    print(f"meaning's own codons at conserved columns: n={len(own_cons)} "
          f"ident={np.mean([x[1] == meaning for x in own_cons]):.3f}")
    keys = ["n", "bodies", "obs_pos", "exp_pos", "disp", "gc3_with", "gc3_without", "gc3_p", "use_dist", "use_pct",
            "pid_with", "pid_without", "cons_n", "cons_share"]
    print("\t".join(["codon"] + keys))
    for cc, m in ((codon, meaning), (ctrl, AA[CI[ctrl]])):
        r = tests(B, cc, m, rng)
        print("\t".join([cc] + [f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in keys]))
    print("contigs with the codon at inner columns (contig, sseqid, residues):")
    for b in B:
        res = [x[1] for x in b["inner"] if x[0] == codon]
        if res:
            print(f"  {b['contig']}\t{b['sseqid']}\t{''.join(res)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
