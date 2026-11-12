"""Is the sequence after a TTA "stop" translated? Codon-position substitution pattern between close clade orthologs.

usage: runon_coding.py <panel.tsv> <genes_dir> <genomes_dir> <ortho_dir> <out_prefix> [pairs, e.g. C_a:C_b,...]

For RBH ortholog pairs (g1, g2) between two clade genomes, each gene's class comes from runon_events.tsv with
R = sister (codon at e*+1 relative to its best sister ortholog that ends at a canonical stop):
  TT    both genes have TTA at e*+1 and run on (test set)
  SS    both genes end with a canonical stop at e*+1 ("same" class; non-coding control after the stop)
  RR    both genes run on with a sense codon other than TTA at e*+1 (C-terminal extensions vs the sister; added to
        see whether short run-on tails that may be coding show codon structure at all)
Regions, in gene orientation, anchored at the end codon (TTA or canonical stop) at e*+1:
  up    the 150 nt before it (coding, positive control)
  down  the nt after it: for TT, up to the shorter of the two run-on ORFs (the part that would be translated if TTA
        were Leu), capped at 150 nt; for SS, 150 nt of contig after the stop (non-coding control) and, separately,
        the same length cap as TT
Alignment: Biopython PairwiseAligner, global with free end gaps on the side away from the anchor
(match 2, mismatch -3, open -6, extend -2). Pairs whose up alignment has < 60% identity over the last 30 aligned
columns are dropped (anchor not homologous).
Per region class, over aligned gap-free columns where both bases sit at the same codon position (relative to the
anchor): mismatch fraction at codon positions 1, 2, 3; R3 = p3 / mean(p1, p2); indel events and how many have
length divisible by 3. Bootstrap over gene pairs (1,000, seed 1) for a 95% interval of R3.
"""
import gzip
import warnings
import random
import sys
from collections import defaultdict

from Bio import Align

warnings.filterwarnings("ignore")

COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"TAA", "TAG", "TGA"}


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(fn):
    op = gzip.open if fn.endswith(".gz") else open
    name, buf = None, []
    with op(fn, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf).upper()


def downstream(contig, b, e, strand, off, n):
    """n nt after gene-coordinate offset `off` (0-based, in gene orientation)."""
    if strand == 1:
        s = b - 1 + off
        return contig[s:s + n]
    end0 = e - 1 - off  # contig index of gene-coordinate off
    lo = max(0, end0 - n + 1)
    return rc(contig[lo:end0 + 1])


def aligner(free_left):
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.match_score, a.mismatch_score = 2, -3
    a.open_gap_score, a.extend_gap_score = -6, -2
    if free_left:
        a.target_left_open_gap_score = a.target_left_extend_gap_score = 0
        a.query_left_open_gap_score = a.query_left_extend_gap_score = 0
    else:
        a.target_right_open_gap_score = a.target_right_extend_gap_score = 0
        a.query_right_open_gap_score = a.query_right_extend_gap_score = 0
    return a


def columns(aln):
    """-> list of (i1 or None, i2 or None) aligned columns."""
    cols = []
    c1, c2 = aln.coordinates
    for k in range(len(c1) - 1):
        a0, a1, b0, b1 = c1[k], c1[k + 1], c2[k], c2[k + 1]
        if a1 - a0 and b1 - b0:
            cols += [(a0 + j, b0 + j) for j in range(a1 - a0)]
        elif a1 - a0:
            cols += [(a0 + j, None) for j in range(a1 - a0)]
        else:
            cols += [(None, b0 + j) for j in range(b1 - b0)]
    return cols


def tally(s1, s2, cols, frame_of, trim_ends=True):
    """Per-codon-position mismatch counts and internal indels. frame_of(i) -> codon position 0/1/2."""
    mm, tot = [0, 0, 0], [0, 0, 0]
    # internal region: between first and last match column
    both = [k for k, (i, j) in enumerate(cols) if i is not None and j is not None]
    if not both:
        return mm, tot, []
    lo, hi = both[0], both[-1]
    indels = []
    run = 0
    for k in range(lo, hi + 1):
        i, j = cols[k]
        if i is not None and j is not None:
            if run:
                indels.append(run)
                run = 0
            p1, p2 = frame_of(i, 1), frame_of(j, 2)
            if p1 != p2:
                continue
            tot[p1] += 1
            if s1[i] != s2[j]:
                mm[p1] += 1
        else:
            run += 1
    return mm, tot, indels


def main(panel_fn, gdir, gendir, odir, outp, pairs_arg=None):
    panel = [l.rstrip("\n").split("\t") for l in open(panel_fn)]
    clade = [r[0] for r in panel if r[2] == "clade"]
    acc = {r[0]: r[1] for r in panel}
    want = None
    if pairs_arg:
        want = {tuple(sorted(p.split(":"))) for p in pairs_arg.split(",")}
    genes = {}
    for lab in clade:
        with open(f"{gdir}/{lab}.genes.tsv") as f:
            next(f)
            for line in f:
                gid, ctg, b, e, st, pb, pe, nt = line.rstrip("\n").split("\t")
                genes[f"{lab}|{gid}"] = (ctg, int(b), int(e), int(st), nt)
    cls = {}
    with open(f"{odir}/runon_events.tsv") as f:
        h = next(f).rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(h, line.rstrip("\n").split("\t")))
            if r["R"] != "sister" or r["B_role"] != "clade":
                continue
            g = f"{r['B']}|{r['gene']}"
            es, ql = int(r["e_star"]), int(r["ql"])
            if r["class"] == "runon" and r["w_p1"] == "TTA":
                cls[g] = ("TT", es, ql)
            elif r["class"] == "same":
                cls[g] = ("SS", es, ql)
            elif r["class"] == "runon":
                cls[g] = ("RR", es, ql)
    contigs = {}
    for lab in clade:
        contigs[lab] = dict(read_fasta(f"{gendir}/{acc[lab]}.fna.gz"))
    upal, downal = aligner(True), aligner(False)
    res = defaultdict(list)  # region -> list of per-pair (mm, tot, indels, pair)
    seen = set()
    nchk = nbad = 0
    with open(f"{odir}/rbh.tsv") as f:
        next(f)
        for line in f:
            q, s = line.split("\t")[:2]
            a, b = q.split("|")[0], s.split("|")[0]
            if a not in clade or b not in clade or a > b:
                continue
            if want is not None and (a, b) not in want:
                continue
            if q not in cls or s not in cls or (q, s) in seen:
                continue
            seen.add((q, s))
            c1, c2 = cls[q], cls[s]
            if c1[0] != c2[0]:
                continue
            kind = c1[0]
            regs = {}
            seqs = []
            for g, (k, es, ql) in ((q, c1), (s, c2)):
                lab = g.split("|")[0]
                ctg, gb, ge, st, nt = genes[g]
                cod = nt[3 * es:3 * es + 3]
                if kind == "TT":
                    assert cod == "TTA", (g, cod)
                elif kind == "SS":
                    assert cod in STOPS, (g, cod)
                else:
                    assert cod not in STOPS and cod != "TTA", (g, cod)
                # sanity: contig extraction agrees with gene nt
                nchk += 1
                if downstream(contigs[lab][ctg], gb, ge, st, 0, len(nt)) != nt:
                    nbad += 1
                up = nt[max(0, 3 * es - 150):3 * es]
                dn_contig = downstream(contigs[lab][ctg], gb, ge, st, 3 * es + 3, 150)
                runon_nt = len(nt) - (3 * es + 3) - 3  # sense nt of the run-on ORF after the TTA, stop excluded
                seqs.append((up, dn_contig, runon_nt))
            (u1, d1, r1), (u2, d2, r2) = seqs
            if len(u1) < 90 or len(u2) < 90:
                continue
            aln = upal.align(u1, u2)[0]
            cols = columns(aln)
            last = [(i, j) for i, j in cols if i is not None and j is not None][-30:]
            ident = sum(u1[i] == u2[j] for i, j in last) / max(1, len(last))
            # the anchor must be aligned: the last column pairs the last base of both
            if ident < 0.6 or cols[-1] != (len(u1) - 1, len(u2) - 1):
                continue
            mm, tot, ind = tally(u1, u2, cols, lambda i, w: (i - (len(u1) if w == 1 else len(u2))) % 3)
            res[f"{kind}_up"].append((mm, tot, ind, (q, s)))
            if kind in ("TT", "RR"):
                L = min(r1, r2, 150)
                variants = [(f"{kind}_down_runon", d1[:L], d2[:L])]
            else:
                variants = [("SS_down150", d1, d2)]
            for name, x1, x2 in variants:
                if len(x1) < 15 or len(x2) < 15:
                    continue
                aln = downal.align(x1, x2)[0]
                cols = columns(aln)
                mm, tot, ind = tally(x1, x2, cols, lambda i, w: i % 3)
                res[name].append((mm, tot, ind, (q, s)))
    print(f"contig extraction check: {nbad} mismatches of {nchk}", file=sys.stderr)
    rng = random.Random(1)
    with open(f"{outp}.tsv", "w") as out:
        out.write("region\tpairs\tcols_p1\tcols_p2\tcols_p3\tp1\tp2\tp3\tR3\tR3_lo95\tR3_hi95\tindels\tindels_mod3_0\n")
        for name in sorted(res):
            rows = res[name]

            def agg(rr):
                mm = [sum(r[0][k] for r in rr) for k in range(3)]
                tt = [sum(r[1][k] for r in rr) for k in range(3)]
                p = [mm[k] / tt[k] if tt[k] else float("nan") for k in range(3)]
                den = (p[0] + p[1]) / 2
                return p, tt, (p[2] / den if den else float("inf"))

            p, tt, r3 = agg(rows)
            boots = []
            for _ in range(1000):
                bs = [rows[rng.randrange(len(rows))] for _ in rows]
                boots.append(agg(bs)[2])
            boots.sort()
            ind = [x for r in rows for x in r[2]]
            out.write(f"{name}\t{len(rows)}\t{tt[0]}\t{tt[1]}\t{tt[2]}\t{p[0]:.4f}\t{p[1]:.4f}\t{p[2]:.4f}\t{r3:.3f}\t"
                      f"{boots[25]:.3f}\t{boots[974]:.3f}\t{len(ind)}\t{sum(1 for x in ind if x % 3 == 0)}\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
