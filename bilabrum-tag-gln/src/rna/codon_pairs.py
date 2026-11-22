"""Genome codon vs RNA-read codon at every aligned codon of the 76-protein tblastn panel (between-strain changes).

usage: codon_pairs.py GENOME HITS.tsv OUTPREFIX ALN.paf [ALN.paf ...]
HITS is the reviewer's bil_human76.tsv (tblastn outfmt 6 with qseq/sseq). Each genomic codon is taken once, from the
highest-bitscore HSP covering it (as stops_in_hsps.py does), only where the query has a residue. Reads: primary
alignments with mapq >= 10; a read counts at a codon if all three bases align without an indel inside the codon.
A site is called with a read codon when >= 2 reads cover it and >= 2/3 agree; the rest are 'mixed'.
Question: does genomic CAG turn into TAG between strains as often as Leu CTG turns into TTG (a synonymous first-position
C->T), and much more often than CAA->TAA or CGA->TGA (C->T that makes TAA or TGA)? And the reverse T->C changes.
Writes OUTPREFIX.codons.tsv.gz (per site) and OUTPREFIX.pairs.txt.
"""
import bisect
import os
import gzip
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "review"))
from stops_in_hsps import read_fasta, rc, walk  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from site_reads import cs_walk  # noqa: E402

TESTS = [  # (label, genome codon, read codon)
    ("C->T, Gln CAG to TAG", "CAG", "TAG"),
    ("C->T, Leu CTG to TTG (synonymous)", "CTG", "TTG"),
    ("C->T, Leu CTA to TTA (synonymous)", "CTA", "TTA"),
    ("C->T, Gln CAA to TAA (TAA stop)", "CAA", "TAA"),
    ("C->T, Arg CGA to TGA (TGA stop)", "CGA", "TGA"),
    ("C->T, His CAC/CAT to TAC/TAT (H->Y)", ("CAC", "CAT"), ("TAC", "TAT")),
    ("C->T, Arg CGG to TGG (R->W)", "CGG", "TGG"),
    ("C->T, Pro CCN to TCN (P->S)", ("CCA", "CCC", "CCG", "CCT"), ("TCA", "TCC", "TCG", "TCT")),
    ("T->C, TAG to CAG", "TAG", "CAG"),
    ("T->C, Leu TTG to CTG (synonymous)", "TTG", "CTG"),
    ("T->C, Leu TTA to CTA (synonymous)", "TTA", "CTA"),
    ("T->C, Tyr TAC/TAT to CAC/CAT (Y->H)", ("TAC", "TAT"), ("CAC", "CAT")),
    ("T->C, Trp TGG to CGG (W->R)", "TGG", "CGG"),
    ("3rd pos, Gln CAG to CAA", "CAG", "CAA"),
    ("3rd pos, Gln CAA to CAG", "CAA", "CAG"),
]


def main(genome, hits, out, *pafs):
    seqs = read_fasta(genome)
    H = []
    for line in open(hits):
        f = line.rstrip("\n").split("\t")
        H.append(dict(q=f[0], s=f[1].split("|")[1] if "|" in f[1] else f[1], qstart=int(f[6]), sstart=int(f[8]),
                      bits=float(f[11]), sframe=int(f[13]), qseq=f[15], sseq=f[16]))
    best = {}
    for i, h in enumerate(H):
        for col, qp, a, b, cod, p in walk(h, seqs[h["s"]]):
            if a == "-" or len(cod) != 3:
                continue
            lo = p if h["sframe"] > 0 else p - 2
            key = (h["s"], lo, h["sframe"] > 0)
            if key not in best or H[best[key][0]]["bits"] < h["bits"]:
                best[key] = (i, cod, a, qp)
    keys = sorted(best)
    by_c = defaultdict(list)
    for k, key in enumerate(keys):
        by_c[key[0]].append((key[1], k))
    starts = {c: [x[0] for x in v] for c, v in by_c.items()}
    reads = [Counter() for _ in keys]
    for paf in pafs:
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if int(f[11]) < 10 or "tp:A:P" not in f[12:]:
                continue
            c, ts, te = f[5], int(f[7]), int(f[8])
            if c not in by_c:
                continue
            a, b = bisect.bisect_left(starts[c], ts), bisect.bisect_right(starts[c], te - 3)
            if a >= b:
                continue
            cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
            cov, ins = {}, set()
            for t, base, isin in cs_walk(cs, ts):
                if isin:
                    ins.add(t)
                else:
                    cov[t] = base
            for lo, k in by_c[c][a:b]:
                bs = [cov.get(lo + d) for d in range(3)]
                if None in bs or lo in ins or lo + 1 in ins:
                    continue
                cod = "".join(bs)
                reads[k][cod if keys[k][2] else rc(cod)] += 1
    pairs = Counter()
    gcount = Counter()
    with gzip.open(out + ".codons.tsv.gz", "wt") as w:
        w.write("contig\tleft0\tstrand\tquery\tqpos\tqres\tgenome\treads\tcall\tread_codons\n")
        for key, n in zip(keys, reads):
            i, cod, a, qp = best[key]
            tot = sum(n.values())
            call = ""
            if tot >= 2:
                top, m = n.most_common(1)[0]
                call = top if 3 * m >= 2 * tot else "mixed"
                gcount[cod] += 1
                pairs[(cod, call)] += 1
            w.write(f"{key[0]}\t{key[1]}\t{'+' if key[2] else '-'}\t{H[i]['q']}\t{qp}\t{a}\t{cod}\t{tot}\t{call}\t"
                    + ",".join(f"{x}:{v}" for x, v in n.most_common()) + "\n")
    with open(out + ".pairs.txt", "w") as w:
        w.write(f"panel codons {len(keys)}; called (>= 2 reads) {sum(gcount.values())}; "
                f"read codon differs from genome at {sum(v for (g, r), v in pairs.items() if g != r and r != 'mixed')}"
                f"; mixed {sum(v for (g, r), v in pairs.items() if r == 'mixed')}\n")
        w.write("test\tgenome codon sites called\tchanged to target\tper 1,000\n")
        for lab, g, r in TESTS:
            g = (g,) if isinstance(g, str) else g
            r = (r,) if isinstance(r, str) else r
            n = sum(gcount[x] for x in g)
            k = sum(pairs[(x, y)] for x in g for y in r)
            w.write(f"{lab}\t{n}\t{k}\t{1000 * k / n if n else float('nan'):.1f}\n")
        stops = ("TAA", "TAG", "TGA")
        w.write("\nnon-stop genome codon -> stop in reads: " + ", ".join(
            f"{s} {sum(v for (g, r), v in pairs.items() if r == s and g not in stops)}" for s in stops) + "\n")
        w.write("genome stop -> read codon: " + "; ".join(
            f"{s}: " + ", ".join(f"{r} {v}" for (g, r), v in sorted(pairs.items(), key=lambda x: -x[1]) if g == s)
            for s in stops) + "\n")
    print(open(out + ".pairs.txt").read(), end="")


if __name__ == "__main__":
    main(*sys.argv[1:])
