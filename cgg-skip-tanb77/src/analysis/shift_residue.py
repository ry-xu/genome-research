"""Which homolog residue sits at a +1 break? For each clean +1 break of a genome, translate the old frame up to the
break codon and the new frame from 4 nt after its start, join them with one X, align that to the subject around spos
(local, BLOSUM62), and report the subject residue opposite X, plus the 6 nt from the break codon on.
usage: python -I src/analysis/shift_residue.py <outdir> <name> <fasta> [codon]  -> TSV on stdout, residue tally on stderr
Arg opposite X at CGG-A breaks: CGG(A) read as Arg; Gly: the C is skipped and GGA read (or a P-site slip onto GGA)."""
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from candidate import read_fasta  # noqa: E402
from summarize import flag  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
FLANK = 20   # codons each side


def sprot(ids):
    out, cur = {}, None
    for line in open(ROOT / "data/sprot_p50.faa"):
        if line[0] == ">":
            cur = line[1:].split()[0]
            cur = cur if cur in ids else None
            if cur:
                out[cur] = []
        elif cur:
            out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def main():
    outdir, name, fasta = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    want = sys.argv[4] if len(sys.argv) > 4 else None
    fs = list(csv.DictReader(gzip.open(outdir / f"{name}.fs.tsv.gz", "rt"), delimiter="\t"))
    flag(fs)
    rows = [r for r in fs if r["clean"] and r["type"] == "+1" and (want is None or r["old_codon"] == want)]
    seqs = read_fasta(fasta)
    subj = sprot({r["sseqid"] for r in rows})
    al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-11, extend_gap_score=-1)
    tally = Counter()
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["gpos", "strand", "sseqid", "spos", "old_codon", "next6", "up", "down", "sres", "score", "k_lo", "k_hi"])
    for r in rows:
        g = seqs[r["contig"]]
        ctx = r["ctx"].upper()
        s = str(Seq(g).reverse_complement()) if r["strand"] == "-" else g
        i = s.find(ctx)
        if i < 0 or s.find(ctx, i + 1) >= 0:
            continue
        site = i + sum(c.islower() for c in r["ctx"])   # ctx: lower-case upstream, the break codon first upper-case
        up = str(Seq(s[site - 3 * FLANK:site]).translate())
        down = str(Seq(s[site + 4:site + 4 + 3 * FLANK]).translate())
        q = up + "X" + down
        sp = int(r["spos"])
        sq = subj.get(r["sseqid"], "")
        lo = max(sp - 1 - 2 * FLANK, 0)
        win = sq[lo:sp - 1 + 2 * FLANK]
        if not win or "*" in up[-8:] or "*" in down[:8]:
            continue
        a = al.align(q, win)[0]
        res = "."
        for (qs, qe), (ts, te) in zip(*a.aligned):
            if qs <= FLANK < qe:
                res = win[ts + FLANK - qs]
        tally[res] += 1
        w.writerow([r["gpos"], r["strand"], r["sseqid"], sp, r["old_codon"], s[site:site + 6], up[-6:], down[:6], res,
                    f"{a.score:.0f}", r["k_lo"], r["k_hi"]])
    print(" ".join(f"{k}:{v}" for k, v in tally.most_common()), file=sys.stderr)


if __name__ == "__main__":
    main()
