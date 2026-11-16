"""In-frame Arg codon use in frame-continuous hits (n_fs 0, >= 20 codons from each end), per genome.
usage: python -I src/analysis/inframe_cgg.py <set.tsv> <outdir> [gca ...]  -> TSV on stdout (all scanned genomes of the set,
or only the listed ones)
Columns: codons, then per 1,000 codons CGG, CGA, CGT, CGC, AGA, AGG, and the raw CGG count."""
import csv
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from candidate import read_fasta  # noqa: E402

from paths import GENOMES as GDIR  # noqa: E402
COMP = str.maketrans("ACGTN", "TGCAN")
ARG = ["CGG", "CGA", "CGT", "CGC", "AGA", "AGG"]


def main():
    setf, outdir = sys.argv[1], Path(sys.argv[2])
    only = set(sys.argv[3:])
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["gca", "family", "genus", "codons"] + [f"{c}_per_k" for c in ARG] + ["CGG_n"])
    for r in csv.DictReader(open(setf), delimiter="\t"):
        a = r["gca"]
        if (only and a not in only) or not (outdir / f"{a}.hits.tsv.gz").exists():
            continue
        seqs = read_fasta(str(GDIR / a[:8] / f"{a}.fna.gz"))
        st = {"+": seqs, "-": {k: v.translate(COMP)[::-1] for k, v in seqs.items()}}
        n, cnt = 0, dict.fromkeys(ARG, 0)
        for h in csv.DictReader(gzip.open(outdir / f"{a}.hits.tsv.gz", "rt"), delimiter="\t"):
            if h["n_fs"] != "0":
                continue
            L = len(seqs[h["contig"]])
            g0, g1 = int(h["gstart"]), int(h["gend"])
            lo, hi = min(g0, g1) - 1, max(g0, g1)
            s = st[h["strand"]][h["contig"]]
            p0 = lo if h["strand"] == "+" else L - hi
            for p in range(p0 + 60, p0 + (hi - lo) - 60, 3):
                c = s[p:p + 3]
                n += 1
                if c in cnt:
                    cnt[c] += 1
        w.writerow([a, r["family"], r["genus"], n] + [f"{1000 * cnt[c] / max(n, 1):.3f}" for c in ARG] + [cnt["CGG"]])


if __name__ == "__main__":
    main()
