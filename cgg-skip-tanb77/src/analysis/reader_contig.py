"""Where does a reader-keeping TANB77 genome's CGG-reading tRNA sit? Per genome: ARAGORN (pyaragorn) Arg tRNAs with their
contig, contig length and GC, and per contig the clean +1 CGG breaks (scan table) and all tRNAs.
usage: python -I src/analysis/reader_contig.py <outdir> <gca> ..."""
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

import pyaragorn

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from candidate import read_fasta  # noqa: E402
from summarize import flag  # noqa: E402

from paths import GENOMES as GDIR  # noqa: E402


def main():
    outdir = Path(sys.argv[1])
    finder = pyaragorn.RNAFinder(11, tmrna=False, linear=True)
    for a in sys.argv[2:]:
        seqs = read_fasta(str(GDIR / a[:8] / f"{a}.fna.gz"))
        fs = list(csv.DictReader(gzip.open(outdir / f"{a}.fs.tsv.gz", "rt"), delimiter="\t"))
        flag(fs)
        cgg = Counter(r["contig"] for r in fs if r["clean"] and r["type"] == "+1" and r["old_codon"] == "CGG")
        p1 = Counter(r["contig"] for r in fs if r["clean"] and r["type"] == "+1")
        trn = {}
        for k, s in seqs.items():
            for g in finder.find_rna(s.encode() if isinstance(s, str) and False else s):
                trn.setdefault(k, []).append((g.amino_acid, g.anticodon.upper(), g.begin, g.end, g.strand))
        tot = sum(len(s) for s in seqs.values())
        gc = sum(s.upper().count("G") + s.upper().count("C") for s in seqs.values()) / tot
        print(f"== {a}: {len(seqs)} contigs, {tot:,} nt, GC {gc:.3f}, tRNAs {sum(map(len, trn.values()))}, "
              f"clean +1 {sum(p1.values())}, at CGG {sum(cgg.values())}")
        for k, s in sorted(seqs.items(), key=lambda x: -len(x[1])):
            arg = [t for t in trn.get(k, []) if t[0].startswith("Arg") or t[1] in ("CCG", "TCG", "ACG", "TCT", "CCT")]
            u = s.upper()
            cgc = (u.count("G") + u.count("C")) / len(u)
            if arg or cgg[k] or p1[k]:
                print(f"  {k}\tlen {len(s):,}\tGC {cgc:.3f}\ttRNAs {len(trn.get(k, []))}\tclean +1 {p1[k]}\tCGG {cgg[k]}\t"
                      + " ".join(f"{t[0]}-{t[1]}@{t[2]}" for t in arg))


if __name__ == "__main__":
    main()
