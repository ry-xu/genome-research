"""Arm C: stop-codon recognition motifs of every release factor found by the GTDB scan.

Each RF1, RF2 and aRF1 hit in results/gtdb/rf_hits.tsv is aligned (local, BLOSUM62) to its E. coli or M. jannaschii
reference and the residues at the reference's recognition-motif positions are read off:
  RF1   PxT     E. coli RF1 188-190 (PAT)        canonical: P.T
  RF2   SPF     E. coli RF2 205-207 (SPF)        canonical: SP[FY]
  aRF1  NIKS    M. jannaschii aRF1 56-59         canonical: NIKS
        YxCxxxF M. jannaschii aRF1 121-127       canonical: Y.C...F
plus a 21-residue window centred on the motif (for the Sphaeropleales mtRF1a-like change near PxT, reported only).
A genome whose RF copies all lack the canonical motif goes to Arm B (PREREG.md, Arm C extra entry).

usage: uv run --no-project --with biopython python -I src/rf_motifs.py [rf_hits.tsv] [out.tsv]
"""
import gzip
import re
import sys
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parent.parent
MARKERS = ROOT.parent / "tag-bacteria" / "data" / "reference" / "markers.faa.gz"
REFS = {"RF1": "RF1|RF1_ECOLI", "RF2": "RF2|RF2_ECOLI", "aRF1": "aRF1|RF1_METJA"}
# (name, first reference position 1-based, length, canonical regex)
MOTIFS = {"RF1": [("PxT", 188, 3, r"P.T")], "RF2": [("SPF", 205, 3, r"SP[FY]")],
          "aRF1": [("NIKS", 56, 4, r"NIKS"), ("YxCxxxF", 121, 7, r"Y.C...F")]}
HALF = 10


def read_faa(fn):
    seqs, name = {}, None
    for line in (gzip.open(fn, "rt") if str(fn).endswith(".gz") else open(fn)):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def aligner():
    a = Align.PairwiseAligner()
    a.mode = "local"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -10, -0.5
    return a


def ref_to_query(aln):
    """-> {reference index (0-based): query index (0-based)} over aligned columns."""
    m = {}
    for (rs, re_), (qs, qe) in zip(*aln.aligned):
        for k in range(re_ - rs):
            m[rs + k] = qs + k
    return m


def read_off(q, m, start, length):
    return "".join(q[m[i]] if i in m else "-" for i in range(start - 1, start - 1 + length))


def main(hits_fn, out_fn):
    refs = {k: read_faa(MARKERS)[v] for k, v in REFS.items()}
    for k, mot in MOTIFS.items():  # the reference must show its own canonical motif at the stated position
        for name, st, ln, rx in mot:
            assert re.fullmatch(rx, refs[k][st - 1:st - 1 + ln]), (k, name)
    al = aligner()
    out = open(out_fn, "w")
    out.write("acc\tgene\tclass\tpident\tlen\tmotif\tresidues\tcanonical\twindow\tref_cover\n")
    seen = set()
    for line in open(hits_fn):
        acc, gene, cls, pid, ln, seq = line.rstrip("\n").split("\t")[:6]
        if cls not in REFS or (acc, gene) in seen:
            continue
        seen.add((acc, gene))
        r = refs[cls]
        aln = al.align(r, seq)[0]
        m = ref_to_query(aln)
        cover = len(m) / len(r)
        for name, st, ml, rx in MOTIFS[cls]:
            res = read_off(seq, m, st, ml)
            win = read_off(seq, m, st - HALF, 2 * HALF + ml)
            out.write(f"{acc}\t{gene}\t{cls}\t{pid}\t{ln}\t{name}\t{res}\t{int(bool(re.fullmatch(rx, res)))}\t{win}\t"
                      f"{cover:.6f}\n")
    out.close()


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if a else ROOT / "results" / "gtdb" / "rf_hits.tsv",
         a[1] if len(a) > 1 else ROOT / "results" / "gtdb" / "rf_motifs.tsv")
