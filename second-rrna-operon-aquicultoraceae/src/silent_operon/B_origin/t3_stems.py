"""Processing stems: do the sequences flanking each 16S and 23S gene pair into a long stem (the RNase III
processing stem of bacterial pre-rRNA)?  Best local complementary alignment of the 5' flank window against the
reverse complement of the 3' flank window (WC +2, GU +1, mismatch -3, gap -5), against 20 mononucleotide shuffles.
Windows: [gene_start-W5, gene_start+40] and [gene_end-40, gene_end+W3] in the gene's orientation.
Usage: python -I src/t3_stems.py  -> out/t3_stems.tsv
"""
import csv, random
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

random.seed(1)
import sys
W = int(sys.argv[1]) if len(sys.argv) > 1 else 250
OUT = sys.argv[2] if len(sys.argv) > 2 else "out/t3_stems.tsv"
FAMONLY = len(sys.argv) > 3
gen = SeqIO.index("data/aq65.fna", "fasta")
ctlg = SeqIO.index("data/ctl/ctl.fna", "fasta")

PAIR = {("A", "T"): 2, ("T", "A"): 2, ("G", "C"): 2, ("C", "G"): 2, ("G", "T"): 1, ("T", "G"): 1}


def sw(a, b):
    """a 5'->3' (5' flank); b 5'->3' (3' flank). Pair a[i] with b reversed."""
    br = b[::-1]
    n, m = len(a), len(br)
    S = np.full((len(a), len(br)), -3, dtype=np.int32)
    A = np.frombuffer(a.encode(), dtype=np.uint8)[:, None]
    B = np.frombuffer(br.encode(), dtype=np.uint8)[None, :]
    for (x, y), v in PAIR.items():
        S[(A == ord(x)) & (B == ord(y))] = v
    H = np.zeros((n + 1, m + 1), dtype=np.int32)
    L = np.zeros((n + 1, m + 1), dtype=np.int32)  # paired positions in path
    best = (0, 0)
    for i in range(1, n + 1):
        diag = H[i - 1, :-1] + S[i - 1]
        up = H[i - 1, 1:] - 5
        row = np.maximum(0, np.maximum(diag, up))
        # left gaps sequentially
        Ld = L[i - 1, :-1] + 1; Lu = L[i - 1, 1:]
        lrow = np.where(row == diag, Ld, np.where(row == up, Lu, 0))
        for j in range(1, m + 1):
            v = row[j - 1]
            g = H[i, j - 1] - 5
            if g > v:
                v = g; lrow[j - 1] = L[i, j - 1]
            if v <= 0:
                v = 0; lrow[j - 1] = 0
            H[i, j] = v; L[i, j] = lrow[j - 1]
        k = int(np.argmax(H[i]))
        if H[i, k] > best[0]:
            best = (int(H[i, k]), int(L[i, k]))
    return best


def window(rec, start, end, strand):
    L = len(rec.seq)
    if strand == "+":
        a0, a1, b0, b1 = start - W, start + 40, end - 40, end + W
    else:
        a0, a1, b0, b1 = end - 40, end + W, start - W, start + 40
    if min(a0, b0) < 1 or max(a1, b1) > L:
        return None, None
    a = rec.seq[a0 - 1:a1]; b = rec.seq[b0 - 1:b1]
    if strand == "-":
        a, b = a.reverse_complement(), b.reverse_complement()
    return str(a).upper(), str(b).upper()


def shuf(s):
    l = list(s); random.shuffle(l); return "".join(l)


items = []
for r in csv.DictReader(open("out/loci.tsv"), delimiter="\t"):
    if r["family"] == "Aquicultoraceae" and float(r["frac"]) >= 0.85:
        items.append((r["id"], r["kind"], r["genus"], r["label"], gen[r["contig"]], int(r["start"]), int(r["end"]), r["strand"]))
for fn, kind in ([] if FAMONLY else [("out/ctl_ssu.fna", "16S"), ("out/ctl_lsu.fna", "23S")]):
    for rec in SeqIO.parse(fn, "fasta"):
        org, _, rest = rec.id.split("_", 2)
        if org not in ("ec", "tt", "tb", "vv"):
            continue
        ctg, pos = rest.rsplit("_", 1)
        st = pos[-1]; s = int(pos[:-1]); e = s + len(rec.seq) - 1
        if len(rec.seq) < (1400 if kind == "16S" else 2700):
            continue
        items.append((rec.id, kind, org, "ctl", ctlg[f"{org}|{ctg}"], s, e, st))

out = open(OUT, "w")
out.write("id\tkind\tgenus\tlabel\tscore\tpaired\tnull_mean\tnull_max\tz\n")
seen = set()
for gid, kind, genus, lab, rec, s, e, st in items:
    a, b = window(rec, s, e, st)
    if a is None:
        out.write(f"{gid}\t{kind}\t{genus}\t{lab}\tNA\tNA\tNA\tNA\tNA\n"); continue
    key = (a, b)
    if key in seen:  # identical copies
        continue
    seen.add(key)
    sc, pl = sw(a, b)
    null = [sw(shuf(a), shuf(b))[0] for _ in range(20)]
    mu = float(np.mean(null)); sd = float(np.std(null)) or 1
    out.write(f"{gid}\t{kind}\t{genus}\t{lab}\t{sc}\t{pl}\t{mu:.1f}\t{max(null)}\t{(sc - mu) / sd:.1f}\n")
    out.flush()
