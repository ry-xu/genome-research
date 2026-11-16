"""Turnover of the skipped C between related genomes, read from orthologous sequence (not from breaks).
For each genome A: windows of 121 nt (gene strand) centred on
  cgg  - the skipped C of each clean +1 CGG break;
  ggn  - the first G of in-frame GGN codons not preceded by C, in frame-continuous hits (n_fs 0), same number as cgg;
  oth  - the first base of other in-frame codons in those hits, same number.
Each window is aligned (blastn -task dc-megablast) to every other genome of the set. Per (window, B) the best HSP that
covers query 31..91 with identity >= 70% is read at the centre:
  cgg: B has a base opposite the C ('kept', 'kept_C' if that base is C) or a gap ('lost');
  ggn/oth: B has extra base(s) between query 60 and 61 ('gain', with the inserted bases) or not ('none'); a gap
  opposite query 61 is 'del'.
  net: B's length minus A's over query 46..76 (+-15 nt): -1 for a lost C wherever the aligner put the gap.
usage: python -I src/analysis/site_turnover.py windows <set.tsv> <outdir> <genome_dir> <out.fa> [seed]
       python -I src/analysis/site_turnover.py classify <blast.tsv> <set.tsv> -> TSV on stdout
blast.tsv: blastn -task dc-megablast -outfmt "6 qseqid sseqid pident qstart qend bitscore qseq sseq" against a db whose
sequence ids are <genome acc>__<contig>.
"""
import csv
import gzip
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from candidate import read_fasta  # noqa: E402
from summarize import flag  # noqa: E402

H = 60
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def windows(setf, outdir, gdir, outfa, seed=1):
    rnd = random.Random(seed)
    out = open(outfa, "w")
    for r in csv.DictReader(open(setf), delimiter="\t"):
        a = r["gca"]
        fsf = outdir / f"{a}.fs.tsv.gz"
        if not fsf.exists():
            continue
        seqs = read_fasta(str(gdir / a[:8] / f"{a}.fna.gz"))
        strands = {"+": {k: v.upper() for k, v in seqs.items()}}
        strands["-"] = {k: rc(v) for k, v in strands["+"].items()}
        fs = list(csv.DictReader(gzip.open(fsf, "rt"), delimiter="\t"))
        flag(fs)
        n = 0
        for x in fs:
            if not (x["clean"] and x["type"] == "+1" and x["old_codon"] == "CGG"):
                continue
            s = strands[x["strand"]][x["contig"]]
            c = x["ctx"].upper()
            i = s.find(c)
            if i < 0 or s.find(c, i + 1) >= 0:
                continue
            p = i + sum(ch.islower() for ch in x["ctx"])
            w = s[p - H:p + H + 1]
            if len(w) == 2 * H + 1 and w[H:H + 3] == "CGG":
                out.write(f">{a}|cgg|{x['contig']}|{x['strand']}|{p}\n{w}\n")
                n += 1
        ggn, oth = [], []
        for h in csv.DictReader(gzip.open(outdir / f"{a}.hits.tsv.gz", "rt"), delimiter="\t"):
            if h["n_fs"] != "0":
                continue
            L = len(strands["+"][h["contig"]])
            g0, g1 = int(h["gstart"]), int(h["gend"])
            lo, hi = min(g0, g1) - 1, max(g0, g1)          # 0-based half-open on + strand
            s = strands[h["strand"]][h["contig"]]
            st = lo if h["strand"] == "+" else L - hi        # hit start on its own strand, in frame
            for p in range(st + 3 * 20, st + (hi - lo) - 3 * 20, 3):
                cod = s[p:p + 3]
                if cod[:2] == "GG" and s[p - 1] != "C":
                    ggn.append((h["contig"], h["strand"], p))
                elif cod[:2] != "GG" and cod != "CGG":
                    oth.append((h["contig"], h["strand"], p))
        for lab, pool in (("ggn", ggn), ("oth", oth)):
            for ctg, sd, p in rnd.sample(pool, min(n, len(pool))):
                w = strands[sd][ctg][p - H:p + H + 1]
                if len(w) == 2 * H + 1:
                    out.write(f">{a}|{lab}|{ctg}|{sd}|{p}\n{w}\n")


def classify(blast, setf):
    meta = {r["gca"]: r for r in csv.DictReader(open(setf), delimiter="\t")}
    best = {}
    for r in csv.reader(open(blast), delimiter="\t"):
        q, sid, pid, qs, qe, bits, qseq, sseq = r[0], r[1], float(r[2]), int(r[3]), int(r[4]), float(r[5]), r[6], r[7]
        b = sid.split("__")[0]
        a = q.split("|")[0]
        if b == a or qs > 31 or qe < 91 or pid < 70:
            continue
        k = (q, b)
        if k not in best or bits > best[k][0]:
            best[k] = (bits, pid, qs, qseq, sseq)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["a", "b", "relation", "kind", "pid", "state", "detail", "net", "q"])
    for (q, b), (bits, pid, qs, qseq, sseq) in best.items():
        a, kind = q.split("|")[:2]
        qi, j = qs - 1, 0                                    # qi: 0-based query index of column j
        while j < len(qseq) and not (qseq[j] != "-" and qi == H):
            if qseq[j] != "-":
                qi += 1
            j += 1
        if j >= len(qseq):
            continue
        net, qi2 = 0, qs - 1
        for qc, sc in zip(qseq, sseq):
            if H - 15 <= qi2 <= H + 15:
                net += (qc == "-") - (sc == "-")
            if qc != "-":
                qi2 += 1
        ins, k = "", j - 1                                   # bases in B just before the centre column
        while k >= 0 and qseq[k] == "-":
            ins = sseq[k] + ins
            k -= 1
        if kind == "cgg":
            state = "lost" if sseq[j] == "-" else ("kept_C" if sseq[j] == "C" else "kept")
            detail = sseq[j] + ("+" + ins if ins else "")
        else:
            state = "del" if sseq[j] == "-" else ("gain" if ins else "none")
            detail = ins
        ma, mb = meta.get(a, {}), meta.get(b, {})
        rel = ("species" if ma.get("species") == mb.get("species") else "genus" if ma.get("genus") == mb.get("genus")
               else "family" if ma.get("family") == mb.get("family") else "order")
        w.writerow([a, b, rel, kind, f"{pid:.1f}", state, detail, net, q])


if __name__ == "__main__":
    if sys.argv[1] == "windows":
        windows(sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4]), sys.argv[5],
                int(sys.argv[6]) if len(sys.argv) > 6 else 1)
    else:
        classify(sys.argv[2], sys.argv[3])
