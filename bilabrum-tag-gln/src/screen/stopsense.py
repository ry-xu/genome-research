"""Stop codons read as sense: what reference residue each in-frame codon aligns to, per assembly.

Input is new-stops-euk's stored nuclear scans (read-only): LABEL.ends.tsv holds, for every best DIAMOND hit that
reaches within 2 residues of its Swiss-Prot reference's C-terminus, the in-frame nucleotide body of the aligned span;
LABEL.codons.npz holds in-frame codon counts over all best-hit spans (N_in, 64).
Bodies are re-aligned here (DIAMOND blastx, the same settings new-stops-euk used, plus btop) in batches of units, and
every aligned column >= EDGE columns from both alignment ends is tallied as (codon, reference residue).

usage:
  stopsense.py scan UNITS.tsv OUTDIR THREADS BATCH   (UNITS: unit, stratum, prefix, name, group; resumable)
Writes OUTDIR/UNIT.npz (M: 64 x 26 counts, codon x reference letter; N_in; bodies, hits) and
OUTDIR/UNIT.stops.tsv.gz (one row per inner aligned stop: contig, sseqid, codon, residue, column, alignment length,
pident).
"""
import gzip
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

TOK = re.compile(r"(\d+)|(..)")  # btop tokens: a run of identities, or a query/subject pair ('*' = query stop)
REF = "lab/work/new-stops-euk/data/ref/euk_ref.dmnd"
B = "TCAG"
CODONS = [a + b + c for a in B for b in B for c in B]
CI = {c: i for i, c in enumerate(CODONS)}
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE1 = dict(zip(CODONS, AA))
LET = "ACDEFGHIKLMNPQRSTVWYUXBZJO"
LI = {c: i for i, c in enumerate(LET)}
STOPS = ("TAA", "TAG", "TGA")
EDGE = 10


def bodies(prefix):
    with open(prefix + ".ends.tsv") as f:
        next(f)
        for k, line in enumerate(f):
            x = line.rstrip("\n").split("\t")
            yield k, x[0], x[1], x[7]


def tally(unit, prefix, hits, out):
    body = {k: (contig, s, b) for k, contig, s, b in bodies(prefix)}
    M = np.zeros((64, len(LET)), np.int64)
    rows, nh = [], 0
    for h in hits:
        q, s, qs, qe, ss, se, sl, pid, fr, btop = h
        k = int(q.split("|", 1)[1])
        contig, s0, b = body[k]
        qs, fr = int(qs), int(fr)
        if fr != 1 or (qs - 1) % 3:
            continue  # bodies are in frame; skip off-frame or reverse hits
        nh += 1
        cols = []  # (codon, subject letter) per aligned column, None for gap columns
        pos = qs - 1
        for n, p in TOK.findall(btop):
            if n:
                for _ in range(int(n)):
                    c = b[pos:pos + 3]
                    cols.append((c, CODE1.get(c, "X")))
                    pos += 3
            elif p[0] == "-":
                cols.append(None)  # gap in query: subject residue only
            elif p[1] == "-":
                cols.append(None)
                pos += 3
            else:
                cols.append((b[pos:pos + 3], p[1]))
                pos += 3
        L = len(cols)
        for j in range(EDGE, L - EDGE):
            x = cols[j]
            if x is None or x[0] not in CI:
                continue
            c, a = x
            M[CI[c], LI.get(a, LI["X"])] += 1
            if c in STOPS:
                rows.append(f"{contig}\t{s}\t{c}\t{a}\t{j}\t{L}\t{pid}\n")
    N_in = np.load(prefix + ".codons.npz")["N_in"]
    np.savez_compressed(os.path.join(out, unit + ".npz"), M=M, N_in=N_in, bodies=len(body), hits=nh)
    with gzip.open(os.path.join(out, unit + ".stops.tsv.gz"), "wt") as f:
        f.write("contig\tsseqid\tcodon\tresidue\tcolumn\talnlen\tpident\n")
        f.writelines(rows)


def scan(units_tsv, out, threads, batch):
    os.makedirs(out, exist_ok=True)
    units = [l.rstrip("\n").split("\t") for l in open(units_tsv) if l.strip() and not l.startswith("unit\t")]
    todo = [u for u in units if not os.path.exists(os.path.join(out, u[0] + ".npz"))]
    for i in range(0, len(todo), batch):
        part = todo[i:i + batch]
        with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
            fa = os.path.join(d, "q.fa")
            with open(fa, "w") as f:
                for u in part:
                    for k, _, _, b in bodies(u[2]):
                        if len(b) >= 3 * (2 * EDGE + 1):
                            f.write(f">{u[0]}|{k}\n{b}\n")
            hf = os.path.join(d, "h.tsv")
            subprocess.run(["diamond", "blastx", "-q", fa, "-d", REF, "-o", hf, "--fast", "-k", "1", "--max-hsps", "1",
                            "-e", "1e-5", "--threads", str(threads), "--quiet", "--outfmt", "6", "qseqid", "sseqid",
                            "qstart", "qend", "sstart", "send", "slen", "pident", "qframe", "btop"], check=True)
            by = {u[0]: [] for u in part}
            for line in open(hf):
                h = line.rstrip("\n").split("\t")
                by[h[0].split("|", 1)[0]].append(h)
            for u in part:
                tally(u[0], u[2], by[u[0]], out)
        print(f"batch {i // batch + 1}: {len(part)} units done", file=sys.stderr, flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "scan":
        scan(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
