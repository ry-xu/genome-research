"""Frame breaks in homolog alignments: DIAMOND blastx -F 15 of whole genomes against clustered Swiss-Prot.

usage: [THREADS=4] [GENCODE=11] [FS_DB=data/sprot_p50.dmnd] python -I src/fs_scan.py <outdir> <name>=<fasta[.gz]> ...
  Genomes are searched in one DIAMOND batch (contig ids prefixed with <name>|), all with one genetic code (GENCODE,
  NCBI table number; batch genomes by code). Per genome writes
  <outdir>/<name>.hits.tsv.gz  one row per kept hit (locus): contig strand gstart gend sseqid sstart send slen bits
                               alen n_fs
  <outdir>/<name>.fs.tsv.gz    one row per frame break (see FS_COLS)
CIGAR semantics (checked on 13,686 alignments, scratch/bench): M and I consume 3 nt, D none, '/' steps the query
back 1 nt (genome lacks a nt: a -1 shift is needed to stay in frame), '\\' steps it forward 1 nt (genome has an extra
nt: a +1 shift, the Euplotes / prfB type).
"""
import bisect
import gzip
import os
import random
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB = Path(os.environ.get("FS_DB", ROOT / "data/sprot_p50.dmnd"))
DBFA = Path(os.environ.get("FS_DBFA", str(DB).replace(".dmnd", ".faa")))
THREADS = os.environ.get("THREADS", "4")
EVALUE = "1e-10"
KHITS = "10"      # DIAMOND hits kept per query range; culled to one below (whole-genome -k 1 loses ~18% of loci)
WIN = 15          # aligned columns each side used for flank identity
CTX = 20          # nt of context each side of the break
GENCODE = os.environ.get("GENCODE", "11")
CODE = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODE[a + b + c] = _aa[16 * i + 4 * j + k]
# the reassignments used by lineages in the scan (NCBI tables)
CODE.update({"4": {"TGA": "W"}, "10": {"TGA": "C"}, "15": {"TAG": "Q"}, "25": {"TGA": "G"}, "16": {"TAG": "L"},
             "6": {"TAA": "Q", "TAG": "Q"}}.get(GENCODE, {}))
STOPS = {c for c, a in CODE.items() if a == "*"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FS_COLS = ["contig", "strand", "gpos", "type", "sseqid", "spos", "slen", "bits", "aln_aa", "left_id", "right_id",
           "left_n", "right_n", "hp", "hp_base", "stop_old", "old_codon", "ctx",
           "k_lo", "k_hi", "site", "site_codon", "site_stop", "win_hp", "win_hp_base", "site_ctx", "a_stop", "stop_motif",
           "win_seq"]
WINSEQ = 30       # nt of gene strand from 3 codons before the tie window's left edge, for motif enrichment
NBG = 2           # random background windows of WINSEQ nt per aligned locus (same gene strand)
ASTOP_MAX = 30    # codons past the tie window searched for the first stop in the abandoned frame
KSLIDE = 8        # codons the break may slide each way when finding its tie window
SITE_CTX = 12     # nt each side of the refined site


def culled(path):
    """DIAMOND lines culled to one hit per range: a hit is dropped when a higher-scoring hit on the same query, either
    strand, covers at least half its query range (DIAMOND's own --range-culling rule). Run with -k 1, DIAMOND's
    target limit drops loci on whole-genome queries (E. coli: 3,359 loci; -k 10 culled here: 4,062, holding 3,342 of
    the 3,359 with the same start and subject)."""
    by = {}
    for i, line in enumerate(open(path)):
        r = line.rstrip("\n").split("\t")
        a, b = sorted((int(r[2]), int(r[3])))
        by.setdefault(r[0], []).append((-float(r[7]), i, a, b, line))
    out = []
    for hits in by.values():
        hits.sort()
        starts, kept, maxlen = [], [], 0
        for _, i, a, b, line in hits:
            lo = bisect.bisect_left(starts, a - maxlen)
            hi = bisect.bisect_right(starts, b)
            if all(min(b, kept[j][1]) - max(a, kept[j][0]) + 1 < (b - a + 1) / 2 for j in range(lo, hi)):
                j = bisect.bisect_right(starts, a)
                starts.insert(j, a)
                kept.insert(j, (a, b))
                maxlen = max(maxlen, b - a + 1)
                out.append((i, line))
    return [line for _, line in sorted(out)]


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


_subj = None


def subjects():
    global _subj
    if _subj is None:
        _subj = read_fasta(DBFA)
    return _subj


def homopolymer(s, p):
    """longest single-base run touching s[p-3 .. p+2] (DIAMOND places a break within about a codon of the true site)"""
    best, base = 0, ""
    for q in range(p - 3, p + 3):
        if q < 0 or q >= len(s):
            continue
        b = s[q]
        i = q
        while i > 0 and s[i - 1] == b:
            i -= 1
        j = q
        while j + 1 < len(s) and s[j + 1] == b:
            j += 1
        if j - i + 1 > best:
            best, base = j - i + 1, b
    return best, base


def walk(seq, strand, qstart, qend, sseq, sstart, cigar):
    """Return the gene-strand sequence, the breaks [(gene-strand nt index, type, subject pos, M-column index)], the
    M-column matches (1 = identical aa), the end coordinate, and the set of (codon index, subject index) M columns."""
    L = len(seq)
    s = seq if strand == "+" else rc(seq)
    pos = qstart - 1 if strand == "+" else L - qstart
    spos = sstart - 1
    cols = []          # 1 = identical aa, 0 = not, one per M column
    mset = set()
    breaks = []
    for n, op in re.findall(r"(\d+)([MID/\\])", cigar):
        n = int(n)
        if op == "M":
            for _ in range(n):
                aa = CODE.get(s[pos:pos + 3], "X")
                cols.append(1 if spos < len(sseq) and aa == sseq[spos] else 0)
                mset.add((pos, spos))
                pos += 3
                spos += 1
        elif op == "I":
            pos += 3 * n
        elif op == "D":
            spos += n
        elif op == "/":
            for _ in range(n):
                breaks.append((pos, "-1", spos, len(cols)))
                pos -= 1
        elif op == "\\":
            for _ in range(n):
                breaks.append((pos, "+1", spos, len(cols)))
                pos += 1
    end = pos if strand == "+" else L - pos + 1  # last nt consumed, forward coordinate (1-based)
    return s, breaks, cols, end, mset


_B62 = None


def b62(a, b):
    global _B62
    if _B62 is None:
        from Bio.Align import substitution_matrices
        _B62 = substitution_matrices.load("BLOSUM62")
    try:
        return _B62[a][b]
    except (KeyError, IndexError):
        return -1


def refine(s, p, d, sb, sseq, mset):
    """Tie window of a break. The break sits before subject residue sb, with frame A codons at p + 3j upstream and
    frame B codons at p + d + 3j downstream (d = +1 or -1). Sliding it to k puts residues j < k on frame A and j >= k on
    frame B; only residues aligned without gaps on the side being re-read are allowed. Returns (k_lo, k_hi) of the
    placements whose BLOSUM62 score ties the best."""
    def aa(i):
        return CODE.get(s[i:i + 3], "X") if 0 <= i and i + 3 <= len(s) else "X"
    best, ks = 0, {0: 0}
    sc = 0
    for k in range(1, KSLIDE + 1):         # slide right: residue sb+k-1 moves from frame B to frame A
        j = k - 1
        if (p + d + 3 * j, sb + j) not in mset or sb + j >= len(sseq):
            break
        sc += b62(aa(p + 3 * j), sseq[sb + j]) - b62(aa(p + d + 3 * j), sseq[sb + j])
        ks[k] = sc
    sc = 0
    for k in range(-1, -KSLIDE - 1, -1):   # slide left: residue sb+k moves from frame A to frame B
        j = k
        if (p + 3 * j, sb + j) not in mset or sb + j < 0:
            break
        sc += b62(aa(p + d + 3 * j), sseq[sb + j]) - b62(aa(p + 3 * j), sseq[sb + j])
        ks[k] = sc
    best = max(ks.values())
    tie = [k for k, v in ks.items() if v == best]
    return min(tie), max(tie)


def main():
    outdir = Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)
    genomes = [a.split("=", 1) for a in sys.argv[2:]]
    seqs = {}
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as td:
        q = Path(td) / "q.fna"
        with open(q, "w") as fh:
            for name, path in genomes:
                g = read_fasta(path)
                seqs[name] = g
                for c, s in g.items():
                    fh.write(f">{name}|{c}\n{s}\n")
        out = Path(td) / "hits.tsv"
        subprocess.run(["diamond", "blastx", "-q", str(q), "-d", str(DB), "-o", str(out), "--threads", THREADS,
                        "-F", "15", "--query-gencode", GENCODE, "-c", "1", "--range-culling", "-k", KHITS, "-e", EVALUE, "--quiet", "--outfmt", "6",
                        "qseqid", "sseqid", "qstart", "qend", "sstart", "send", "slen", "bitscore", "qframe",
                        "cigar"], check=True)
        sub = subjects()
        fh_h, fh_f, fh_b = {}, {}, {}
        for name, _ in genomes:
            fh_h[name] = gzip.open(outdir / f"{name}.hits.tsv.gz.part", "wt")
            fh_h[name].write("contig\tstrand\tgstart\tgend\tsseqid\tsstart\tsend\tslen\tbits\talen\tn_fs\n")
            fh_f[name] = gzip.open(outdir / f"{name}.fs.tsv.gz.part", "wt")
            fh_f[name].write("\t".join(FS_COLS) + "\n")
            fh_b[name] = gzip.open(outdir / f"{name}.bg.tsv.gz.part", "wt")
            fh_b[name].write("win_seq\n")
        bad = 0
        for line in culled(out):
            qid, sid, qs, qe, ss, se, slen, bits, fr, cig = line.rstrip("\n").split("\t")
            name, contig = qid.split("|", 1)
            seq = seqs[name][contig]
            strand = "+" if int(fr) > 0 else "-"
            sseq = sub.get(sid, "")
            s, breaks, cols, end, mset = walk(seq, strand, int(qs), int(qe), sseq, int(ss), cig)
            if end != int(qe):
                bad += 1
            g0, g1 = sorted((int(qs), int(qe)))
            # background windows: random codon-aligned starts inside the aligned span, gene strand
            a0 = int(qs) - 1 if strand == "+" else len(seq) - int(qs)
            span = abs(int(qe) - int(qs)) + 1
            rng = random.Random(f"{qid}:{qs}")
            for _ in range(NBG if span > WINSEQ + 3 else 0):
                b0 = a0 + 3 * rng.randrange((span - WINSEQ) // 3)
                fh_b[name].write(s[b0:b0 + WINSEQ] + "\n")
            fh_h[name].write(f"{contig}\t{strand}\t{g0}\t{g1}\t{sid}\t{ss}\t{se}\t{slen}\t{bits}\t{len(cols)}\t"
                             f"{len(breaks)}\n")
            for p, typ, spos, ci in breaks:
                left = cols[max(0, ci - WIN):ci]
                right = cols[ci:ci + WIN]
                hp, hb = homopolymer(s, p)
                # the abandoned frame: codons continuing from p in the old frame; a stop within 2 codons either side
                old = [s[p + 3 * k:p + 3 * k + 3] for k in range(-2, 3)]
                stop_old = ",".join(str(k - 2) for k, c in enumerate(old) if c in STOPS)
                ctx = s[max(0, p - CTX):p].lower() + s[p:p + CTX]
                gpos = p + 1 if strand == "+" else len(seq) - p
                k_lo, k_hi = refine(s, p, 1 if typ == "+1" else -1, spos, sseq, mset)
                site = p + 3 * k_hi       # rightmost tie: frame A would read its next codon here
                sc = s[site:site + 3]
                whp, whb = 0, ""
                for q in range(p + 3 * k_lo - 1, site + 4, 3):
                    h, b = homopolymer(s, q + 1)
                    if h > whp:
                        whp, whb = h, b
                # first stop in the abandoned frame A from the left edge of the tie window; its distance (codons) past the
                # rightmost tie, and the two frame-A codons before it plus the stop (the E- and P-site codons of a
                # Euplotes / prfB type site)
                a_stop, motif = "nan", ""
                for m in range(k_lo, k_hi + ASTOP_MAX + 1):
                    c = s[p + 3 * m:p + 3 * m + 3]
                    if len(c) < 3:
                        break
                    if c in STOPS:
                        a_stop = m - k_hi
                        motif = s[max(0, p + 3 * m - 6):p + 3 * m] + "-" + c
                        break
                fh_f[name].write("\t".join(map(str, [
                    contig, strand, gpos, typ, sid, spos + 1, slen, bits, len(cols),
                    f"{sum(left) / len(left):.2f}" if left else "nan", f"{sum(right) / len(right):.2f}" if right else "nan",
                    len(left), len(right), hp, hb, stop_old, s[p:p + 3], ctx,
                    k_lo, k_hi, site + 1 if strand == "+" else len(seq) - site, sc, int(sc in STOPS), whp, whb,
                    s[max(0, site - SITE_CTX):site].lower() + s[site:site + SITE_CTX], a_stop, motif,
                    s[max(0, p + 3 * k_lo - 9):max(0, p + 3 * k_lo - 9) + WINSEQ]])) + "\n")
        for name, _ in genomes:
            fh_h[name].close()
            fh_f[name].close()
            fh_b[name].close()
            for k in ("bg", "hits", "fs"):  # fs last: its presence marks the genome done
                os.replace(outdir / f"{name}.{k}.tsv.gz.part", outdir / f"{name}.{k}.tsv.gz")
        if bad:
            print(f"WARNING: {bad} hits whose CIGAR walk did not end at qend", file=sys.stderr)


if __name__ == "__main__":
    main()
