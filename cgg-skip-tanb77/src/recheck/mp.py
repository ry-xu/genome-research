"""Parse miniprot PAF (protein -> genome, -S no splicing) into per-codon alignment columns and frameshift events.

miniprot CIGAR (README): M n = n codons vs n aa; I n = n aa with no genome codon; D n = n genome codons facing
a protein gap; F n = n genome nt skipped (frameshift); G n = n genome nt (1-2) aligned to one aa (frameshift).
Frame change d = (genome nt consumed - 3 * aa consumed) mod 3:  F1, G1 -> d=1 (genome has one extra nt: the
ribosome must skip one nt = +1 frameshift);  F2, G2 -> d=2 (genome lacks one nt: -1 frameshift).

Coordinates: 0-based on the + strand. A codon locus is (contig, strand, start) with start = lowest + coordinate.
"Oriented" offsets run 5'->3' along the coding strand from the alignment start.
"""
import re
from collections import defaultdict

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
STD = dict(zip(CODONS, AA))
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.upper())
    if name: seqs[name] = "".join(buf)
    return seqs


def parse_paf(paf, genome, prots, win=7):
    """Yield one dict per alignment with M columns and frameshift events (oriented offsets)."""
    for line in open(paf):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 12: continue
        qn, qlen, qs, qe, strand, tn, tlen, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8])
        tags = {t[:2]: t[5:] for t in f[12:]}
        if "cg" not in tags: continue
        AS = int(tags.get("AS", 0))
        seg = genome[tn][ts:te]
        if strand == "-": seg = rc(seg)
        prot = prots[qn]
        cols = []   # (type, off, qpos)  type in M, I, D, F, G
        events = []  # (off, d, oplen, op)
        evq = []
        off, q = 0, qs
        for n, op in re.findall(r"(\d+)([MIDFGNUV])", tags["cg"]):
            n = int(n)
            if op == "M":
                for _ in range(n):
                    cols.append(("M", off, q)); off += 3; q += 1
            elif op == "I":
                for _ in range(n):
                    cols.append(("I", off, q)); q += 1
            elif op == "D":
                for _ in range(n):
                    cols.append(("D", off, None)); off += 3
            elif op == "F":
                evq.append(q); events.append((off, n % 3, n, "F")); cols.append(("F", off, None)); off += n
            elif op == "G":
                evq.append(q); events.append((off, (n - 3) % 3, n, "G")); cols.append(("G", off, q)); off += n; q += 1
            else:  # introns should not occur with -S
                cols.append(("N", off, None)); off += n
        assert off == te - ts, (qn, off, te - ts)
        # per-M-column local identity (+-win M columns, excluding self) and nearby gap flag (+-3 columns)
        Mi = [i for i, c in enumerate(cols) if c[0] == "M"]
        fa = []
        for i in Mi:
            _, o, qq = cols[i]
            fa.append(STD.get(seg[o:o + 3], "X"))
        ident = [1 if fa[k] == prot[cols[Mi[k]][2]] else 0 for k in range(len(Mi))]
        cum = [0]
        for x in ident: cum.append(cum[-1] + x)
        recs = []
        gapcol = [c[0] != "M" for c in cols]
        gcum = [0]
        for x in gapcol: gcum.append(gcum[-1] + x)
        for k, i in enumerate(Mi):
            lo, hi = max(0, k - win), min(len(Mi), k + win + 1)
            nwin = hi - lo - 1
            wid = (cum[hi] - cum[lo] - ident[k]) / nwin if nwin > 0 else 0.0
            glo, ghi = max(0, i - 3), min(len(cols), i + 4)
            neargap = (gcum[ghi] - gcum[glo]) > 0
            _, o, qq = cols[i]
            recs.append((o, seg[o:o + 3], prot[qq], wid, neargap, k, len(Mi) - 1 - k))
        Dcols = [(c[1], seg[c[1]:c[1] + 3]) for c in cols if c[0] == "D"]
        yield dict(qn=qn, ref=qn.split("|")[0], qlen=qlen, qs=qs, qe=qe, strand=strand, contig=tn, ts=ts, te=te,
                   AS=AS, seg=seg, recs=recs, events=events, evq=evq, Dcols=Dcols, nM=len(Mi))


def gcoord(a, off, length=3):
    """+ strand 0-based start of an oriented interval [off, off+length)."""
    return a["ts"] + off if a["strand"] == "+" else a["te"] - off - length


def best_per_ref(alns):
    """For each (locus, reference genome) keep the column from the highest-scoring alignment."""
    best = {}
    for ai, a in enumerate(alns):
        for (o, codon, raa, wid, neargap, kfrom, kto) in a["recs"]:
            key = (a["contig"], a["strand"], gcoord(a, o), a["ref"])
            cur = best.get(key)
            if cur is None or a["AS"] > cur[0]:
                best[key] = (a["AS"], ai, o, codon, raa, wid, neargap, kfrom, kto)
    return best
