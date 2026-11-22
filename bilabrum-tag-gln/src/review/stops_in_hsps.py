"""Independent stop-in-HSP scan from tblastn outfmt 6 (standard code, -db_gencode 1).

usage: python3 -I stops_in_hsps.py GENOME.fna HITS.tsv OUTPREFIX

HITS columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen sframe slen qseq sseq
Writes OUTPREFIX.sites.tsv (one row per distinct genomic stop site, from the highest-bitscore HSP covering it),
OUTPREFIX.ends.tsv (downstream terminator after HSPs reaching the query C-terminus),
OUTPREFIX.hsps.tsv (per-locus best HSP with GC3 and TAG/TAA/TGA inner counts) and prints a summary.
"""
import sys
from collections import Counter, defaultdict

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = ("TAA", "TAG", "TGA")
B = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[i] for i, (a, b, c) in enumerate((a, b, c) for a in B for b in B for c in B)}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def rc(s):
    return s.translate(COMP)[::-1]


def codon_at(seq, strand, p):
    """p = 0-based forward coordinate of the codon's first base in reading direction."""
    if strand > 0:
        return seq[p:p + 3]
    return rc(seq[p - 2:p + 1])


def walk(h, seq):
    """yield (col, qpos, qres, sres, codon, nt_pos_first_base) for each alignment column with a subject residue."""
    strand = 1 if h["sframe"] > 0 else -1
    p = h["sstart"] - 1
    q = h["qstart"]
    for col, (a, b) in enumerate(zip(h["qseq"], h["sseq"])):
        if b != "-":
            yield col, (q if a != "-" else None), a, b, codon_at(seq, strand, p), p
            p += 3 * strand
        if a != "-":
            q += 1


def main(genome, hits, out):
    seqs = read_fasta(genome)
    H = []
    for line in open(hits):
        f = line.rstrip("\n").split("\t")
        h = dict(q=f[0], s=f[1], pid=float(f[2]), qstart=int(f[6]), qend=int(f[7]), sstart=int(f[8]), send=int(f[9]),
                 evalue=float(f[10]), bits=float(f[11]), qlen=int(f[12]), sframe=int(f[13]), qseq=f[15], sseq=f[16])
        h["s"] = h["s"].split("|")[1] if "|" in h["s"] else h["s"]
        H.append(h)
    # per genomic codon position, the best-scoring HSP that covers it (deduplicates paralogous queries)
    best = {}
    chk = Counter()
    for i, h in enumerate(H):
        seq = seqs[h["s"]]
        for col, qp, a, b, cod, p in walk(h, seq):
            if b in "ACDEFGHIKLMNPQRSTVWY*":
                chk[CODE.get(cod, "?") == b] += 1
            key = (h["s"], p, h["sframe"] > 0)
            if key not in best or H[best[key]]["bits"] < h["bits"]:
                best[key] = i
    owned = Counter(best.values())
    sites = []
    ncod = Counter()           # codon -> count at owned, ungapped columns
    pair = defaultdict(Counter)  # codon -> query residue counts
    hsp_rows = []
    for i, h in enumerate(H):
        if owned[i] == 0:
            continue
        seq = seqs[h["s"]]
        cols = list(walk(h, seq))
        nsub = len(cols)
        g3 = n3 = 0
        inner = Counter()
        for k, (col, qp, a, b, cod, p) in enumerate(cols):
            if best[(h["s"], p, h["sframe"] > 0)] != i:
                continue
            if len(cod) == 3 and cod[2] in "GCAT":
                n3 += 1
                g3 += cod[2] in "GC"
            if a != "-":
                ncod[cod] += 1
                pair[cod][a] += 1
            if cod in STOPS:
                dist = min(k, nsub - 1 - k)
                lo, hi = max(0, col - 12), col + 13
                win_q, win_s = h["qseq"][lo:hi], h["sseq"][lo:hi]
                idn = sum(x == y and x != "-" for x, y in zip(win_q, win_s))
                gaps_q = h["qseq"][max(0, col - 6):col + 7].count("-")
                start = p if h["sframe"] > 0 else p - 2
                sites.append(dict(contig=h["s"], pos1=start + 1, strand="+" if h["sframe"] > 0 else "-", codon=cod,
                                  qres=a, query=h["q"].split("|")[2] if "|" in h["q"] else h["q"], qpos=qp,
                                  hsp_q=f"{h['qstart']}-{h['qend']}/{h['qlen']}", pid=h["pid"], bits=h["bits"],
                                  edge_dist=dist, win_ident=idn, qgaps6=gaps_q,
                                  ctx_nt=(seq[start - 12:start + 15] if h["sframe"] > 0 else rc(seq[start - 12:start + 15])),
                                  q_win=win_q, s_win=win_s))
                if dist >= 10:
                    inner[cod] += 1
        hsp_rows.append((h["q"].split("|")[2] if "|" in h["q"] else h["q"], h["s"], h["sstart"], h["send"],
                         h["qstart"], h["qend"], h["qlen"], h["pid"], h["bits"], nsub, round(g3 / max(n3, 1), 3),
                         inner["TAG"], inner["TAA"], inner["TGA"]))
    with open(out + ".sites.tsv", "w") as fh:
        keys = list(sites[0].keys()) if sites else []
        fh.write("\t".join(keys) + "\n")
        for s in sites:
            fh.write("\t".join(str(s[k]) for k in keys) + "\n")
    with open(out + ".pairs.tsv", "w") as fh:
        fh.write("codon\tresidue\tn\n")
        for c in sorted(pair):
            for a, n in sorted(pair[c].items()):
                fh.write(f"{c}\t{a}\t{n}\n")
    with open(out + ".hsps.tsv", "w") as fh:
        fh.write("query\tcontig\tsstart\tsend\tqstart\tqend\tqlen\tpid\tbits\tncod\tgc3\tTAG_inner\tTAA_inner\tTGA_inner\n")
        for r in hsp_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    # gene ends: HSPs (any bitscore rank) reaching within 3 residues of the query C-terminus; per locus keep best.
    ends = {}
    for i, h in enumerate(H):
        if h["qend"] < h["qlen"] - 3 or h["evalue"] > 1e-10:
            continue
        seq = seqs[h["s"]]
        strand = 1 if h["sframe"] > 0 else -1
        # 0-based reading-direction first base of the codon after the last aligned codon.
        # forward: last codon is 1-based send-2..send, next starts at 0-based send.
        # reverse: last codon is 1-based send..send+2, next codon is 1-based send-3..send-1, first base (highest) 0-based send-2.
        pp = h["send"] if strand > 0 else h["send"] - 2
        key = (h["s"], pp, strand)
        if key in ends and ends[key]["bits"] >= h["bits"]:
            continue
        cods = []
        for k in range(200):
            c = codon_at(seq, strand, pp + 3 * k * strand)
            if len(c) < 3:
                break
            cods.append(c)
        first = next(((k, c) for k, c in enumerate(cods) if c in STOPS), (None, None))
        first_nontag = next(((k, c) for k, c in enumerate(cods) if c in ("TAA", "TGA")), (None, None))
        tags_before = sum(c == "TAG" for c in cods[:first_nontag[0]]) if first_nontag[0] is not None else None
        ends[key] = dict(query=h["q"].split("|")[2] if "|" in h["q"] else h["q"], contig=h["s"], send=h["send"],
                         strand="+" if strand > 0 else "-", qend=h["qend"], qlen=h["qlen"], bits=h["bits"],
                         first_stop=first[1], first_stop_dist=first[0], first_TAA_TGA=first_nontag[1],
                         first_TAA_TGA_dist=first_nontag[0], tags_before=tags_before,
                         down=" ".join(cods[:12]))
    with open(out + ".ends.tsv", "w") as fh:
        keys = ["query", "contig", "send", "strand", "qend", "qlen", "bits", "first_stop", "first_stop_dist",
                "first_TAA_TGA", "first_TAA_TGA_dist", "tags_before", "down"]
        fh.write("\t".join(keys) + "\n")
        for e in ends.values():
            fh.write("\t".join(str(e[k]) for k in keys) + "\n")
    # summary
    print(f"coordinate check: genome codon translates to the HSP residue {chk[True]}, mismatches {chk[False]}")
    print(f"HSPs {len(H)}; owning HSPs {len(hsp_rows)}; owned codon columns {sum(ncod.values())}")
    for c in STOPS:
        allc = [s for s in sites if s["codon"] == c]
        inn = [s for s in allc if s["edge_dist"] >= 10]
        print(f"{c}: all {len(allc)} (in {len({s['query'] for s in allc})} queries); inner>=10 {len(inn)} "
              f"(in {len({s['query'] for s in inn})} queries); inner residues {Counter(s['qres'] for s in inn).most_common()}")
        print(f"    all residues {Counter(s['qres'] for s in allc).most_common()}")
    for c in ("CAG", "CAA", "TGG", "TAC", "TAT"):
        n = sum(pair[c].values())
        print(f"{c} n={n} top {pair[c].most_common(5)}")
    tot = sum(ncod.values())
    for c in ("TAG", "CAG", "CAA"):
        print(f"{c} per 1000 aligned codons: {1000 * ncod[c] / tot:.2f}")
    ec = Counter(e["first_stop"] for e in ends.values())
    print(f"gene ends (HSP qend>=qlen-3, e<=1e-10), {len(ends)} loci; first in-frame stop codon: {dict(ec)}")
    near = [e for e in ends.values() if e["first_TAA_TGA_dist"] is not None and e["first_TAA_TGA_dist"] <= 20]
    print(f"  loci with TAA/TGA within 20 codons: {len(near)}; terminator {dict(Counter(e['first_TAA_TGA'] for e in near))}; "
          f"TAGs read before it {sum(e['tags_before'] for e in near)} in {sum(e['tags_before'] > 0 for e in near)} loci")
    near0 = [e for e in ends.values() if e["first_stop_dist"] is not None and e["first_stop_dist"] <= 3]
    print(f"  loci with any stop within 0-3 codons after HSP: {len(near0)}; {dict(Counter(e['first_stop'] for e in near0))}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
