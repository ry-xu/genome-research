"""Read codons at genomic TAG sites from RNA-seq reads mapped with minimap2 (PAF, --cs=long).

usage: site_reads.py GENOME SITES.tsv ALN.paf OUTPREFIX
SITES.tsv is the reviewer's bil.sites.tsv (contig, pos1, strand, ...): pos1 is the 1-based leftmost base of the
codon on the + strand for both strands, so the codon is pos1..pos1+2 (reverse-complemented for - sites).
Only primary alignments (tp:A:P) with mapq >= 10 count. A read covers a site if all three bases are aligned with no
indel inside the codon. Writes OUTPREFIX.sites.tsv (per-site codon counts) and OUTPREFIX.summary.txt.
"""
import sys
from collections import Counter, defaultdict

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(p):
    seqs, name, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def cs_walk(cs, tstart):
    """Yield (target_pos0, read_base or None for deletion, ins_after) along the forward target strand."""
    i, t = 0, tstart
    while i < len(cs):
        op = cs[i]
        j = i + 1
        while j < len(cs) and cs[j] not in "=*+-~:":
            j += 1
        body = cs[i + 1:j]
        if op == "=":
            for b in body:
                yield t, b.upper(), False
                t += 1
        elif op == "*":
            yield t, body[1].upper(), False
            t += 1
        elif op == "-":
            for _ in body:
                yield t, None, False
                t += 1
        elif op == "+":
            yield t - 1, "", True
        else:
            raise ValueError(f"cs op {op}")
        i = j


def main(genome, sites_tsv, paf, out):
    seqs = read_fasta(genome)
    sites = []
    for line in open(sites_tsv):
        f = line.rstrip("\n").split("\t")
        if f[0] == "contig":
            hdr = f
            continue
        r = dict(zip(hdr, f))
        p = int(r["pos1"])
        lo = p - 1  # pos1 is the leftmost base on the + strand for both strands (stops_in_hsps.py: start = p - 2 on -)
        g = seqs[r["contig"]][lo:lo + 3]
        g = g if r["strand"] == "+" else rc(g)
        assert g == r["codon"], (r["contig"], p, r["strand"], g, r["codon"])
        sites.append((r["contig"], lo, r["strand"], r))
    by_contig = defaultdict(list)
    for k, (c, lo, s, r) in enumerate(sites):
        by_contig[c].append((lo, k))
    counts = [Counter() for _ in sites]
    ident = Counter()
    nreads = Counter()
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        tags = dict(x.split(":", 2)[0::2] for x in f[12:])
        if tags.get("tp") != "P" or int(f[11]) < 10:
            continue
        tname, tstart, tend = f[5], int(f[7]), int(f[8])
        kind = "chrom" if tname.startswith("OZ") else "unplaced"
        nreads[kind] += 1
        ident[(kind, "match")] += int(f[9])
        ident[(kind, "alen")] += int(f[10])
        hits = [(lo, k) for lo, k in by_contig.get(tname, []) if lo >= tstart and lo + 3 <= tend]
        if not hits:
            continue
        cov = {}
        bad = set()
        for t, b, ins in cs_walk(tags["cs"], tstart):
            if ins:
                bad.add(t)
            else:
                cov[t] = b
        for lo, k in hits:
            bs = [cov.get(lo + d) for d in range(3)]
            if None in bs or any(lo + d in bad for d in range(2)):
                counts[k]["indel"] += 1
                continue
            cod = "".join(bs)
            counts[k][cod if sites[k][2] == "+" else rc(cod)] += 1
    with open(out + ".sites.tsv", "w") as w:
        w.write("contig\tpos1\tstrand\tquery\tqpos\tqres\twin_ident\tedge_dist\treads\tTAG\tCAG\tother\n")
        for (c, lo, s, r), n in zip(sites, counts):
            tot = sum(v for k, v in n.items() if k != "indel")
            oth = ",".join(f"{k}:{v}" for k, v in n.most_common() if k not in ("TAG", "CAG"))
            w.write(f"{c}\t{r['pos1']}\t{s}\t{r['query']}\t{r['qpos']}\t{r['qres']}\t{r['win_ident']}\t"
                    f"{r['edge_dist']}\t{tot}\t{n['TAG']}\t{n['CAG']}\t{oth}\n")
    with open(out + ".summary.txt", "w") as w:
        for kind in ("chrom", "unplaced"):
            a = ident[(kind, "alen")]
            w.write(f"{kind}: {nreads[kind]} primary alignments mapq>=10, identity "
                    f"{ident[(kind, 'match')] / a if a else 0:.4f} over {a} aligned bases\n")
        cov_sites = [n for n in counts if sum(v for k, v in n.items() if k != "indel")]
        agg = Counter()
        for n in counts:
            agg.update(n)
        w.write(f"sites {len(sites)}; covered {len(cov_sites)}; read codons: {dict(agg.most_common(12))}\n")
        maj = Counter()
        for n in cov_sites:
            m = Counter({k: v for k, v in n.items() if k != "indel"}).most_common(1)[0][0]
            maj[m] += 1
        w.write(f"majority codon per covered site: {dict(maj.most_common())}\n")
    print(open(out + ".summary.txt").read(), end="")


if __name__ == "__main__":
    main(*sys.argv[1:5])
