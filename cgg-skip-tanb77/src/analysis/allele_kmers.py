"""Read support for the extra C at +1 CGG breaks: the assembled 25-mer (with the C) against the same 25-mer without it.
usage: allele_kmers.py sites <sites.tsv> <fasta> <out_kmers.tsv>   -> one row per (site, allele, kmer)
       allele_kmers.py count <kmers.tsv> <reads.txt> > counts.tsv   reads.txt: read sequences, one per line
sites.tsv: reads_check-style rows (label, qseqid, strand, sseqid, shift, ctx; ctx[15] is the break codon's first
base). Alleles: with = 25-mer centred on that base; without = the same span with the base removed; left/right =
25-mers 60 nt up- and downstream (depth). A real C gives with ~ flank depth and without ~ 0; an assembly insertion
gives the reverse.
Pipe the reads through `grep -F -f` with every kmer and its reverse complement first: a pure-Python scan of a whole
run is too slow."""
import csv
import gzip
import sys

K, FL = 25, 60
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    out, name = {}, None
    for line in op(path, "rt"):
        if line[0] == ">":
            name = line[1:].split()[0]
            out[name] = []
        else:
            out[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}


def sites(spath, fpath, opath):
    seqs = read_fasta(fpath)
    w = csv.writer(open(opath, "w"), delimiter="\t", lineterminator="\n")
    w.writerow(["site", "label", "sseqid", "allele", "kmer"])
    for i, r in enumerate(csv.DictReader(open(spath), delimiter="\t")):
        s = seqs[r["qseqid"]]
        s = rc(s) if r["strand"] == "-" else s
        c = r["ctx"].upper()
        p = s.find(c)
        if p < 0 or s.find(c, p + 1) >= 0:
            continue
        p += 15
        h = K // 2
        al = {"with": s[p - h:p + h + 1], "without": s[p - h:p] + s[p + 1:p + h + 2],
              "left": s[p - FL - h:p - FL + h + 1], "right": s[p + FL - h:p + FL + h + 1]}
        for a, k in al.items():
            w.writerow([i, r["label"], r["sseqid"], a, k])


def count(kpath, rpath):
    rows = list(csv.DictReader(open(kpath), delimiter="\t"))
    want = {}
    for r in rows:
        want[r["kmer"]] = 0
        want[rc(r["kmer"])] = 0
    for line in open(rpath):
        x = line.strip()
        for j in range(len(x) - K + 1):
            k = x[j:j + K]
            if k in want:
                want[k] += 1
    by = {}
    for r in rows:
        by.setdefault((r["site"], r["label"], r["sseqid"]), {})[r["allele"]] = want[r["kmer"]] + want[rc(r["kmer"])]
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["site", "label", "sseqid", "with", "without", "left", "right", "with_share"])
    for (i, lab, sid), c in by.items():
        n = c["with"] + c["without"]
        w.writerow([i, lab, sid, c["with"], c["without"], c.get("left", ""), c.get("right", ""),
                    f"{c['with'] / n:.3f}" if n else ""])


if __name__ == "__main__":
    if sys.argv[1] == "sites":
        sites(*sys.argv[2:5])
    else:
        count(*sys.argv[2:4])
