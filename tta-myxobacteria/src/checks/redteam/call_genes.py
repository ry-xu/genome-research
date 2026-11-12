"""Red-team gene calls: Pyrodigal single mode, table 11 (TTA read as Leu), one genome per call.

usage: call_genes.py <genome.fna.gz> <label> <outdir>
writes <outdir>/<label>.genes.tsv  (gid, contig, begin, end, strand, partial_begin, partial_end, nt incl. stop)
       <outdir>/<label>.faa        (complete-or-partial proteins >= 30 aa, stop removed; id = label|gid)
       <outdir>/<label>.contigs.tsv (contig, length, gc)
"""
import gzip
import sys

import pyrodigal

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(fn):
    op = gzip.open if fn.endswith(".gz") else open
    name, buf = None, []
    with op(fn, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf).upper()


def main(fn, label, outdir):
    contigs = list(read_fasta(fn))
    gf = pyrodigal.GeneFinder()
    gf.train(*[s.encode() for _, s in contigs], translation_table=11)
    with open(f"{outdir}/{label}.genes.tsv", "w") as g, open(f"{outdir}/{label}.faa", "w") as p, \
            open(f"{outdir}/{label}.contigs.tsv", "w") as c:
        g.write("gid\tcontig\tbegin\tend\tstrand\tpb\tpe\tnt\n")
        for n, s in contigs:
            gc = (s.count("G") + s.count("C")) / max(1, len(s) - s.count("N"))
            c.write(f"{n}\t{len(s)}\t{gc:.4f}\n")
            for i, gene in enumerate(gf.find_genes(s.encode())):
                nt = s[gene.begin - 1:gene.end]
                if gene.strand == -1:
                    nt = rc(nt)
                gid = f"{n}_{gene.begin}_{gene.end}_{gene.strand}"
                g.write(f"{gid}\t{n}\t{gene.begin}\t{gene.end}\t{gene.strand}\t{int(gene.partial_begin)}\t"
                        f"{int(gene.partial_end)}\t{nt}\n")
                prot = gene.translate(translation_table=11).rstrip("*")
                if len(prot) >= 30:
                    p.write(f">{label}|{gid}\n{prot}\n")


if __name__ == "__main__":
    main(*sys.argv[1:4])
