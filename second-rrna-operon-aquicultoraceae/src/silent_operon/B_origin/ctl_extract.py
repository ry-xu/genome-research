"""Extract control 16S/23S genes (E. coli K-12, V. vulnificus MO6-24/O, T. bispora DSM 43833) from nhmmer tables.
Writes out/ctl_lsu.fna, out/ctl_ssu.fna; ids  <org>_<kind>_<contig>_<start><strand>"""
from Bio import SeqIO
seqs = SeqIO.to_dict(SeqIO.parse("data/ctl/ctl.fna", "fasta"))
for kind, tbl, out in [("23S", "out/ctl_lsu.tbl", "out/ctl_lsu.fna"), ("16S", "out/ctl_ssu.tbl", "out/ctl_ssu.fna")]:
    with open(out, "w") as fh:
        for line in open(tbl):
            if line.startswith("#"):
                continue
            f = line.split()
            a, b, st = int(f[6]), int(f[7]), f[11]
            lo, hi = min(a, b), max(a, b)
            s = seqs[f[0]].seq[lo - 1:hi]
            if st == "-":
                s = s.reverse_complement()
            org, ctg = f[0].split("|")
            fh.write(f">{org}_{kind}_{ctg}_{lo}{st}\n{s}\n")
