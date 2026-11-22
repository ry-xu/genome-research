"""Per human-position codons of one query in one genome from tblastn outfmt 6 (all HSPs on the best contig).
usage: codons.py GENOME HITS QUERYSUBSTR CONTIG > tsv (qpos, qres, sres, codon, contig_pos)"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stops_in_hsps import read_fasta, walk
genome, hits, q, contig = sys.argv[1:5]
seqs = read_fasta(genome)
out = {}
for line in open(hits):
    f = line.rstrip("\n").split("\t")
    if q not in f[0]:
        continue
    s = f[1].split("|")[1] if "|" in f[1] else f[1]
    if s != contig:
        continue
    h = dict(qstart=int(f[6]), sstart=int(f[8]), sframe=int(f[13]), qseq=f[15], sseq=f[16], bits=float(f[11]))
    for col, qp, a, b, cod, p in walk(h, seqs[s]):
        if qp is not None and qp not in out:
            out[qp] = (a, b, cod, p + 1)
for k in sorted(out):
    print(k, *out[k], sep="\t")
