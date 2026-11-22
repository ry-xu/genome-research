"""Sanity check: every codon row in <acc>.columns.tsv matches the genome at its stated coordinate/strand.
usage: verify_coords.py <genomes_dir> <columns.tsv...>"""
import sys, csv, os
from Bio import SeqIO
from Bio.Seq import Seq
gd = sys.argv[1]
for f in sys.argv[2:]:
    rows = list(csv.DictReader(open(f), delimiter='\t'))
    acc = rows[0]['acc']
    G = str(next(SeqIO.parse(os.path.join(gd, acc + '.fna'), 'fasta')).seq).upper()
    bad = 0
    for r in rows:
        x = int(r['coord'])
        c = G[x - 1:x + 2] if r['strand'] == '1' else str(Seq(G[x - 3:x]).reverse_complement())
        bad += c != r['codon']
    print(acc, 'rows', len(rows), 'mismatches', bad)
