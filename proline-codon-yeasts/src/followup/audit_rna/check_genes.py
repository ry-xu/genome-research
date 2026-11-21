"""(a) rebuild each gene from genes.tsv rows, check codon_pos order, translate, align to S. cerevisiae proteins.
usage: check_genes.py ref.fa genes.tsv table sc_dir"""
import sys, csv, glob
from collections import OrderedDict
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
fa, genes, table, scdir = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
g = ''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper()
COMP = str.maketrans('ACGTN', 'TGCAN')
rows = OrderedDict()
bad = 0
for r in csv.DictReader(open(genes), delimiter='\t'):
    rows.setdefault(r['gene'], []).append(r)
sc = {}
for f in glob.glob(scdir + '/sc_*.fa'):
    sc[f.split('sc_')[1][:-3]] = ''.join(l.strip() for l in open(f) if not l.startswith('>'))
al = Align.PairwiseAligner(); al.mode = 'local'; al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score = -11; al.extend_gap_score = -1
for name, rs in rows.items():
    nt = ''
    for k, r in enumerate(rs):
        if int(r['codon_pos']) != k % 3 + 1: bad += 1
        b = g[int(r['pos1']) - 1]
        nt += b if r['strand'] == '1' else b.translate(COMP)
    p = str(Seq(nt[:len(nt)//3*3]).translate(table=table))
    stops = p[:-1].count('*')
    line = f'{name}\tlen_nt {len(nt)}\tstrand {rs[0]["strand"]}\tstart {nt[:3]}\tend {p[-1]}\tinternal_stops {stops}'
    for key, sp in sc.items():
        if key in name.lower() or (key == 'cob' and 'cob' in name.lower()):
            a = al.align(p.replace('*', 'X'), sp)[0]
            # reverse-strand control
            prc = str(Seq(nt[::-1].translate(COMP)[:len(nt)//3*3]).translate(table=table)).replace('*', 'X')
            arc = al.score(prc, sp)
            ident = sum(1 for x, y in zip(*[s for s in a]) if x == y and x != '-')
            alen = sum(1 for x, y in zip(*[s for s in a]) if x != '-' and y != '-')
            line += f'\tvs_Sc_{key}: score {a.score:.0f} (revcomp {arc:.0f}) id {ident}/{alen} q{a.coordinates[0][0]}-{a.coordinates[0][-1]} Sc{a.coordinates[1][0]}-{a.coordinates[1][-1]}'
    print(line)
print('codon_pos order mismatches', bad)
