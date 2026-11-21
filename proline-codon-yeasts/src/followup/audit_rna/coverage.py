"""(c) per-gene coverage with the caller's own pileup; exon vs intron-ORF part (Sc-aligned span for cox1/cob/cox2/cox3).
usage: coverage.py ref.fa genes.tsv bam sc_dir rdd.tsv"""
import sys, csv, glob, statistics
sys.path.insert(0, 'src/followup')  # run from the folder root
import rdd
from collections import OrderedDict
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
fa, genes, bam, scdir, rddtsv = sys.argv[1:6]
g = rdd.fasta1(fa); COMP = str.maketrans('ACGTN', 'TGCAN')
rows = OrderedDict()
for r in csv.DictReader(open(genes), delimiter='\t'):
    rows.setdefault(r['gene'], []).append((int(r['pos1']), int(r['strand']), int(r['codon_pos'])))
sc = {f.split('sc_')[1][:-3]: ''.join(l.strip() for l in open(f) if not l.startswith('>')) for f in glob.glob(scdir + '/sc_*.fa')}
al = Align.PairwiseAligner(); al.mode = 'local'; al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score = -11; al.extend_gap_score = -1
pr = rdd.pileup(fa, bam)
dep = {p: sum(c.values()) for p, c in pr.items()}
# positions on both strands are dropped by the caller; same here
strands = {}
for name, rs in rows.items():
    for p, s, _ in rs: strands.setdefault(p, set()).add(s)
calls = {int(r['pos']) for r in csv.DictReader(open(rddtsv), delimiter='\t')}
tot = tot20 = 0; ex_tot = ex20 = 0; intr20 = 0
print('gene\tnt\tdepth>=20\tfrac\tmedian_depth\texon_nt(Sc-aligned)\texon_depth>=20\tRDD_in_exon\tRDD_outside')
for name, rs in rows.items():
    nt = ''.join(g[p - 1] if s == 1 else g[p - 1].translate(COMP) for p, s, _ in rs)
    exon = None
    for key, sp in sc.items():
        if key in name.lower():
            prot = str(Seq(nt[:len(nt)//3*3]).translate(table=4)).replace('*', 'X')
            a = al.align(prot, sp)[0]
            exon = set(p for p, _, _ in rs[3 * a.coordinates[0][0]: 3 * a.coordinates[0][-1]])
    ds = [dep.get(p, 0) for p, _, _ in rs]
    n20 = sum(d >= 20 for d in ds); tot += len(rs); tot20 += n20
    e20 = sum(1 for p, _, _ in rs if exon and p in exon and dep.get(p, 0) >= 20)
    re = sum(1 for p, _, _ in rs if p in calls and exon and p in exon); ro = sum(1 for p, _, _ in rs if p in calls) - re
    if exon: ex_tot += len(exon); ex20 += e20; intr20 += n20 - e20
    print(f'{name}\t{len(rs)}\t{n20}\t{n20/len(rs):.2f}\t{statistics.median(ds)}\t{len(exon) if exon else "-"}\t{e20 if exon else "-"}\t{re if exon else "-"}\t{ro}')
print(f'TOTAL coding nt {tot}; depth>=20 {tot20} ({tot20/tot:.3f}); cox1/cob/cox2/cox3 Sc-aligned nt {ex_tot}, of which depth>=20 {ex20}; '
      f'depth>=20 in the non-Sc-aligned parts of those ORFs (intron ORF / UTR) {intr20}')
