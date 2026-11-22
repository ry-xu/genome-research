"""Per-codon tally of relatives' consensus residues at conserved columns, with Fisher tests.
usage: codon_stats.py <out_prefix> <columns.tsv files...>
Conserved column: >= MINRES relatives with a residue (of ~10-11) and
  strict: all of them identical; lenient: top residue >= 80%.
For codon X (table-4 aa A; synonyms S = other table-4 codons of A; for stop codons S = all sense codons):
  p_dev  : one-sided Fisher, X's columns carry residue Y more often than S's columns (Y = most common non-A residue)
  p_lossA: one-sided Fisher, X's columns carry A less often than S's columns
"""
import sys, collections, csv, math
from scipy.stats import fisher_exact
from Bio.Data import CodonTable

t4 = dict(CodonTable.unambiguous_dna_by_id[4].forward_table)
for s in CodonTable.unambiguous_dna_by_id[4].stop_codons:
    t4[s] = '*'
CODONS = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
MINRES = 6
prefix = sys.argv[1]
NTEST = 64 * 8
summary = open(prefix + '.codon_table.tsv', 'w')
summary.write('acc\tdef\tcodon\tt4\tn_in_genes\tn_cons\tn_t4aa\ttopY\tnY\tsyn_n\tsyn_t4aa\tsyn_Y\tp_dev\tp_lossA\tflag\tbreakdown\n')
usage = open(prefix + '.usage_at_aa_columns.tsv', 'w')
usage.write('acc\tdef\taa\tcodon\tcount_at_cons_cols_of_aa\ttotal_cons_cols_of_aa\n')

for f in sys.argv[2:]:
    rows = list(csv.DictReader(open(f), delimiter='\t'))
    acc = rows[0]['acc']
    n_in = collections.Counter(r['codon'] for r in rows)
    for dname, fn in [('strict', lambda r: int(r['nres']) >= MINRES and int(r['ntop']) == int(r['nres'])),
                      ('lenient', lambda r: int(r['nres']) >= MINRES and int(r['ntop']) >= 0.8 * int(r['nres']))]:
        # exclude initiator columns (>=2 relatives start their protein here): start codons are not sense-codon evidence
        cons = [r for r in rows if fn(r) and r['top'] != '-' and int(r['nstart']) < 2]
        tally = collections.defaultdict(collections.Counter)  # codon -> residue counts
        for r in cons:
            tally[r['codon']][r['top']] += 1
        bycol = collections.defaultdict(collections.Counter)  # consensus aa -> query codon counts
        for r in cons:
            bycol[r['top']][r['codon']] += 1
        for aa in sorted(bycol):
            tot = sum(bycol[aa].values())
            for c in CODONS:
                if bycol[aa][c]:
                    usage.write(f'{acc}\t{dname}\t{aa}\t{c}\t{bycol[aa][c]}\t{tot}\n')
        flagged = []
        for c in CODONS:
            A = t4[c]
            syn = [s for s in CODONS if t4[s] == A and s != c] if A != '*' else [s for s in CODONS if t4[s] != '*']
            T = tally[c]
            n = sum(T.values())
            nA = T.get(A, 0)
            nonA = [(k, v) for k, v in T.most_common() if k != A]
            Y, nY = nonA[0] if nonA else ('-', 0)
            S = collections.Counter()
            for s in syn:
                S.update(tally[s])
            sn = sum(S.values()); sA = S.get(A, 0); sY = S.get(Y, 0)
            p_dev = fisher_exact([[nY, n - nY], [sY, sn - sY]], alternative='greater')[1] if n and Y != '-' else 1.0
            p_loss = fisher_exact([[nA, n - nA], [sA, sn - sA]], alternative='less')[1] if n and A != '*' else float('nan')
            flag = ''
            if p_dev < 0.05 / NTEST and nY > nA:
                flag = 'DEVIATION'
            elif p_dev < 0.05 / NTEST:
                flag = 'partial'
            elif n < 5:
                flag = 'too_rare'
            br = ' '.join(f'{k}:{v}' for k, v in T.most_common())
            summary.write(f'{acc}\t{dname}\t{c}\t{A}\t{n_in[c]}\t{n}\t{nA}\t{Y}\t{nY}\t{sn}\t{sA}\t{sY}\t{p_dev:.2e}\t{p_loss:.2e}\t{flag}\t{br}\n')
            if flag in ('DEVIATION', 'partial'):
                flagged.append(f'{c}({A}): {Y} {nY}/{n} [A {nA}/{n}; syn {Y} {sY}/{sn}, syn {A} {sA}/{sn}] p_dev={p_dev:.1e} p_lossA={p_loss:.1e} {flag}')
        print(acc, dname, 'conserved codon columns', len(cons))
        for x in flagged:
            print('   ', x)
