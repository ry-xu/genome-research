"""Which amino acid's 'consensus-residue profile' best explains each codon's tally?
For codon X, tally T_X = relatives' consensus residues at conserved columns where the query has X.
Profile of amino acid a = pooled tallies of all table-4 codons of a except X (plus pseudocount 0.5/residue).
logL(X|a) = sum_r T_X[r] * log p_a(r). Reports best a, table-4 a, and log-likelihood difference.
usage: profile_fit.py <codon_table.tsv from codon_stats.py> <def strict|lenient> [codon ...]"""
import sys, csv, collections, math
from Bio.Data import CodonTable
t4 = dict(CodonTable.unambiguous_dna_by_id[4].forward_table)
for s in CodonTable.unambiguous_dna_by_id[4].stop_codons:
    t4[s] = '*'
AAS = 'ACDEFGHIKLMNPQRSTVWY'
f, dname = sys.argv[1], sys.argv[2]
only = set(sys.argv[3:])
tal = collections.defaultdict(dict)
for r in csv.DictReader(open(f), delimiter='\t'):
    if r['def'] != dname:
        continue
    c = collections.Counter()
    for kv in r['breakdown'].split():
        k, v = kv.split(':'); c[k] = int(v)
    tal[r['acc']][r['codon']] = c
print('acc\tcodon\tt4\tn\tbest\tlogL_best-logL_t4\tsecond\tlogL_best-logL_second')
for acc in sorted(tal):
    T = tal[acc]
    for X in sorted(T):
        if only and X not in only:
            continue
        n = sum(T[X].values())
        if n == 0:
            continue
        ll = {}
        for a in AAS:
            prof = collections.Counter()
            for c, aa in t4.items():
                if aa == a and c != X:
                    prof.update(T.get(c, {}))
            tot = sum(prof.values())
            if tot < 20:
                continue
            ll[a] = sum(v * math.log((prof.get(r, 0) + 0.5) / (tot + 0.5 * 20)) for r, v in T[X].items())
        ranked = sorted(ll, key=lambda a: -ll[a])
        A = t4[X]
        d4 = ll[ranked[0]] - ll[A] if A in ll else float('nan')
        print(f'{acc}\t{X}\t{A}\t{n}\t{ranked[0]}\t{d4:.1f}\t{ranked[1]}\t{ll[ranked[0]] - ll[ranked[1]]:.1f}')
