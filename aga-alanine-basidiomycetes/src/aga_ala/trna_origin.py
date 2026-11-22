"""Where does each genome's tRNA(UCU) come from? Identity to the same genome's other tRNAs.
usage: python trna_origin.py aragorn.tsv [acc_list]
For each tRNA with anticodon TCT, global-aligns it (match 1, mismatch -1, gap -2) to every other tRNA of its genome and
prints identity (matches / aligned columns) to the best tRNA-Ala, best other tRNA-Arg, and the best tRNA overall."""
import sys, csv
from collections import defaultdict
from Bio import Align

al = Align.PairwiseAligner(mode='global', match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1)


def ident(a, b):
    aln = al.align(a, b)[0]
    s1, s2 = aln[0], aln[1]
    cols = [(x, y) for x, y in zip(s1, s2) if not (x == '-' and y == '-')]
    return sum(x == y for x, y in cols) / len(cols)


def main(path, accs=None):
    keep = set(open(accs).read().split()) if accs else None
    by = defaultdict(list)
    for r in csv.DictReader(open(path), delimiter='\t'):
        if keep is None or r['acc'] in keep:
            by[r['acc']].append(r)
    print('acc\tpair3_70\tbest_ala\tbest_arg_other\tbest_any\tbest_any_aa_ac')
    for acc, rs in sorted(by.items()):
        for q in rs:
            if q['anticodon'].upper() != 'TCT':
                continue
            res = []
            for r in rs:
                if r is q:
                    continue
                res.append((ident(q['seq'], r['seq']), r['aa'], r['anticodon']))
            ala = max([x for x in res if x[1] == 'Ala'], default=(0, '', ''))
            arg = max([x for x in res if x[1] == 'Arg' and x[2].upper() != 'TCT'], default=(0, '', ''))
            top = max(res, default=(0, '', ''))
            print(f"{acc}\t{q['pair3_70']}\t{ala[0]:.2f}\t{arg[0]:.2f}\t{top[0]:.2f}\t{top[1]}:{top[2]}")


if __name__ == '__main__':
    main(*sys.argv[1:])
