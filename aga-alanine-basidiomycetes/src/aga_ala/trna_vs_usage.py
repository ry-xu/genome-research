"""Join per-genome AGA placement (codon_contrast.py output) with ARAGORN tRNA(UCU) acceptor stems.
usage: python trna_vs_usage.py contrast.tsv aragorn.tsv out.tsv
A genome reads AGA as Ala-like if AGA sits at >= 1 strict Ala column and 0 strict Arg columns; Arg-like if at >= 3 Arg
columns and 0 Ala columns. For each genome, reports its tRNA(UCU) count and how many have G3:U70 (ARAGORN 'G:T') and A73.
Prints the 2x2 of G3:U70 tRNA(UCU) presence by class."""
import sys, csv
from collections import defaultdict


def main(contrast, aragorn, out):
    use = {}
    with open(contrast) as f:
        next(f)
        for r in csv.DictReader(f, delimiter='\t'):
            if r['acc'] == 'pooled':
                continue
            use[r['acc']] = (r['genus'], int(r['codon@Y']), int(r['codon@STD']))
    tr = defaultdict(lambda: [0, 0, 0])  # UCU tRNAs, with G3:U70, with G3:U70 and A73
    with open(aragorn) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['anticodon'] != 'TCT':
                continue
            t = tr[r['acc']]
            t[0] += 1
            if r['pair3_70'] == 'G:T':
                t[1] += 1
                t[2] += r['N73'] == 'A'
    n = {('ala', True): 0, ('ala', False): 0, ('arg', True): 0, ('arg', False): 0}
    with open(out, 'w') as o:
        o.write('acc\tgenus\tAGA_at_Ala\tAGA_at_Arg\tclass\ttRNA_UCU\tG3U70\tG3U70_A73\n')
        for acc, (genus, a, r) in sorted(use.items(), key=lambda x: (-x[1][1], x[1][2])):
            cls = 'ala' if a >= 1 and r == 0 else 'arg' if r >= 3 and a == 0 else 'other'
            t = tr.get(acc, [0, 0, 0])
            if cls != 'other':
                n[(cls, t[1] > 0)] += 1
            o.write(f'{acc}\t{genus}\t{a}\t{r}\t{cls}\t{t[0]}\t{t[1]}\t{t[2]}\n')
    print(f"AGA only at Ala columns: {n[('ala', True)]} of {n[('ala', True)] + n[('ala', False)]} have a G3:U70 tRNA(UCU)")
    print(f"AGA at >=3 Arg, 0 Ala:   {n[('arg', True)]} of {n[('arg', True)] + n[('arg', False)]}")


if __name__ == '__main__':
    main(*sys.argv[1:])
