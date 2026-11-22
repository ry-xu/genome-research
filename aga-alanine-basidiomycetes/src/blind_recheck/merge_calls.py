"""Merge Codetta (method 1) and ortholog tally (method 2) per codon per genome.
usage: merge_calls.py <codetta_vs_t4.tsv> <all.codon_table.tsv> <out.tsv>
Class per codon:
  t4_both      : Codetta calls table-4 aa, or method 2 lenient n>=5 and top residue = table-4 aa
  DEVIATION    : either method calls something else (listed for manual review)
  too_rare     : Codetta '?' AND method 2 lenient n_cons < 5
  t4_m2_only   : Codetta '?' but method 2 lenient n>=5 with table-4 aa on top and no flag
"""
import sys, csv, collections
cod = {(r['acc'], r['codon']): r for r in csv.DictReader(open(sys.argv[1]), delimiter='\t')}
m2 = {}
for r in csv.DictReader(open(sys.argv[2]), delimiter='\t'):
    m2[(r['acc'], r['codon'], r['def'])] = r
out = open(sys.argv[3], 'w')
out.write('acc\tcodon\tt4\tcodetta\tN\tn_in_genes\tstrict\tlenient\tm2_flag_lenient\tclass\n')
rare = collections.defaultdict(list)
for (acc, c), r in sorted(cod.items()):
    s, l = m2[(acc, c, 'strict')], m2[(acc, c, 'lenient')]
    def desc(x):
        return f"{x['n_cons']}:{x['breakdown'] or '-'}"
    t4 = r['t4']
    ln = int(l['n_cons'])
    ltop = l['breakdown'].split(':')[0] if l['breakdown'] else '-'
    if r['codetta'] not in ('?', t4) or l['flag'] in ('DEVIATION', 'partial') or s['flag'] in ('DEVIATION', 'partial'):
        cls = 'DEVIATION_or_flag'
    elif r['codetta'] == t4:
        cls = 't4_both' if (ln >= 5 and ltop == t4) else 't4_codetta'
    elif ln < 5:
        cls = 'too_rare'
    elif ltop == t4 or (t4 == '*'):
        cls = 't4_m2_only'
    else:
        cls = 'm2_unclear'
    out.write(f"{acc}\t{c}\t{t4}\t{r['codetta']}\t{r['N_aligned']}\t{l['n_in_genes']}\t{desc(s)}\t{desc(l)}\t{l['flag']}\t{cls}\n")
    if cls in ('too_rare', 'm2_unclear', 't4_m2_only'):
        rare[(acc, cls)].append(f"{c}(N={r['N_aligned']},n={ln},in={l['n_in_genes']})")
for (acc, cls), v in sorted(rare.items()):
    print(acc, cls, ' '.join(v))
