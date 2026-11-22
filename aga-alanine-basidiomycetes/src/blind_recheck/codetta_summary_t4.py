"""Summarize Codetta outputs against NCBI table 4.
usage: codetta_summary_t4.py <out_tsv> <codetta .genetic_code.out files...>"""
import sys, os
from Bio.Data import CodonTable
t4 = dict(CodonTable.unambiguous_dna_by_id[4].forward_table)
for s in CodonTable.unambiguous_dna_by_id[4].stop_codons:
    t4[s] = '*'
AAS = 'ACDEFGHIKLMNPQRSTVWY?'
out = open(sys.argv[1], 'w')
out.write('acc\tcodon\tt4\tcodetta\tN_aligned\tN_used\tlogP_t4aa\tlogP_call\tstatus\n')
for f in sys.argv[2:]:
    acc = os.path.basename(f).split('.fna')[0]
    lines = open(f).read().splitlines()
    calls, logp = {}, {}
    sec = None
    for ln in lines:
        if ln.startswith('# codon   inference'):
            sec = 'calls'; continue
        if ln.startswith('# codon      logP'):
            sec = 'logp'; continue
        if ln.startswith('#'):
            if sec and ln.strip() == '#': sec = None
            continue
        p = ln.split()
        if sec == 'calls' and len(p) >= 6:
            calls[p[0]] = (p[1], int(p[4]), int(p[5]))
        elif sec == 'logp' and len(p) == 22:
            logp[p[0]] = dict(zip(AAS, map(float, p[1:])))
    final = lines[-1].strip()
    devs, unc = [], []
    for c, (inf, na, nu) in calls.items():
        a = t4[c]
        lp = logp.get(c, {})
        lpa = lp.get(a, float('nan')) if a != '*' else float('nan')
        lpc = lp.get(inf, float('nan'))
        if inf == '?':
            st = 'uncalled'
            unc.append(f'{c}({na})')
        elif inf == a:
            st = 'agrees'
        else:
            st = 'DEVIATES'
            devs.append(f'{c}:{a}->{inf} N={na}')
        out.write(f'{acc}\t{c}\t{a}\t{inf}\t{na}\t{nu}\t{lpa:.2f}\t{lpc:.2f}\t{st}\n')
    print(f'{acc}\tcode={final}\n  deviations: {"; ".join(devs) or "none"}\n  uncalled(N): {" ".join(unc)}')
