# Task 5: what is behind T0 'non-skip' calls; power of p1/cgg counts
import sys, csv, math
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
def binom_cdf(k, n, p): return sum(math.comb(n, i) * p**i * (1-p)**(n-i) for i in range(k+1))
sk = [r for r in rows if r['call'] == '1']
tot_p1 = sum(int(r['p1']) for r in sk); tot_c = sum(int(r['cgg']) for r in sk)
print(f'skipper genomes n={len(sk)} pooled cgg/p1 = {tot_c}/{tot_p1} = {tot_c/tot_p1:.3f}; min p1 among skippers = {min(int(r["p1"]) for r in sk)}; median p1 = {sorted(int(r["p1"]) for r in sk)[len(sk)//2]}')
ps = tot_c / tot_p1
print('genome\tfamily\tgenus\tgc\tloci\tp1\tcgg\tshare\tcall\tP(cgg<=obs|skip rate)\tP(call0-like: share<0.3 | skipper)\tP(cgg>=obs | null 0.02)\tP(cgg>=obs | null 0.05)')
for r in rows:
    if r['call'] == '1': continue
    n, k = int(r['p1']), int(r['cgg'])
    lo = binom_cdf(k, n, ps)
    # prob a true skipper with this p1 would show share < 0.3
    lo30 = binom_cdf(math.ceil(0.3*n) - 1, n, ps) if n else 1.0
    hi2 = 1 - binom_cdf(k-1, n, 0.02) if k else 1.0
    hi5 = 1 - binom_cdf(k-1, n, 0.05) if k else 1.0
    T0 = ' <-- T0 genus' if r['genus'] in ('CAG-492', 'CAG-465') else ''
    print(f"{r['gca']}\t{r['family']}\t{r['genus']}\t{r['gc']}\t{r['loci']}\t{n}\t{k}\t{r['share']}\t{r['call']}\t{lo:.3f}\t{lo30:.3f}\t{hi2:.4f}\t{hi5:.4f}{T0}")
# p1 rate per locus: skippers vs T0
import statistics
rate = [int(r['p1'])/int(r['loci']) for r in sk]
print('skipper p1/loci: min %.4f median %.4f' % (min(rate), statistics.median(rate)))
for r in rows:
    if r['genus'] in ('CAG-492','CAG-465'):
        print(r['genus'], 'p1/loci %.4f' % (int(r['p1'])/int(r['loci'])), 'gc', r['gc'])
# skippers with similar GC (<0.265)
for r in sk:
    if float(r['gc']) < 0.265: print('  low-GC skipper', r['genus'], r['gc'], r['p1'], r['loci'])
