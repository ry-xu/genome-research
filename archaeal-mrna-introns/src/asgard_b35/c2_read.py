import numpy as np, re
L = np.load('c2_W31dna_loci.npy')
ref = -3.60  # L02 window, DNA MFE from c2.log CHECK line (same code)
ex = [4487605, 859933]
for l in open('/Users/ryan/lab/agents/altair/kills/culture_circles/aragorn_introns.txt'):
    m = re.search(r'\[(\d+),(\d+)\].*i\(', l)
    if m: ex += [int(m.group(1)), int(m.group(2))]
print('excluded anchors', len(ex))
def near(a): return any(abs(a - e) <= 100 + 31 for e in ex)
sel = np.array([r for r in L if r[1] <= ref and not near(r[0])])
share = sel[:, 3] / (sel[:, 3] + sel[:, 2])
n = len(sel); hi = int((share >= 0.10).sum())
print(f'loci DNA dG <= {ref}: {n}; with share >= 0.10: {hi} ({hi/n:.4%}); with >=1 removal read: {int((sel[:,3]>=1).sum())}; removal reads {int(sel[:,3].sum())} / through {int(sel[:,2].sum())}')
for r, s in sorted(zip(sel.tolist(), share), key=lambda x: -x[1])[:8]: print('  top', [round(x, 2) for x in r], round(s, 3))
w = 'structure skip is general here' if n >= 5 and hi / n >= 0.05 else ('not general' if n >= 50 and hi / n < 0.01 else 'unclear')
print('C2 word:', w)
# context: rank of L02 DNA dG among all loci, and stronger subsets
print('fraction of all loci at or below L02 DNA dG', round((L[:, 1] <= ref).mean(), 3))
for cut in (-5, -7, -9):
    m = np.array([r[1] <= cut and not near(r[0]) for r in L])
    print(f'  DNA dG<={cut}: loci {m.sum()}, share>=0.10 {int(((L[m,3]/(L[m,3]+L[m,2]))>=0.10).sum())}, removal reads {int(L[m,3].sum())} / through {int(L[m,2].sum())}')
