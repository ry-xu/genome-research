# C1 (RULES.md): duplicate-family concordance at L02, place 0
import csv, collections, math
rows=[r for r in csv.DictReader(open('/Users/ryan/lab/work/archaeal-introns/culture/out/rna_mol_templates.tsv'),delimiter='\t') if r['site']=='L02' and r['place']=='0']
cls=collections.Counter(r['class'] for r in rows); print('classes',dict(cls))
fam=collections.defaultdict(lambda: collections.Counter())
for r in rows: fam[(r['runs'],r['start'],r['end'])][r['class']]+=int(r['fragments'])
tot=sum(int(r['fragments']) for r in rows); rem=sum(int(r['fragments']) for r in rows if r['class']!='through')
p=rem/tot; print('reads',tot,'removal reads',rem,'p',round(p,4))
E=0; obs=0; n2=0
for k,c in fam.items():
    n=sum(c.values())
    if n<2: continue
    n2+=1; E+=1-p**n-(1-p)**n
    if 'through' in c and any(x!='through' for x in c): obs+=1; print('mixed',k,dict(c))
P=sum(math.exp(-E)*E**i/math.factorial(i) for i in range(obs+1))
print(f'families {len(fam)}; with n>=2 {n2}; expected mixed {E:.2f}; observed {obs}; Poisson P(X<=obs) {P:.3g}')
w = 'NOT TESTABLE' if E<5 else ('PCR-stage skip rejected' if obs<=E/4 and P<0.01 else ('not rejected' if obs>=E/2 else 'unclear'))
print('C1 word:', w)
