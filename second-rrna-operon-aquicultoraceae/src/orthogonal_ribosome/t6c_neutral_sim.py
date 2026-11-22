"""T6-C neutral null: same statistics on A. secundus div with 15% random substitutions (sequences from src/t6c_neutral_gen.py).
usage: uv run python src/t6c_neutral_sim.py out/t6/tmp/neutral*/all.sto"""
import collections, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from t5_count import binom_sf

OK = {('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G')}
aln, ss = collections.defaultdict(str), ''
for l in open(sys.argv[1]):
    if l.startswith('#=GC SS_cons'): ss += l.split()[2]
    elif l.strip() and not l.startswith(('#', '//')): n, s = l.split(); aln[n] += s.upper().replace('T', 'U')
st, pairs = {o: [] for o in '<([{'}, []
for i, c in enumerate(ss):
    if c in st: st[c].append(i)
    elif c in '>)]}': pairs.append((st['<([{'['>)]}'.index(c)]].pop(), i))
GEN = {'GCA_036504445.1': 'DASXVE01', 'GCA_902810715.1': 'CADDYS01', 'GCA_027392615.1': 'CADDYS01',
       'GCA_014859275.1': 'Kmv33', 'GCA_035643475.1': 'DASQQM01'}


def comp(a, b):
    A, B = aln[a], aln[b]
    ok = [(i, j) for i, j in pairs if all(x[k] in 'ACGU' for x in (A, B) for k in (i, j))]
    side = [A[k] != B[k] for i, j in ok for k in (i, j)]; f = sum(side) / len(side)
    one = [(i, j) for i, j in ok if A[i] != B[i] or A[j] != B[j]]
    both = [(i, j) for i, j in one if A[i] != B[i] and A[j] != B[j]]
    exp = f ** 2 / (1 - (1 - f) ** 2); p = binom_sf(len(both), len(one), exp)
    kept = sum((B[i], B[j]) in OK for i, j in one) / len(one)
    ratio = len(both) / len(one) / exp
    print(f'{a} vs {b}: pairs {len(ok)}, f {f:.3f}, changed pairs {len(one)}, both {len(both)} (x{ratio:.2f} independence, p {p:.2g}), '
          f'changed pairs still WC/GU {kept:.2f}; all pairs WC/GU in {b} {sum((B[i], B[j]) in OK for i, j in ok) / len(ok):.2f}')
    return kept, ratio, p


print('control (normal copies):')
ctrl = [comp('sec_norm', n) for n in aln if n.endswith('_norm') and n != 'sec_norm']
print('divergent copies:')
res = collections.defaultdict(list)
for n in aln:
    if n.endswith('_div') and n != 'sec_div': res[n].append(comp('sec_div', n))
mt = sum(all(k >= 0.8 and r >= 3 and p <= 0.01 for k, r, p in v) for v in res.values())
dc = sum(all(k < 0.6 for k, r, p in v) for v in res.values())
v = 'n/a (null, no control)' if not ctrl else 'UNSETTLED (control)' if min(k for k, r, p in ctrl) < 0.8 else 'MAINTAINED' if mt >= 3 else 'DECAYING' if dc >= 3 else 'UNSETTLED'
print(f'genera maintained {mt}/{len(res)}, decaying {dc}/{len(res)}\nVERDICT (T6-C): {v}')
