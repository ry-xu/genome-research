"""T6 (RULES_T6.md): divergent 16S RNA-bridge residues vs the normal 23S; PC-T6a intra-16S covariation control.
usage: uv run python src/t6_bridges.py pc     -> PC-T6a only
       uv run python src/t6_bridges.py t6     -> T6 (run only after PC-T6a passes)
"""
import collections, os, subprocess, sys
sys.path.insert(0, os.path.dirname(__file__))
from t5_count import binom_sf

O = 'out/t6'; T = f'{O}/tmp'
CM = 'data/external/rfam/RF00177.cm'
CMALIGN = 'cmalign'
PAIRS_OK = {('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G')}


def fasta(p):
    d, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: ''.join(v).upper().replace('T', 'U') for k, v in d.items()}


def binom_lower(x, n, p):
    return 1 - binom_sf(x + 1, n, p)


def pc():
    sec = fasta('out/t5/sec16.fa')
    f = f'{T}/sec16.fa'; open(f, 'w').write(''.join(f'>{k}\n{v}\n' for k, v in sec.items()))
    sto = f'{T}/sec16.sto'
    subprocess.run([CMALIGN, '--cpu', '1', '-o', sto, CM, f], check=True, stdout=subprocess.DEVNULL)
    aln, ss = collections.defaultdict(str), ''
    for l in open(sto):
        if l.startswith('#=GC SS_cons'): ss += l.split()[2]
        elif l.strip() and not l.startswith('#') and not l.startswith('//'):
            n, s = l.split(); aln[n] += s
    st, pairs = {'<': [], '(': [], '[': [], '{': []}, []
    close = {'>': '<', ')': '(', ']': '[', '}': '{'}
    for i, c in enumerate(ss):
        if c in st: st[c].append(i)
        elif c in close: pairs.append((st[close[c]].pop(), i))
    a, b = aln['norm'].upper().replace('T', 'U'), aln['div'].upper().replace('T', 'U')
    res = lambda s, i: s[i] if s[i] in 'ACGU' else None
    ok = [(i, j) for i, j in pairs if all(res(s, k) for s in (a, b) for k in (i, j))]
    side = [a[k] != b[k] for i, j in ok for k in (i, j)]
    fside = sum(side) / len(side)
    one = [(i, j) for i, j in ok if a[i] != b[i] or a[j] != b[j]]
    both = [(i, j) for i, j in one if a[i] != b[i] and a[j] != b[j]]
    exp = fside ** 2 / (1 - (1 - fside) ** 2)
    p = binom_sf(len(both), len(one), exp)
    pairing = sum((b[i], b[j]) in PAIRS_OK for i, j in both)
    norm_pairing = sum((a[i], a[j]) in PAIRS_OK for i, j in ok) / len(ok)
    share = len(both) / len(one) if one else 0
    passed = share >= 3 * exp and p <= 0.01 and both and pairing / len(both) >= 0.8
    print(f'PC-T6a: SS_cons pairs {len(pairs)}, aligned in both {len(ok)}; per-side change share f {fside:.3f}')
    print(f'  pairs with >= 1 side changed {len(one)}, both changed {len(both)} (share {share:.3f}; independence {exp:.3f}; '
          f'ratio {share / exp:.2f}; binomial p {p:.2e})')
    print(f'  both-changed pairs still WC/GU in div: {pairing}/{len(both)} ({pairing / max(1, len(both)):.2f}); '
          f'norm WC/GU at aligned pairs {norm_pairing:.2f}')
    print('PC-T6a', 'PASS' if passed else 'FAIL')


def align(seqs, tag):
    f = f'{T}/{tag}.in.fa'; open(f, 'w').write(''.join(f'>{k}\n{v}\n' for k, v in seqs.items()))
    out = subprocess.run(['mafft', '--localpair', '--maxiterate', '1000', '--thread', '1', '--quiet', f],
                         capture_output=True, text=True, check=True).stdout
    open(f'{T}/{tag}.aln.fa', 'w').write(out)
    al = {}
    for blk in out.split('>')[1:]:
        h, *s = blk.split('\n'); al[h.split()[0]] = ''.join(s).upper().replace('T', 'U')
    return al


def colmap(al, num):
    """E. coli number -> (ecoli, norm, div) bases or None per side; and per-column iteration."""
    m, k = {}, 0
    for c in range(len(al['ecoli'])):
        if al['ecoli'][c] != '-':
            n = num[k]; k += 1
            m[n] = tuple(None if al[s][c] == '-' else al[s][c] for s in ('ecoli', 'norm', 'div'))
    assert k == len(num)
    return m


def t6():
    con = [l.rstrip('\n').split('\t') for l in open(f'{O}/contacts_16s23s.tsv')][1:]
    pairs = [(int(r[0]), int(r[2])) for r in con]
    B16 = sorted({i for i, _ in pairs}); P23 = sorted({j for _, j in pairs})
    n16 = [(int(l.split('\t')[0]), l.split('\t')[1]) for l in open(f'{O}/ecoli16S_num.tsv')]
    n23 = [(int(l.split('\t')[0]), l.split('\t')[1]) for l in open(f'{O}/ecoli23S_num.tsv')]
    s16 = fasta('out/t5/sec16.fa'); s23 = fasta(f'{T}/sec23_raw.fa')  # 23S loci: plus strand vs E. coli (blastn, all HSPs plus)
    al16 = align({'ecoli': ''.join(b for _, b in n16), **s16}, 't6_16s')
    al23 = align({'ecoli': ''.join(b for _, b in n23), **s23}, 't6_23s')
    m16 = colmap(al16, [n for n, _ in n16]); m23 = colmap(al23, [n for n, _ in n23])
    ass = lambda t: t is not None and all(t)
    chg = lambda t: t[1] != t[2]
    def bg(m, excl):
        cols = [t for n, t in m.items() if n not in excl and ass(t)]
        return sum(map(chg, cols)) / len(cols), len(cols)
    f16, n16c = bg(m16, set(B16)); f23, n23c = bg(m23, set(P23))
    a16 = [i for i in B16 if ass(m16[i])]; c16 = [i for i in a16 if chg(m16[i])]
    print(f'B16 assessable {len(a16)}/{len(B16)}; changed {len(c16)}; background f16 {f16:.3f} over {n16c} columns; '
          f'f23 {f23:.3f} over {n23c} columns; P23 assessable {sum(ass(m23[j]) for j in P23)}/{len(P23)}')
    print('E. coli 16S identity to norm / div over aligned E. coli columns:',
          *(f"{sum(t[0] == t[k] for t in m16.values() if t[k]) / sum(1 for t in m16.values() if t[k]):.3f}" for k in (1, 2)))
    # co-adaptation
    co = [(i, j) for i, j in pairs if i in c16 and ass(m23[j])]
    co_both = [(i, j) for i, j in co if chg(m23[j])]
    p_co = binom_sf(len(co_both), len(co), f23) if co else None
    # bridge changes vs f16, both tails (reported)
    p_hi = binom_sf(len(c16), len(a16), f16); p_lo = binom_lower(len(c16), len(a16), f16)
    print(f'bridge changes {len(c16)}/{len(a16)} vs expected {f16 * len(a16):.1f}: P(>=) {p_hi:.3g}, P(<=) {p_lo:.3g}')
    print(f'co-adaptation: contact pairs with a changed B16 residue and assessable 23S partner {len(co)}; partner changed '
          f'{len(co_both)} (share {len(co_both) / max(1, len(co)):.2f}; f23 {f23:.3f}; binomial p {p_co})')
    if len(a16) < 24: v = 'UNASSESSABLE'
    elif len(c16) <= 1: v = 'COMPATIBLE'
    elif co and len(co_both) / len(co) >= 0.5 and p_co <= 0.01: v = 'CO-ADAPTED TO THE DIVERGENT 23S'
    else: v = 'CHANGED, PARTNERS NOT'
    print('\nper-residue (E. coli no.: ecoli/norm/div; partners E. coli no. ecoli/norm/div)')
    part = collections.defaultdict(list)
    for i, j in pairs: part[i].append(j)
    for i in B16:
        t = m16[i]; fmt = lambda t: '/'.join(x or '-' for x in t) if t else 'n/a'
        flag = '' if not ass(t) else ('  CHANGED' if chg(t) else '')
        print(f'  16S {i}: {fmt(t)}{flag}   ->', ', '.join(f'{j} {fmt(m23[j])}{"*" if ass(m23[j]) and chg(m23[j]) else ""}' for j in part[i]))
    rev = [j for j in P23 if ass(m23[j]) and chg(m23[j])]
    print('reverse view: changed P23 residues', rev, '; their 16S partners changed:',
          {j: [i for i, jj in pairs if jj == j and i in c16] for j in rev})
    # div 16S insertions: runs of >= 20 div bases in columns gapped in E. coli
    ins, run, last = [], 0, None
    k = 0; nums = [n for n, _ in n16]
    for c in range(len(al16['ecoli'])):
        if al16['ecoli'][c] != '-': last = nums[k]; k += 1
        if al16['ecoli'][c] == '-' and al16['div'][c] != '-': run += 1
        elif run:
            if run >= 20: ins.append((last, run))
            run = 0
    print('div 16S insertions >= 20 nt (after E. coli position, length):', ins)
    print('\nVERDICT (T6):', v)


if __name__ == '__main__':
    {'pc': pc, 't6': t6}[sys.argv[1]]()
