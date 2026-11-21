"""Arm T of RULES.md: origin (D), own-genome nearest class (S1) and identity bases (S3) of CUN-reading UAG tRNAs.
usage: python origin.py trnas.tsv control|target
Identity: global alignment (match 1, mismatch -1, gap -2/-1, as trna_origin.py); the anticodon is set to N in every
sequence and N columns are left out of the identity (matches / non-N aligned columns)."""
import sys, csv
from collections import defaultdict, Counter
from statistics import mean
from Bio import Align

al = Align.PairwiseAligner(mode='global', match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1)
WICK = ['BK074538.1', 'BK074532.1', 'BK074542.1']
SPIN = ['BK074262.1']
SIS = {**{a: ['BK074528.1', 'BK074537.1', 'BK074540.1'] for a in WICK},
       'BK074262.1': ['BK074255.1', 'BK074256.1', 'BK074258.1', 'BK074260.1', 'BK074265.1'],
       'NC_005789.1': ['NC_006626.1', 'NC_018056.1', 'NC_025331.1', 'BK073244.1', 'HE983614.1']}
DIPO = None  # filled from accs.txt
CUN = ('TAG', 'TTAG')


def masked(r):
    s, o, n = r['seq'].upper(), int(r['ac_off']) - 1, int(r['ac_len'])
    return s[:o] + 'N' * n + s[o + n:]


def ident(a, b):
    x = al.align(a, b)[0]
    cols = [(p, q) for p, q in zip(x[0], x[1]) if p != 'N' and q != 'N' and not (p == '-' and q == '-')]
    return sum(p == q for p, q in cols) / len(cols)


def ac(r):
    return r['anticodon'].upper()


def feats(r):
    """Sprinzl-like positions: 1-7 (stem5), 66-72 (stem3), 73, anticodon loop 32,33,37,38."""
    s, o, n = r['seq'].upper(), int(r['ac_off']) - 1, int(r['ac_len'])
    f = {f'p{i + 1}': r['stem5'][i] for i in range(min(7, len(r['stem5'])))}
    f.update({f'p{66 + i}': r['stem3'][i] for i in range(min(7, len(r['stem3'])))})
    f['p73'] = r['N73']
    f.update({'p32': s[o - 2:o - 1], 'p33': s[o - 1:o], 'p37': s[o + n:o + n + 1], 'p38': s[o + n + 1:o + n + 2]})
    return f


def load(path):
    by = defaultdict(list)
    for r in csv.DictReader(open(path), delimiter='\t'):
        by[r['acc']].append(r)
    return by


def D(by, q, sisters, old_ac, new_ac):
    old = [ident(masked(q), masked(r)) for s in sisters for r in by[s] if ac(r) in old_ac]
    new = [ident(masked(q), masked(r)) for s in sisters for r in by[s] if ac(r) in new_ac]
    if not old or not new:
        return None, len(old), len(new)
    return mean(old) - mean(new), len(old), len(new)


def S1(by, acc, q):
    res = sorted(((ident(masked(q), masked(r)), r['aa'], ac(r)) for r in by[acc]
                  if r is not q and ac(r) not in CUN and r['aa'] != 'Leu'), reverse=True)
    return res[:3]


def specific_sets(by, accs):
    """(pos, base) carried by >=90% of a class and <20% of all other tRNAs; classes by anticodon, UAG tRNAs left out."""
    pool = [r for a in accs for r in by[a] if ac(r) not in CUN and len(ac(r)) == 3]
    cls = defaultdict(list)
    for r in pool:
        cls[ac(r)].append(feats(r))
    out = {}
    for c, fs in cls.items():
        if len(fs) < 0.8 * len(accs):
            continue
        others = [feats(r) for r in pool if ac(r) != c]
        keep = []
        for p in fs[0]:
            for b in 'ACGT':
                inc = sum(f.get(p) == b for f in fs) / len(fs)
                oth = sum(f.get(p) == b for f in others) / len(others)
                if inc >= 0.9 and oth < 0.2:
                    keep.append((p, b))
        out[c] = keep
    return out


def S3(by, q, sisters, sets):
    fq = feats(q)
    sis = [feats(r) for s in sisters for r in by[s] if ac(r) in CUN]
    rows = {}
    for c, keep in sets.items():
        carried = [(p, b) for p, b in keep if fq.get(p) == b]
        gained = [(p, b) for p, b in carried if sum(f.get(p) == b for f in sis) < len(sis) / 2]
        rows[c] = (len(keep), carried, gained)
    return rows


def calib(by, acc, sisters):
    out = []
    for x in by[acc]:
        a = ac(x)
        if a in CUN or a == 'TGG' or len(a) != 3:
            continue
        d, n1, n2 = D(by, x, sisters, (a,), ('TGG',))
        if d is not None and n1:
            out.append((round(d, 3), x['aa'] + '(' + a + ')'))
    return sorted(out)


def main(path, mode):
    by = load(path)
    if mode == 'control':
        q = [r for r in by['NC_005789.1'] if ac(r) in CUN]
        print('PC-T1 Ashbya UAG tRNAs:', len(q))
        for r in q:
            d, n1, n2 = D(by, r, SIS['NC_005789.1'], CUN, ('TGC',))
            print(f"  D(UAG vs Ala-TGC) = {d if d is None else round(d, 3)}  (n sis UAG {n1}, n sis Ala {n2}); "
                  f"pair3_70 {r['pair3_70']} N73 {r['N73']}")
            sis = [(s, x['pair3_70'], x['anticodon']) for s in SIS['NC_005789.1'] for x in by[s] if ac(x) in CUN]
            print('  sisters UAG 3:70 ->', sis)
            print('  PASS' if d is not None and d > 0.05 and r['pair3_70'] == 'G:T' and all(p != 'G:T' for _, p, _ in sis)
                  else '  FAIL')
            print('  calibration D_X (Ashbya, vs sisters TGG):', calib(by, 'NC_005789.1', SIS['NC_005789.1'])[:3], '...',
                  calib(by, 'NC_005789.1', SIS['NC_005789.1'])[-3:])
        q = [r for r in by['NC_001224.1'] if ac(r) in CUN]
        print('PC-T2 S. cerevisiae UAG tRNAs:', len(q))
        for r in q:
            top = S1(by, 'NC_001224.1', r)
            print('  S1 top3:', [(round(s, 3), a, c) for s, a, c in top], 'PASS' if top[0][1] == 'His' else 'FAIL')
        return
    accs = open('data/tree/accs.txt').read().split()
    sets = specific_sets(by, accs)
    print('Pro-specific (TGG) bases:', sets.get('TGG'))
    print('classes with a specific set:', {c: len(k) for c, k in sets.items()})
    for acc in WICK + SPIN:
        sisters = SIS[acc]
        for q in [r for r in by[acc] if ac(r) in CUN]:
            d, n1, n2 = D(by, q, sisters, CUN, ('TGG',))
            print(f'\n== {acc} UAG tRNA ({len(q["seq"])} nt, 3:70 {q["pair3_70"]}, N73 {q["N73"]})')
            print(f'  D(UAG vs Pro) = {round(d, 3)}  (n sis UAG {n1}, n sis Pro {n2})')
            print('  id to own Pro:', [round(ident(masked(q), masked(r)), 3) for r in by[acc] if ac(r) == 'TGG'])
            print('  S1 top3:', [(round(s, 3), a, c) for s, a, c in S1(by, acc, q)])
            cal = calib(by, acc, sisters)
            print(f'  calibration D_X n={len(cal)} min {cal[0]} median {cal[len(cal) // 2]} max {cal[-1]}')
            rows = S3(by, q, sisters, sets)
            pro = rows.get('TGG')
            print('  S3 Pro-specific carried', pro[1] if pro else None, 'gained vs sisters UAG', pro[2] if pro else None)
            dec = [(c, len(v[2])) for c, v in rows.items() if c != 'TGG' and v[2]]
            print('  S3 decoy gains (other classes):', dec)
    # relatives: UAG tRNA identity to own Pro, for comparison with targets
    print('\nall 49: UAG tRNA identity to own tRNA-Pro (masked)')
    for acc in accs:
        for q in [r for r in by[acc] if ac(r) in CUN]:
            p = [round(ident(masked(q), masked(r)), 3) for r in by[acc] if ac(r) == 'TGG']
            print(' ', acc, '*' if acc in WICK + SPIN else ' ', p, q['pair3_70'], q['N73'], len(q['seq']))


if __name__ == '__main__':
    main(*sys.argv[1:])
