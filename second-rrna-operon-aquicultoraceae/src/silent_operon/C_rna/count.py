"""Count divergent vs normal 23S and 16S read pairs per species group (RULES.md), plus the rRNA-share control (C1).

usage: python3 -I src/count.py PAF CLASSES.tsv [TOTAL_PAIRS]
PAF: minimap2 -c lines (fastq-dump --split-spot style: both mates share one name; minimap2 may append /1 /2), already
filtered to >= 60 aligned columns and identity >= 0.85. Counting uses identity >= 0.98 only (the earlier method);
the 0.85 lines feed the rRNA-share control. A pair counts once, at the locus of its higher-AS counted mate.
Prints key/value lines, then a verdict per group (gene 23S and 16S).
"""
import collections, math, sys

MIN_ID = 0.98
D23 = {'AS': 545 / 577}            # DNA div:norm 23S pairs, SRR6051694 (Deneb)
K_PREC = {'AS': {'div': 1114 / 244, 'norm': 1082 / 292}}  # T5-A2 precursor-only body:spacer
SEC16_AS = {'div': 1892 - 150 - 150, 'norm': 1691 - 150 - 150}  # GCA_009993245.1 16S lengths (as T5: s1 - s0)


def betainc(x, a, b):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    lb = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x)
    def cf(x, a, b):
        c, d = 1.0, 1 - (a + b) * x / (a + 1); d = 1 / (d if abs(d) > 1e-300 else 1e-300); h = d
        for m in range(1, 500):
            for num in (m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m)), -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))):
                d = 1 + num * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300)
                c = 1 + num / c; c = c if abs(c) > 1e-300 else 1e-300; h *= d * c
        return h
    if x < (a + 1) / (a + b + 2): return math.exp(lb) * cf(x, a, b) / a
    return 1 - math.exp(lb) * cf(1 - x, b, a) / b


def beta_ppf(q, a, b):
    lo, hi = 0.0, 1.0
    for _ in range(100):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if betainc(m, a, b) < q else (lo, m)
    return (lo + hi) / 2


def cp(k, n):
    lo = 0.0 if k == 0 else beta_ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else beta_ppf(0.975, k + 1, n - k)
    return lo, hi


def binom_sf(x, n, p):
    return 1.0 if x <= 0 else 0.0 if x > n else betainc(p, x, n - x + 1)


def verdict(d, n, D, gate=30):
    """r = (d/n)/D with two-sided 95% CP interval; RULES.md categories."""
    if d + n == 0: return 'NO READS', math.nan, (math.nan, math.nan)
    conv = lambda p: math.inf if p >= 1 else p / (1 - p) / D
    lo, hi = cp(d, d + n); r = conv(d / (d + n))
    if n < gate: v = 'UNINFORMATIVE (norm < 30)'
    elif d >= 10 and r >= 0.25: v = 'ACCUMULATES'
    elif conv(hi) < 0.10: v = 'SILENT'
    else: v = 'LOW / UNRESOLVED'
    return v, r, (conv(lo), conv(hi))


def main(paf, cls_path, total=None):
    cls = {}
    for l in open(cls_path):
        if l.startswith('target\t'): continue
        t, g, gene, c, s0, s1 = l.rstrip('\n').split('\t'); cls[t] = (g, gene, c, int(s0), int(s1))
    out = collections.Counter()
    tlen = {}

    def pair(lines):
        out['pairs_rrna85'] += 1
        mates = [x for x in lines if x[5] >= MIN_ID]
        if not mates: return
        out['pairs_rrna98'] += 1
        b = mates[0]
        for x in mates[1:]:
            if x[3] > b[3]: b = x
        t = b[0]; nm0 = b[4] == 0
        if t not in cls:
            out['decoy98'] += 1; return
        g, gene, c, s0, s1 = cls[t]
        k = f'{g}_{gene}_{c}'
        if gene == '23S':
            out[k] += 1; out[k + '_nm0'] += nm0; return
        if any(x[0] == t and x[2] >= s1 + 20 for x in mates): out[k + '_spacer'] += 1
        if all(x[0] == t and x[1] >= s0 - 10 and x[2] <= s1 + 10 for x in mates):
            out[k + '_body'] += 1; out[k + '_body_nm0'] += nm0
        else:
            out[k + '_edge'] += 1

    cur, mates = None, []
    for l in open(paf):
        f = l.rstrip('\n').split('\t')
        if 'tp:A:P' not in f[12:]: continue
        if int(f[10]) < 60 or int(f[9]) / int(f[10]) < 0.85: continue
        as_ = next(int(x[5:]) for x in f[12:] if x.startswith('AS:i:'))
        nm = next(int(x[5:]) for x in f[12:] if x.startswith('NM:i:'))
        q = f[0][:-2] if f[0][-2:] in ('/1', '/2') else f[0]
        if q != cur:
            if mates: pair(mates)
            cur, mates = q, []
        mates.append((f[5], int(f[7]), int(f[8]), as_, nm, int(f[9]) / int(f[10])))
        tlen[f[5]] = int(f[6])
    if mates: pair(mates)
    for k in sorted(out): print(k, out[k])
    if total:
        total = int(total)
        print(f'total_pairs {total}')
        print(f'C1 rRNA share (>=85% id): {out["pairs_rrna85"] / total:.4f}; (>=98% id): {out["pairs_rrna98"] / total:.4f}')
    for g in sorted({v[0] for v in cls.values()}):
        L = {}
        for t, (gg, gene, c, s0, s1) in cls.items():
            if gg == g: L.setdefault((gene, c), []).append(s1 - s0)
        for tag in ('', '_nm0'):
            d, n = out[f'{g}_23S_div{tag}'], out[f'{g}_23S_norm{tag}']
            D = D23.get(g) or (max(L[('23S', 'div')]) / max(L[('23S', 'norm')]))
            v, r, (lo, hi) = verdict(d, n, D)
            print(f'{g} 23S{tag}: div {d} norm {n} D {D:.3f} r23 {r:.3f} (95% CP {lo:.3f}-{hi:.3f}) -> {v}')
            d, n = out[f'{g}_16S_div_body{tag}'], out[f'{g}_16S_norm_body{tag}']
            sec = SEC16_AS if g == 'AS' else {c: max(L[('16S', c)]) for c in ('div', 'norm')}
            v, r, (lo, hi) = verdict(d, n, sec['div'] / sec['norm'])
            if v == 'ACCUMULATES':
                kp = K_PREC.get(g)
                if not kp: v = 'LOW / UNRESOLVED (precursor guard not calibrated)'
                else:
                    pg = binom_sf(out[f'{g}_16S_div_body'], out[f'{g}_16S_div_body'] + out[f'{g}_16S_div_spacer'], kp['div'] / (1 + kp['div']))
                    if pg > 0.01: v = f'LOW / UNRESOLVED (precursor guard p {pg:.3g})'
            print(f'{g} 16S{tag}: div body {d} norm body {n} r16 {r:.3f} (95% CP {lo:.3f}-{hi:.3f}) -> {v}; '
                  f'spacer div {out[f"{g}_16S_div_spacer"]} norm {out[f"{g}_16S_norm_spacer"]}')


if __name__ == '__main__':
    main(*sys.argv[1:])
