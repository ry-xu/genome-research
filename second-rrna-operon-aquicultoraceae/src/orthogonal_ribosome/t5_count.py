"""T5 counts (RULES_T5.md): divergent vs normal 16S read pairs, plus the 23S replicate (PC-T5b).

usage: uv run python src/orthogonal_ribosome/t5_count.py PAF CLASS16.tsv [--sim]
PAF holds minimap2 -c lines already filtered to identity >= 0.98 and >= 60 aligned columns (fastq-dump --split-spot:
both mates share the read name). A pair counts once, at the locus of its higher-AS mate.
"""
import collections, math, sys

CLS23 = 'data/expression_23s/loci_class_all.tsv'  # repo copy (lab: work/frozen-core/results/lsu/bac/cca/expr/)
K_PREC = {'div': 1114 / 244, 'norm': 1082 / 292}  # T5-A2: body : spacer pairs from precursor-only reads (out/t5/prec)
CG = ['PFKS', 'PFIW', 'MNXI', 'JAACTV', 'PFTP', 'PFNG', 'PEXG', 'PFFU']  # the 8 Crystal Geyser genomes' contig prefixes


def betainc(x, a, b):
    """Regularized incomplete beta I_x(a, b) (continued fraction), no scipy."""
    if True:
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
    """Beta quantile by bisection on betainc."""
    lo, hi = 0.0, 1.0
    for _ in range(100):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if betainc(m, a, b) < q else (lo, m)
    return (lo + hi) / 2


def binom_sf(x, n, p):
    """P(X >= x), X ~ Binomial(n, p)."""
    return 1.0 if x <= 0 else 0.0 if x > n else betainc(p, x, n - x + 1)


def cp(k, n):
    lo = 0.0 if k == 0 else beta_ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else beta_ppf(0.975, k + 1, n - k)
    return lo, hi


def main(paf, cls16_path):
    c16 = {}
    for l in open(cls16_path):
        n, c, s0, s1 = l.rstrip('\n').split('\t'); c16[n] = (c, int(s0), int(s1))
    c23 = {}
    for l in open(CLS23):
        f = l.split('\t')
        if any(x in f[0].split('|')[1] for x in CG): c23[f[0]] = f[1]
    out = collections.Counter()

    def pair(mates):
        b = mates[0]
        for x in mates[1:]:
            if x[3] > b[3]: b = x  # higher-AS mate (first wins ties)
        t = b[0]; nm0 = b[4] == 0
        if t in c16:
            c, s0, s1 = c16[t]
            if c not in ('div', 'norm'):
                out[f'16S_{c}'] += 1; return
            if any(x[0] == t and x[2] >= s1 + 20 for x in mates): out[f'16S_{c}_spacer'] += 1
            if all(x[0] == t and x[1] >= s0 - 10 and x[2] <= s1 + 10 for x in mates):
                out[f'16S_{c}_body'] += 1; out[f'16S_{c}_body_nm0'] += nm0
            else:
                out[f'16S_{c}_edge'] += 1
        elif t in c23:
            out[f'23S_{c23[t]}'] += 1; out[f'23S_{c23[t]}_nm0'] += nm0
        else:
            out['other_target'] += 1

    # streamed: minimap2 writes a pair's records consecutively, so a pair is a run of lines with one name
    cur, mates, seen = None, [], 0
    for l in open(paf):
        f = l.rstrip('\n').split('\t')
        if 'tp:A:P' not in f[12:]: continue
        as_ = next(int(x[5:]) for x in f[12:] if x.startswith('AS:i:'))
        nm = next(int(x[5:]) for x in f[12:] if x.startswith('NM:i:'))
        q = f[0][:-2] if f[0][-2:] in ("/1", "/2") else f[0]  # minimap2 tags interleaved same-name mates /1 /2
        if q != cur:
            if mates: pair(mates); seen += 1
            cur, mates = q, []
        mates.append((f[5], int(f[7]), int(f[8]), as_, nm))
    if mates: pair(mates); seen += 1
    out['pairs_seen'] = seen
    for k in sorted(out): print(k, out[k])
    sec = {c: s1 - s0 for n, (c, s0, s1) in c16.items() if n.startswith('GCA_009993245.1|') and c in ('div', 'norm')}
    verdicts = []
    for tag in ('', '_nm0'):
        d, n = out[f'16S_div_body{tag}'], out[f'16S_norm_body{tag}']
        if d + n == 0:
            print(f'verdict{tag}: no reads'); verdicts.append('NO READS'); continue
        conv = lambda p: math.inf if p >= 1 else p / (1 - p) * sec['norm'] / sec['div']
        lo, hi = cp(d, d + n); r = conv(d / (d + n))
        # precursor guard on all-read counts (spacer pairs are not split by NM): body pairs beyond unprocessed transcript
        ds = out['16S_div_spacer']; kd = K_PREC['div']; pg = binom_sf(out['16S_div_body'], out['16S_div_body'] + ds, kd / (1 + kd)) if kd else None
        if n < 30: v = 'UNINFORMATIVE (G1: norm < 30)'
        elif d >= 10 and r >= 0.25: v = 'ACCUMULATES' if pg is not None and pg <= 0.01 else 'LOW / UNRESOLVED (precursor guard)'
        elif conv(hi) < 0.10: v = 'DOES NOT ACCUMULATE'
        else: v = 'LOW / UNRESOLVED'
        verdicts.append(v)
        corr = {c: out[f'16S_{c}_body'] - (K_PREC[c] or 0) * out[f'16S_{c}_spacer'] for c in ('div', 'norm')}
        rc = corr['div'] / sec['div'] / (corr['norm'] / sec['norm']) if corr['norm'] > 0 else math.nan
        ralt = r * sec['div'] / (sec['div'] - 195)  # sensitivity: div insertions (91+70+34 nt) excised from mature RNA
        print(f'verdict{tag}: div {d} norm {n} r {r:.3f} (95% two-sided CP {conv(lo):.3f}-{conv(hi):.3f}) -> {v}')
        print(f'  precursor guard p {pg}; precursor-corrected body div {corr["div"]:.1f} norm {corr["norm"]:.1f} r {rc:.3f}; '
              f'r if insertions excised {ralt:.3f}; body:spacer div {out["16S_div_body"]}:{ds} norm {out["16S_norm_body"]}:{out["16S_norm_spacer"]}')
    print('VERDICT (T5):', verdicts[0] if len(set(verdicts)) == 1 else f'UNRESOLVED (primary {verdicts[0]}, NM=0 {verdicts[1]})')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
