"""PC2 (planted orthogonal pair in Aquicultor GCA_018333515.1, 20 replicates at 3% and 1%) and NC1 (20 random
unplanted tails on each of the 6 divergent-16S carriers). RULES.md positive controls; canonical tail GAUCACCUCCUUU."""
import sys, json
from multiprocessing import Pool
sys.path.insert(0, 'src')
import numpy as np
import sdscore as S

PC2_ACC = 'GCA_018333515.1'
CARRIERS = ['GCA_036504445.1', 'GCA_009993245.1', 'GCA_902810715.1', 'GCA_014859275.1', 'GCA_035643475.1', 'GCA_027392615.1']


def rand_tail(rng):
    while True:
        t = ''.join(rng.choice(list('ACGU'), 13))
        if 5 <= S.gc(t) <= 8 and not S.has_core(t):
            return t


def rrc(s):  # RNA reverse complement
    return s.translate(str.maketrans('ACGU', 'UGCA'))[::-1]


def score(acc, seqs, ws, tail, tag):
    wl = [w for _, w in ws]
    can, _ = S.duplex(S.CANON, wl)
    dg, mid = S.duplex(tail, wl)
    e, n, npos = S.stat(dg, mid, can, exclusive=True)
    dec = S.decoy_stats(acc, seqs, ws, can, S.gc(tail), S.has_core(tail), exclusive=True,
                        cache_dir=f'out/cache/{tag}')
    r = S.call(e, n, npos, dec)
    excl = (dg <= S.STRONG) & (can > S.CAN_WEAK)
    return r, excl


def pc2(job):
    frac, rep = job
    rng = np.random.default_rng(1000 + rep)
    tail = rand_tail(rng)
    seqs = S.genome(PC2_ACC)
    ws = S.windows(PC2_ACC, seqs)
    k = max(1, round(frac * len(ws)))
    pick = rng.choice(len(ws), k, replace=False)
    motif = rrc(tail[2:10])
    ws = [(g, w[:6] + motif + w[14:]) if i in set(pick) else (g, w) for i, (g, w) in enumerate(ws)]
    r, excl = score(PC2_ACC, seqs, ws, tail, f'pc2_{frac}_{rep}')
    r.update(kind='PC2', frac=frac, rep=rep, acc=PC2_ACC, tail=tail, planted=k, planted_excl=float(excl[pick].mean()))
    return r


def nc1(job):
    acc, rep = job
    rng = np.random.default_rng(2000 + rep)
    tail = rand_tail(rng)
    seqs = S.genome(acc)
    ws = S.windows(acc, seqs)
    r, _ = score(acc, seqs, ws, tail, 'nc1')  # unplanted windows: decoys depend only on genome and GC, so shared
    r.update(kind='NC1', rep=rep, acc=acc, tail=tail)
    return r


if __name__ == '__main__':
    what = sys.argv[1]
    fracs = [float(x) for x in sys.argv[2].split(',')] if what == 'pc2' and len(sys.argv) > 2 else [0.03, 0.01]
    jobs = [(f, r) for f in fracs for r in range(20)] if what == 'pc2' else [(a, r) for a in CARRIERS for r in range(20)]
    with Pool(1) as p:
        res = p.map(pc2 if what == 'pc2' else nc1, jobs, chunksize=1)
    with open(f'out/{what}_result.jsonl', 'a' if what == 'pc2' and len(sys.argv) > 2 else 'w') as o:
        for r in res:
            o.write(json.dumps({k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in r.items()}) + '\n')
    if what == 'pc2':
        for f in fracs:
            rr = [r for r in res if r['frac'] == f]
            print(f, 'calls', sum(r['call'] for r in rr), '/20; median planted exclusive', np.median([r['planted_excl'] for r in rr]))
    else:
        print('NC1 calls', sum(r['call'] for r in res), '/', len(res))
