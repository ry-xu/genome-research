"""A3 power check: PC2 exactly as frozen, but the planted tail is a tested carrier tail (A. secundus, DASXVE01), so the
detection limit of T1 is known for the tails T1 actually tests. A7: plants each carrier tail in its own carrier genome."""
import sys, json
from multiprocessing import Pool
sys.path.insert(0, 'src')
import numpy as np
import sdscore as S
import pc2_nc1 as P

TAILS = {'GCA_009993245.1': 'UUUCAAUUGGAUA', 'GCA_036504445.1': 'UAAGGCUCUUCUU'}

def job(j):
    carrier, frac, rep = j
    tail = TAILS[carrier]
    rng = np.random.default_rng(5000 + rep)
    seqs = S.genome(carrier); ws = S.windows(carrier, seqs)
    k = round(frac * len(ws)); pick = rng.choice(len(ws), k, replace=False); ps = set(pick)
    motif = P.rrc(tail[2:10])
    ws = [(g, w[:6] + motif + w[14:]) if i in ps else (g, w) for i, (g, w) in enumerate(ws)]
    r, excl = P.score(carrier, seqs, ws, tail, f'pc2s_{carrier}_{frac}_{rep}')
    r.update(carrier=carrier, frac=frac, rep=rep, tail=tail, planted=k, planted_excl=float(excl[pick].mean()))
    return {k2: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k2, v in r.items()}

if __name__ == '__main__':
    jobs = [(c, f, r) for c in TAILS for f in (0.0, 0.03, 0.05, 0.10) for r in range(10)]
    with Pool(1) as p:
        res = p.map(job, jobs, chunksize=1)
    with open('out/pc2_self_result.jsonl', 'w') as o:
        for r in res: o.write(json.dumps(r) + '\n')
    for c in TAILS:
        for f in (0.0, 0.03, 0.05, 0.10):
            rr = [r for r in res if r['carrier'] == c and r['frac'] == f]
            print(c, TAILS[c], f, 'calls', sum(r['call'] for r in rr), '/', len(rr), 'median z', round(float(np.median([r['z'] for r in rr])), 2), 'median planted excl', round(float(np.median([r['planted_excl'] for r in rr])), 2))
