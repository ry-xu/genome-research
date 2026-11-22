"""A5 power check for the 4 windows-fallback carriers: in non-carrier GCA_018333515.1, one 13-nt window of the carrier's
70-nt block (chosen at random, seed 7000+rep) gets its positions 3-10 planted (reverse complement at -14..-7) in a share
of eligible genes; the T1 block statistic (max E over the 58 windows) is then called against 1,000 decoy blocks.
Deviation (A5): decoy blocks are scored once on the unplanted genome per block, not per replicate (cost)."""
import sys, json, os
sys.path.insert(0, 'src')
import numpy as np
import sdscore as S
import pc2_nc1 as P
import t1_run as T

BLOCKS = ['GCA_902810715.1', 'GCA_014859275.1', 'GCA_035643475.1', 'GCA_027392615.1']


def job(j):
    carrier, frac, rep = j
    kind, blk = T.tested(carrier); assert kind == 'block'
    acc = P.PC2_ACC; seqs = S.genome(acc); ws = S.windows(acc, seqs); wl = [w for _, w in ws]
    can0, _ = S.duplex(S.CANON, wl); weak0 = np.where(can0 > S.CAN_WEAK)[0]
    dec = T.decoy_blocks(acc, seqs, [wl[i] for i in weak0], len(wl), blk)  # unplanted genome, cached per block
    rng = np.random.default_rng(7000 + rep)
    i0 = int(rng.integers(0, 58)); motif = P.rrc(blk[i0 + 2:i0 + 10])
    k = round(frac * len(wl)); pick = rng.choice(len(wl), k, replace=False); ps = set(pick)
    wlp = [w[:6] + motif + w[14:] if i in ps else w for i, w in enumerate(wl)]
    can, _ = S.duplex(S.CANON, wlp); weak = np.where(can > S.CAN_WEAK)[0]
    e, n, npos = T.block_stats(blk, [wlp[i] for i in weak], len(wlp))
    r = S.call(e, n, npos, dec)
    r.update(carrier=carrier, frac=frac, rep=rep, planted_window=i0, planted=k)
    return {a: (float(b) if isinstance(b, (np.floating, np.integer)) else b) for a, b in r.items()}


if __name__ == '__main__':
    jobs = [(c, f, r) for c in BLOCKS for f in (0.0, 0.03, 0.05, 0.10) for r in range(5)]
    out = 'out/pc2_block_result.jsonl'  # appended per job, so a stopped run resumes where it left off
    res = [json.loads(l) for l in open(out)] if os.path.exists(out) else []
    done = {(r['carrier'], r['frac'], r['rep']) for r in res}
    with open(out, 'a') as o:
        for j in jobs:
            if j in done: continue
            r = job(j); res.append(r); o.write(json.dumps(r) + '\n'); o.flush()
            print(j, 'z', round(r['z'], 2), 'p', r['p'], 'call', r['call'], flush=True)
    for c in BLOCKS:
        for f in (0.0, 0.03, 0.05, 0.10):
            rr = [r for r in res if r['carrier'] == c and r['frac'] == f]
            print(c, f, 'calls', sum(r['call'] for r in rr), '/', len(rr), 'median z', round(float(np.median([r['z'] for r in rr])), 2))
