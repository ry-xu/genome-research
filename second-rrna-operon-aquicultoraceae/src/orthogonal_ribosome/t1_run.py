"""T1 channel test on the 6 divergent-16S carriers (RULES.md T1 + A1). Defined tails: E(tail) vs 1,000 GC-matched
decoy tails. Undefined tails (A1: < 7 of 13 columns held): max E over the 58 13-nt windows of the 70-nt block
(last 40 nt of the divergent 16S model hit + 30 nt downstream) vs 1,000 genomic 70-nt decoy blocks (GC fraction
+-0.05, without CCUC/CUCC unless the block has one, both strands), each scored as its own max over 58 windows.
Canonical tail GAUCACCUCCUUU (identical in all Aquicultoraceae normal 16S read). Usage: t1_run.py ACC [ACC...]
(non-carrier controls: t1_run.py --tail-from CARRIER ACC...)."""
import sys, json, os
from multiprocessing import Pool
sys.path.insert(0, 'src')
import numpy as np
import sdscore as S

DIV = {}  # carrier -> (contig, end3, strand) of the divergent 16S model hit (model end 1533)
for l in open('out/t0/ssu.tbl'):
    if l[0] == '#': continue
    f = l.split()
    if int(f[6]) == 1533 and float(f[14]) < 1500 and int(f[5]) == 1 and float(f[14]) > 900:
        acc, contig = f[0].split('|')
        DIV[acc] = (contig, int(f[8]), f[9])
TAILS = {l.split('\t')[0]: l.strip().split('\t') for l in open('out/t0/tails.tsv') if not l.startswith('acc')}


def tested(carrier):
    """('tail', 13-mer) or ('block', 70-mer) for a carrier."""
    _, _, held, vt, _ = TAILS[carrier]
    if int(held) >= 7:
        return 'tail', vt
    contig, b, st = DIV[carrier]
    s = S.genome(carrier)[contig]
    blk = s[b - 40:b + 30] if st == '+' else S.rc(s[b - 31:b + 39])
    assert len(blk) == 70
    return 'block', S.rna(blk)


def excl_stats(tail, wl_weak, n_all):
    d, m = S.duplex(tail, wl_weak)
    hit = d <= S.STRONG
    return hit.sum() / n_all, int(hit.sum()), int((hit & (m >= S.POS_LO) & (m <= S.POS_HI)).sum())


def block_stats(blk, wl_weak, n_all):
    return max((excl_stats(blk[i:i + 13], wl_weak, n_all) for i in range(len(blk) - 12)), key=lambda x: x[0])


def decoy_blocks(acc, seqs, wl_weak, n_all, blk, n=1000, seed=0):
    f = f'out/cache/t1blk/{acc}_{blk}.npy'
    os.makedirs('out/cache/t1blk', exist_ok=True)
    if os.path.exists(f):
        return np.load(f)
    rng = np.random.default_rng(seed)
    g = S.gc(blk) / 70
    bl = S.random_kmers(seqs, n, 70, g - 0.05, g + 0.05, S.has_core(blk), rng, frac=True)
    res = np.array([block_stats(b, wl_weak, n_all) for b in bl], float)
    np.save(f, res)
    return res


def run(job):
    carrier, acc = job
    kind, t = tested(carrier)
    seqs = S.genome(acc)
    ws = S.windows(acc, seqs)
    wl = [w for _, w in ws]
    can, _ = S.duplex(S.CANON, wl)
    weak = np.where(can > S.CAN_WEAK)[0]
    wl_weak = [wl[i] for i in weak]
    if kind == 'tail':
        e, n, npos = excl_stats(t, wl_weak, len(wl))
        dec = S.decoy_stats(acc, seqs, ws, can, S.gc(t), S.has_core(t), exclusive=True, cache_dir='out/cache/t1')
        best = t
        d, _ = S.duplex(t, wl_weak)
        genes = [ws[weak[i]][0] for i in np.where(d <= S.STRONG)[0]]
    else:
        i = max(range(58), key=lambda i: excl_stats(t[i:i + 13], wl_weak, len(wl))[0])
        best = t[i:i + 13]
        e, n, npos = excl_stats(best, wl_weak, len(wl))
        dec = decoy_blocks(acc, seqs, wl_weak, len(wl), t)
        d, _ = S.duplex(best, wl_weak)
        genes = [ws[weak[i]][0] for i in np.where(d <= S.STRONG)[0]]
    r = S.call(e, n, npos, dec)
    r.update(carrier=carrier, acc=acc, kind=kind, tested=t, best_window=best, genes=len(wl), weak_can=len(weak),
             excl_genes=genes)
    return r


if __name__ == '__main__':
    a = sys.argv[1:]
    if a[0] == '--tail-from':
        jobs = [(a[1], x) for x in a[2:]]; out = f'out/t1/ctrl_{a[1]}.jsonl'
    else:
        jobs = [(x, x) for x in a]; out = 'out/t1/carriers.jsonl'
    os.makedirs('out/t1', exist_ok=True)
    with Pool(min(2, len(jobs))) as p:
        res = p.map(run, jobs, chunksize=1)
    with open(out, 'a') as o:
        for r in res:
            o.write(json.dumps({k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in r.items()}) + '\n')
            print(r['carrier'], r['acc'], r['kind'], r['best_window'], f"E={r['E']:.4f} dec={r['dec_mean']:.4f} p={r['p']:.4f} z={r['z']:.2f} pos={r['pos_ratio']:.2f} call={r['call']}")
