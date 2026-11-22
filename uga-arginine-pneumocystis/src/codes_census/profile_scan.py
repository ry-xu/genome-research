"""All-columns codon-meaning scan for census_org.py censuses (complements sitescan_org.py, which uses strict columns).

usage: python profile_scan.py out.tsv census.npz [census2.npz ...]
env: NULL=seed (shuffle each genome's codons among sites with the same declared residue first, as in sitescan_org),
     MIN_N (default 10), MIN_LLR (default 10)
For each genome, n[c, y] = sites with codon c whose column consensus (sitescan_org.load) is residue y, over ALL
aligned columns. For codon c declared X, and every residue Y, the profile q_Y is the consensus distribution under the
genome's other codons declared Y (for Y = X, c itself is left out), smoothed toward the genome's overall consensus
distribution (pseudo-weight ALPHA = 20 sites). LLR(Y) = sum_y n[c, y] log(q_Y[y] / q_X[y]): how much better c's
positions fit Y's codons than its own synonyms'. Writes every (genome, codon) with n >= MIN_N and best LLR >= MIN_LLR
for some Y != X."""
import sys, os, zipfile
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sitescan_org import load, AA, CODONS, shuffle_synonymous, blinded

ALPHA = 20.0
CHUNK = 4_000_000


def npy_chunks(z, name):
    """yield a 1-D .npy member of an open .npz (zip) in CHUNK-element pieces, without loading it whole"""
    with z.open(name) as h:
        v = np.lib.format.read_magic(h)
        shape, _, dt = (np.lib.format.read_array_header_1_0 if v == (1, 0) else np.lib.format.read_array_header_2_0)(h)
        for i in range(0, shape[0], CHUNK):
            n = min(CHUNK, shape[0] - i)
            yield np.frombuffer(h.read(n * dt.itemsize), dt, n)


def lean_counts(paths):
    """Same acc, tax, table, N (as int32) and decl as load() + the bincount below, but one census file and CHUNK sites
    at a time with narrow integer types (load() widens ~85 M sites to int64 and peaked at 7.7 GB). The column
    consensus uses np.add.at over the sites in the same order as load(), so the float sums are the same."""
    acc, tax, table, M, K = [], [], [], 0, 0
    for p in paths:
        d = np.load(p)
        acc += d['acc'].tolist(); tax += d['taxonomy'].tolist(); table += [int(t) for t in d['table']]
        M = len(d['models'])
        with zipfile.ZipFile(p) as z:
            K = max(K, max(int(k.max()) for k in npy_chunks(z, 'site_col.npy')) + 1)
    fam = Counter(t.rsplit(';', 1)[0] for t in tax)
    w = np.array([1 / fam[t.rsplit(';', 1)[0]] for t in tax])
    ri = {a: i for i, a in enumerate(AA)}
    gt = np.array(table)
    tabres = np.full((int(gt.max()) + 1, 64), -1, np.int8)
    for t in set(table):
        tabres[t] = [ri.get(blinded(t)[c], -1) for c in CODONS]

    def sites():  # (sg int64, sc int64 in 0..63, col int64, res int8) per chunk, in load()'s order
        off = 0
        for p in paths:
            with zipfile.ZipFile(p) as z:
                st = [npy_chunks(z, f'site_{x}.npy') for x in ('genome', 'codon', 'model', 'col')]
                for g, c, m, k in zip(*st):
                    sg = g.astype(np.int64) + off
                    sc = c.astype(np.int64) % 64
                    yield sg, sc, m.astype(np.int64) * K + k, tabres[gt[sg], sc]
                off += len(np.load(p)['acc'])

    freq = np.zeros((M * K, 20))
    for sg, sc, col, res in sites():
        ok = res >= 0
        np.add.at(freq, (col[ok], res[ok]), w[sg[ok]])
    cons = freq.argmax(1)
    G = len(acc)
    N = np.zeros(G * 64 * 20, np.int32)
    decl = np.full((G, 64), -1, int)
    for sg, sc, col, res in sites():
        N += np.bincount((sg * 64 + sc) * 20 + cons[col], minlength=G * 64 * 20).astype(np.int32)
        decl[sg, sc] = res  # constant per (genome, codon); -1 = stop
    return acc, tax, table, N.reshape(G, 64, 20), decl


def main(out, *paths):
    min_n, min_llr = int(os.environ.get('MIN_N', 10)), float(os.environ.get('MIN_LLR', 10))
    if os.environ.get('NULL'):  # the shuffle needs every site at once: original (memory-heavy) path
        D = load(paths)
        shuffle_synonymous(D, int(os.environ['NULL']))
        acc, tax, table, sg, sc, col, cons, res = (D[k] for k in ('acc', 'tax', 'table', 'sg', 'sc', 'col', 'cons', 'res'))
        G = len(acc)
        y = cons[col]
        N = np.bincount((sg.astype(np.int64) * 64 + sc) * 20 + y, minlength=G * 64 * 20).reshape(G, 64, 20)
        decl = np.full((G, 64), -1, int)
        decl[sg, sc] = res  # constant per (genome, codon); -1 = stop
    else:
        acc, tax, table, N, decl = lean_counts(paths)
        G = len(acc)
    A = np.zeros((G, 20, 20))
    for c in range(64):
        d = decl[:, c]
        ok = d >= 0
        np.add.at(A, (np.flatnonzero(ok), d[ok]), N[ok, c])
    bg = N.sum(1) + 1.0
    bg /= bg.sum(1, keepdims=True)
    rows = 0
    with open(out, 'w') as fh:
        fh.write('acc\ttable\ttaxonomy\tcodon\tdeclared\tbest\tn\tllr\tllr_per_site\tq_best_n\tq_decl_n\ttop_cons\n')
        for g in range(G):
            for c in range(64):
                X = decl[g, c]
                n = N[g, c].astype(float)
                nt = n.sum()
                if X < 0 or nt < min_n:
                    continue
                prof = A[g].copy()
                prof[X] -= n
                tot = prof.sum(1)
                if tot[X] < 1:
                    continue  # no synonyms in this genome (e.g. AUG): nothing to compare against
                q = (prof + ALPHA * bg[g]) / (tot[:, None] + ALPHA)
                ll = (n * np.log(q)).sum(1)
                llr = ll - ll[X]
                llr[X] = -np.inf
                llr[tot < 1] = -np.inf
                Y = int(np.argmax(llr))
                if llr[Y] < min_llr:
                    continue
                top = ','.join(f'{AA[i]}{int(n[i])}' for i in np.argsort(-n)[:4] if n[i] > 0)
                fh.write(f'{acc[g]}\t{table[g]}\t{tax[g]}\t{CODONS[c]}\t{AA[X]}\t{AA[Y]}\t{int(nt)}\t{llr[Y]:.1f}\t'
                         f'{llr[Y] / nt:.2f}\t{int(tot[Y])}\t{int(tot[X])}\t{top}\n')
                rows += 1
    print(f'{G} genomes, {rows} flagged (n >= {min_n}, LLR >= {min_llr})', file=sys.stderr)


if __name__ == '__main__':
    main(*sys.argv[1:])
