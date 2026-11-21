"""Conserved-site codon-meaning flags for census_org.py censuses, against the STANDARD code.

usage: python sitescan_org.py out_prefix census.npz [census2.npz ...]   (all made with the same profiles)
env: STRICT (default 0.9), MIN_SITES (default 6), NULL=seed (noise null: each genome's site codons are shuffled
     among sites whose codons the declared table reads as the same residue; the column consensus is unchanged)
Column consensus comes from the data, not the HMM (HMMER's entropy weighting leaves only ~2% of columns with a max
emission >= 0.7): each (profile, column)'s residue frequencies over all genomes, each genome's residue read with its
declared table (sense-to-stop erased; stops not counted), each genome weighted 1/(genomes in its family, taxonomy
minus the last rank). A column is strict if one residue has >= STRICT of the weight; a site there votes for it.
Node level (each taxonomy prefix with >= 3 genomes): a codon's sites are the distinct strict sites where >= 2 of the
node's genomes have it (>= 1 if the node has < 5 genomes). Genome level: a genome's own distinct strict sites.
Flag: >= MIN_SITES sites, the standard residue (for a stop: none) has <= 20% of them, one other residue >= 60%.
Node flags are reported at origin nodes only (dropped if the parent has the same flag).
Status: declared  - the top residue is what the declared table says (for a node: >= 80% of the genomes that have
                    the codon at strict sites; the rest show up as genome-level declaration lapses);
        editing   - Embryophyta, or one C<->U change at any codon position gives a standard codon for the top residue;
        candidate - neither.
Writes out_prefix.node.tsv, out_prefix.genome.tsv.
"""
import sys, os
from collections import defaultdict, Counter
import numpy as np
from Bio.Data import CodonTable
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from codes import blinded

AA = 'ACDEFGHIKLMNPQRSTVWY'
B = 'TCAG'
CODONS = [a + b + c for a in B for b in B for c in B]
STD = CodonTable.unambiguous_dna_by_id[1].forward_table
STRICT = float(os.environ.get('STRICT', 0.9))
MIN_SITES = int(os.environ.get('MIN_SITES', 6))


def meaning(t):
    try:
        f = CodonTable.unambiguous_dna_by_id[int(t)].forward_table
    except KeyError:
        f = STD
    return {c: f.get(c, '*') for c in CODONS}


def edit_ok(c, aa):
    for i in range(3):
        for a, b in (('C', 'T'), ('T', 'C')):
            if c[i] == a and STD.get(c[:i] + b + c[i + 1:]) == aa:
                return True
    return False


def load(paths, lean=False):
    """-> dict: acc, tax, table, models, sg, sc, col (model*K+column), res (declared residue index, -1 stop),
    cons [M*K], strict [M*K], w [G]
    lean=True: the same values in narrow types (sg, col int32; sc, res int8), the column sums done CHUNK sites at a
    time in the same order (so the float sums are the same); the default widens every site to int64 (~90 bytes a
    site; 85 M sites peaked at 7.7 GB)."""
    if lean:
        return load_lean(paths)
    acc, tax, table, sg, sc, sm, sk, models = [], [], [], [], [], [], [], None
    for p in paths:
        d = np.load(p)
        off = len(acc)
        acc += d['acc'].tolist(); tax += d['taxonomy'].tolist(); table += [int(t) for t in d['table']]
        sg.append(d['site_genome'] + off); sc.append(d['site_codon']); sm.append(d['site_model']); sk.append(d['site_col'])
        models = d['models'].tolist()
    sg, sc, sm, sk = map(np.concatenate, (sg, sc, sm, sk))
    sc = sc.astype(np.int64) % 64
    fam = Counter(t.rsplit(';', 1)[0] for t in tax)
    w = np.array([1 / fam[t.rsplit(';', 1)[0]] for t in tax])
    ri = {a: i for i, a in enumerate(AA)}
    tabres = {t: np.array([ri.get(blinded(t)[c], -1) for c in CODONS]) for t in set(table)}
    gt = np.array(table)
    res = np.full(len(sg), -1)
    for t, arr in tabres.items():
        m = gt[sg] == t
        res[m] = arr[sc[m]]
    ok = res >= 0
    M = len(models); K = int(sk.max()) + 1
    col = sm.astype(np.int64) * K + sk
    freq = np.zeros((M * K, 20))
    np.add.at(freq, (col[ok], res[ok]), w[sg[ok]])
    tot = freq.sum(1)
    strict = (tot > 0) & (freq.max(1) >= STRICT * np.maximum(tot, 1e-12))
    return dict(acc=acc, tax=tax, table=table, models=models, K=K, sg=sg, sc=sc, col=col, res=res,
                cons=freq.argmax(1), strict=strict, freq=freq, tot=tot)


CHUNK = 4_000_000


def load_lean(paths):
    acc, tax, table, sg, sc, sm, sk, models = [], [], [], [], [], [], [], None
    for p in paths:
        d = np.load(p)
        off = len(acc)
        acc += d['acc'].tolist(); tax += d['taxonomy'].tolist(); table += [int(t) for t in d['table']]
        sg.append((d['site_genome'] + off).astype(np.int32)); sc.append(d['site_codon'])
        sm.append(d['site_model']); sk.append(d['site_col'])
        models = d['models'].tolist()
    sg = np.concatenate(sg)
    sc = (np.concatenate(sc).astype(np.int16) % 64).astype(np.int8)
    sm, sk = np.concatenate(sm), np.concatenate(sk)
    fam = Counter(t.rsplit(';', 1)[0] for t in tax)
    w = np.array([1 / fam[t.rsplit(';', 1)[0]] for t in tax])
    ri = {a: i for i, a in enumerate(AA)}
    ts = sorted(set(table))
    tabres = np.array([[ri.get(blinded(t)[c], -1) for c in CODONS] for t in ts], np.int8)
    gt = np.array([ts.index(t) for t in table], np.int16)
    M = len(models); K = int(sk.max()) + 1
    col = sm.astype(np.int32) * K + sk
    del sm, sk
    res = np.empty(len(sg), np.int8)
    freq = np.zeros((M * K, 20))
    for i in range(0, len(sg), CHUNK):
        s = slice(i, i + CHUNK)
        r = tabres[gt[sg[s]], sc[s]]
        res[s] = r
        ok = r >= 0
        np.add.at(freq, (col[s][ok], r[ok]), w[sg[s][ok]])
    tot = freq.sum(1)
    strict = (tot > 0) & (freq.max(1) >= STRICT * np.maximum(tot, 1e-12))
    return dict(acc=acc, tax=tax, table=table, models=models, K=K, sg=sg, sc=sc, col=col, res=res,
                cons=freq.argmax(1), strict=strict, freq=freq, tot=tot)


def shuffle_synonymous(D, seed):
    """noise null: permute codons among each genome's sites with the same declared residue (stops among stops)"""
    rng = np.random.default_rng(seed)
    sc = D['sc'].copy()
    key = D['sg'] * 21 + (D['res'] + 1)
    order = np.argsort(key, kind='stable')
    k = key[order]
    bounds = np.flatnonzero(np.diff(k)) + 1
    for grp in np.split(order, bounds):
        if len(grp) > 1:
            sc[grp] = sc[rng.permutation(grp)]
    D['sc'] = sc


def main(out, *paths):
    D = load(paths)
    if os.environ.get('NULL'):
        shuffle_synonymous(D, int(os.environ['NULL']))
    acc, tax, table, sg, sc, col, cons, strict = (D[k] for k in ('acc', 'tax', 'table', 'sg', 'sc', 'col', 'cons', 'strict'))
    keep = strict[col]
    print(f'{len(acc)} genomes, {len(sg)} sites, {strict.sum()} strict columns of {(D["tot"] > 0).sum()}, '
          f'{keep.mean():.3f} of sites strict', file=sys.stderr)
    per = defaultdict(lambda: defaultdict(set))
    for g, c, k in zip(sg[keep], sc[keep], col[keep]):
        per[g][c].add(k)
    mean = {t: meaning(t) for t in set(table)}

    def flag(sites):
        if len(sites) < MIN_SITES:
            return None
        cnt = Counter(AA[cons[k]] for k in sites)
        top, n = cnt.most_common(1)[0]
        return top, n / len(sites), cnt, len(sites)

    def status(c, top, gs, lineage):
        cod = CODONS[c]
        if sum(mean[table[i]][cod] == top for i in gs) >= 0.8 * len(gs):
            return 'declared'
        if 'Embryophyta' in lineage or edit_ok(cod, top):
            return 'editing'
        return 'candidate'

    hdr = 'codon\tstd\tdeclared\tsite_aa\tsites\ttop_share\tstd_share\tvotes\tstatus'
    with open(out + '.genome.tsv', 'w') as fh:
        fh.write('acc\ttable\ttaxonomy\t' + hdr + '\n')
        for g in sorted(per):
            for c, sites in sorted(per[g].items()):
                f = flag(sites)
                if f is None:
                    continue
                top, sh, cnt, U = f
                std = STD.get(CODONS[c], '*')
                ssh = cnt.get(std, 0) / U
                if top == std or sh < 0.6 or ssh > 0.2:
                    continue
                fh.write(f"{acc[g]}\t{table[g]}\t{tax[g]}\t{CODONS[c]}\t{std}\t{mean[table[g]][CODONS[c]]}\t{top}\t{U}\t"
                         f"{sh:.2f}\t{ssh:.2f}\t{dict(cnt.most_common(4))}\t{status(c, top, [g], tax[g])}\n")
    nodes = defaultdict(list)
    for i, t in enumerate(tax):
        p = t.split(';')
        for k in range(2, len(p) + 1):
            nodes[';'.join(p[:k])].append(i)
    flags = {}
    for node, ii in nodes.items():
        if len(ii) < 3:
            continue
        need = 2 if len(ii) >= 5 else 1
        for c in range(64):
            sup, gws = Counter(), []
            for i in ii:
                s = per[i].get(c)
                if s:
                    sup.update(s)
                    gws.append(i)
            f = flag({x for x, n in sup.items() if n >= need})
            if f is None:
                continue
            top, sh, cnt, U = f
            std = STD.get(CODONS[c], '*')
            ssh = cnt.get(std, 0) / U
            if top == std or sh < 0.6 or ssh > 0.2:
                continue
            decl = Counter(mean[table[i]][CODONS[c]] for i in gws)
            flags[(node, c)] = (len(ii), len(gws), top, sh, ssh, cnt, U, decl, status(c, top, gws, node))
    with open(out + '.node.tsv', 'w') as fh:
        fh.write('node\tgenomes\tgenomes_with_codon\t' + hdr + '\n')
        for (node, c), f in sorted(flags.items()):
            parent = node.rsplit(';', 1)[0]
            if (parent, c) in flags and flags[(parent, c)][2] == f[2]:
                continue
            g, gw, top, sh, ssh, cnt, U, decl, st = f
            fh.write(f"{node}\t{g}\t{gw}\t{CODONS[c]}\t{STD.get(CODONS[c], '*')}\t{dict(decl)}\t{top}\t{U}\t{sh:.2f}\t"
                     f"{ssh:.2f}\t{dict(cnt.most_common(4))}\t{st}\n")


if __name__ == '__main__':
    main(*sys.argv[1:])
