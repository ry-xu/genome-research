"""Shared scoring for RULES.md T1/PC1/PC2/NC1.

Windows: the 20 nt upstream of each eligible gene's start (complete start, >= 300 nt CDS, >= 20 nt upstream on the
contig), as RNA. Duplex: ViennaRNA duplexfold(tail, window), 37 C defaults.
Per gene: strong if dG(tail) <= -7.0; exclusive if strong and dG(canonical) > -4.5.
Decoy tails: random 13-mers of the genome (both strands, ACGT only), GC count within +-1 of the tested tail,
without CCUC/CUCC unless the tested tail has one. Decoy blocks (windows fallback): random 70-nt genomic blocks,
GC fraction within +-0.05, scored as the max over their 58 13-nt windows.
"""
import gzip
import os
import random

import numpy as np
import pyarrow.parquet as pq
import RNA

ATLAS = 'data/external/codon-atlas/'
STRONG, CAN_WEAK = -7.0, -4.5
POS_LO, POS_HI = -13, -8  # RULES A2 (was -17, -4)
COMP = str.maketrans('ACGTN', 'TGCAN')
CANON = 'GAUCACCUCCUUU'


def rc(s):
    return s.translate(COMP)[::-1]


def rna(s):
    return s.upper().replace('T', 'U')


def has_core(t):
    return 'CCUC' in t or 'CUCC' in t


def gc(s):
    return sum(c in 'GC' for c in s)


def genome(acc):
    idx = {l.split('\t')[0]: l.split('\t')[1] for l in open(ATLAS + 'index/genomes.tsv')}
    seqs, n = {}, None
    for l in gzip.open(ATLAS + idx[acc], 'rt'):
        if l[0] == '>':
            n = l[1:].split()[0]; seqs[n] = []
        else:
            seqs[n].append(l.strip().upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def windows(acc, seqs=None):
    """[(gene_id, window_rna)] for eligible genes."""
    seqs = seqs or genome(acc)
    g = pq.read_table(f'{ATLAS}products/{acc[:8]}/{acc}/genes.parquet').to_pandas()
    out = []
    for r in g.itertuples():
        if r.partial5 or r.cds_codons * 3 < 300:
            continue
        s = seqs[r.contig]
        if r.strand == 1:
            if r.begin - 1 < 20: continue
            w = s[r.begin - 21:r.begin - 1]
        else:
            if len(s) - r.end < 20: continue
            w = rc(s[r.end:r.end + 20])
        if len(w) == 20 and set(w) <= set('ACGT'):
            out.append((r.gene_id, rna(w)))
    return out


def duplex(tail, ws):
    """dG array and midpoint (relative to start, -20..-1) of the window-side paired span."""
    dg = np.empty(len(ws)); mid = np.empty(len(ws))
    for k, w in enumerate(ws):
        d = RNA.duplexfold(tail, w)
        dg[k] = d.energy
        s2 = d.structure.split('&')[1] if '&' in d.structure else ''
        if s2:
            mid[k] = (d.j + (d.j + len(s2) - 1)) / 2 - 21
        else:
            mid[k] = np.nan
    return dg, mid


def stat(dg_t, mid_t, dg_can, exclusive=True):
    """(E, n_hit, n_hit_in_position) for one tail."""
    hit = dg_t <= STRONG
    if exclusive:
        hit &= dg_can > CAN_WEAK
    inpos = hit & (mid_t >= POS_LO) & (mid_t <= POS_HI)
    return hit.mean(), int(hit.sum()), int(inpos.sum())


def random_kmers(seqs, n, k, gc_lo, gc_hi, allow_core, rng, frac=False):
    """n random k-mers (RNA) from the genome with GC (count, or fraction if frac) in [gc_lo, gc_hi]."""
    contigs = [(c, s) for c, s in seqs.items() if len(s) >= k + 1]
    lens = np.array([len(s) for _, s in contigs], float)
    p = lens / lens.sum()
    out, tries = [], 0
    while len(out) < n:
        tries += 1
        if tries > 5_000_000:
            raise RuntimeError('cannot draw decoys')
        c, s = contigs[rng.choice(len(contigs), p=p)]
        i = rng.integers(0, len(s) - k)
        w = s[i:i + k]
        if not set(w) <= set('ACGT'):
            continue
        if rng.random() < 0.5:
            w = rc(w)
        w = rna(w)
        g = gc(w) / k if frac else gc(w)
        if not (gc_lo <= g <= gc_hi):
            continue
        if not allow_core and has_core(w):
            continue
        out.append(w)
    return out


def decoy_stats(acc, seqs, ws, dg_can, tail_gc, allow_core, exclusive, n=1000, seed=0, cache_dir='out/cache'):
    """E, n_hit, n_inpos for n GC-matched decoy tails; cached per (acc, gc, core, exclusive, n, seed)."""
    os.makedirs(cache_dir, exist_ok=True)
    f = f'{cache_dir}/{acc}_gc{tail_gc}_core{int(allow_core)}_ex{int(exclusive)}_n{n}_s{seed}.npy'
    if os.path.exists(f):
        return np.load(f)
    rng = np.random.default_rng(seed)
    tails = random_kmers(seqs, n, 13, tail_gc - 1, tail_gc + 1, allow_core, rng)
    wl = [w for _, w in ws]
    res = []
    for t in tails:
        if exclusive:  # only weak-canonical genes can be exclusive
            sel = np.where(dg_can > CAN_WEAK)[0]
            d, m = duplex(t, [wl[i] for i in sel])
            hit = d <= STRONG
            inpos = hit & (m >= POS_LO) & (m <= POS_HI)
            res.append((hit.sum() / len(wl), hit.sum(), inpos.sum()))
        else:
            d, m = duplex(t, wl)
            res.append(stat(d, m, dg_can, exclusive=False))
    res = np.array(res, float)
    np.save(f, res)
    return res


def call(e_obs, n_obs, npos_obs, dec):
    """RULES T1 genome call: p <= 0.01 (>= 99% of decoys below), z >= 3, positional ratio >= 1.15 (A2)."""
    de = dec[:, 0]
    p = (1 + (de >= e_obs).sum()) / (1 + len(de))
    sd = de.std()
    z = (e_obs - de.mean()) / sd if sd > 0 else (np.inf if e_obs > de.mean() else 0.0)
    pos_obs = npos_obs / n_obs if n_obs else 0.0
    pos_dec = dec[:, 2].sum() / dec[:, 1].sum() if dec[:, 1].sum() else np.nan
    ratio = pos_obs / pos_dec if pos_dec and pos_dec > 0 else (np.inf if pos_obs > 0 else np.nan)  # A6
    ok = p <= 0.01 and z >= 3 and n_obs > 0 and ratio >= 1.15  # RULES A2, A6 (gate on the reported p)
    return dict(E=e_obs, n=n_obs, n_pos=npos_obs, p=p, z=z, pos_ratio=ratio, dec_mean=de.mean(), call=bool(ok))
