"""Arm A score: is a sense codon piled up at the 3' end of called genes?

From codonpos.profile()['pos'] (65 codons x NB bins):
  n3  = count of codon c at d <= 100 (codons from the gene end; the read-on tail of a gene that really ends at c
        lies here, and so does c itself)
  nI  = count of c at d > 150 (interior; d 101-150 is a buffer, the first 20 codons are excluded)
  r   = (n3 / N3) / (nI / NI), with N3, NI the totals over the 61 sense codons in each region
  rho = the null ratio for codon c at this GC (median over ordinary genomes; C-terminal amino-acid composition makes
        it differ from 1, e.g. Lys/Arg/Glu codons sit near C-termini)
  exp3 = nI * (N3 / NI) * rho      expected 3' count if c behaves like ordinary genomes
  excess = n3 - exp3              ~ genes ending at c + c in the read-on tails
  lr = log2((n3 + 1) / (exp3 + 1))
  p   = Poisson upper tail P(N >= n3 | exp3)
  depl = (nI / NI) / median interior rate of c at this GC (interior depletion; 1 = ordinary)
"""
import math

import numpy as np

import codonpos as cp

W3 = slice(0, 10)       # d 1-100
INT = slice(12, 15)     # d >= 151
SENSE = np.array(cp.SENSE)
GC0, NGC = 0.15, 26     # GC bins of 0.025 from 0.15 to 0.80 (same edges as before, extended both ways: the old top bin
                        # pooled GC 0.675-0.754, where AT-rich codon rates vary ~9-fold; audit 2026-11-10)
EXP_FLOOR = 2.0         # p-value expectation floor: a codon absent from interiors (unassigned, e.g. M. capricolum CGG)
                        # otherwise turns a handful of 3' hits into p ~ 1e-300


def gc_bin(gc):
    assert 0.1 < gc < 0.9, f"not a GC fraction: {gc}"
    return int(min(NGC - 1, max(0, (gc - GC0) // 0.025)))


def regions(pos):
    n3 = pos[:64, W3].sum(axis=1).astype(float)
    nI = pos[:64, INT].sum(axis=1).astype(float)
    return n3, nI


def null_table(profiles):
    """profiles: list of (gc, pos). -> rho[NGC, 64], irate[NGC, 64] medians (nan where < 5 genomes), with fallback
    to the nearest filled GC bin."""
    by = [[] for _ in range(NGC)]
    for gc, pos in profiles:
        n3, nI = regions(pos)
        N3, NI = n3[SENSE].sum(), nI[SENSE].sum()
        if N3 < 1000 or NI < 1000:
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            r = (n3 / N3) / (nI / NI)
        by[gc_bin(gc)].append((r, nI / NI))
    rho = np.full((NGC, 64), np.nan)
    irate = np.full((NGC, 64), np.nan)
    for b in range(NGC):
        if len(by[b]) >= 5:
            R = np.array([x[0] for x in by[b]])
            I = np.array([x[1] for x in by[b]])
            rho[b] = np.nanmedian(np.where(np.isfinite(R), R, np.nan), axis=0)
            irate[b] = np.median(I, axis=0)
    filled = [b for b in range(NGC) if np.isfinite(rho[b]).any()]
    for b in range(NGC):
        if b not in filled and filled:
            nb = min(filled, key=lambda x: abs(x - b))
            rho[b], irate[b] = rho[nb], irate[nb]
    return rho, irate


def poisson_sf_log10(k, mu):
    """log10 P(N >= k), N ~ Poisson(mu); normal/Chernoff-free exact sum via logs for moderate k."""
    if k <= 0:
        return 0.0
    if mu <= 0:
        return -300.0
    # sum_{j>=k} exp(-mu) mu^j / j!, computed from the first term with a ratio series
    lt = -mu + k * math.log(mu) - math.lgamma(k + 1)
    s, term, j = 1.0, 1.0, k
    while True:
        j += 1
        term *= mu / j
        s += term
        if term < 1e-12 * s or j > k + 100000:
            break
    return (lt + math.log(s)) / math.log(10)


def score(pos, gc, rho, irate, sense=cp.SENSE):
    """-> list of dicts, one per codon in `sense` (the codons that are sense under the table the genes were called
    with: the 61 of table 11, plus TGA for table 4, TAG for table 16). A real stop that is sense under the calling
    table has no null rho (never interior in ordinary genomes) and is scored with rho = 1."""
    sense = list(sense)
    n3, nI = regions(pos)
    N3, NI = n3[sense].sum(), nI[sense].sum()
    b = gc_bin(gc)
    rows = []
    for c in sense:
        rh = rho[b, c] if np.isfinite(rho[b, c]) else 1.0
        exp3 = nI[c] * (N3 / NI) * rh if NI > 0 else 0.0
        ir = irate[b, c] if np.isfinite(irate[b, c]) and irate[b, c] > 0 else np.nan
        rows.append(dict(codon=cp.CODONS[c], aa=cp.AA[c], n3=int(n3[c]), nI=int(nI[c]), exp3=exp3,
                         excess=n3[c] - exp3, lr=math.log2((n3[c] + 1) / (exp3 + 1)),
                         log10p=poisson_sf_log10(int(n3[c]), max(exp3, EXP_FLOOR)),
                         depl=(nI[c] / NI) / ir if NI > 0 and np.isfinite(ir) else np.nan))
    return rows
