"""tRNA calls: ARAGORN (pyaragorn) exactly as codon-hunt's genome-wide screen ran it (src/trna_screen.py:
RNAFinder(11, tmrna=False, linear=True); acceptor stem from aragorn_trnas.acceptor, T-loop anchored), so new calls
line up with that screen's 191,958 GTDB species (imports: codon_hunt_trnas).

trna.parquet columns (codon-hunt TRNA_COLS minus acc): contig, begin, end, strand, aa, anticodon, energy,
stem5, stem3, pair3_70 (e.g. "G:T" = G3:U70), N73, CCA, stem_pairs (Watson-Crick/GU pairs of 7), seq
"""
import threading

import pyaragorn

PAIRS = {("G", "C"), ("C", "G"), ("A", "T"), ("T", "A"), ("G", "T"), ("T", "G")}
TAIL_RANK = {"CCA": 0, "": 1, "CC": 2, "C": 3}
# ARAGORN's C code keeps global state: two find_rna calls in parallel threads swapped tRNAs between genomes
# (found in review, 2026-11-11), so calls are serialized. codon-hunt ran it in separate processes.
_LOCK = threading.Lock()
COLS = ["contig", "begin", "end", "strand", "aa", "anticodon", "energy", "stem5", "stem3", "pair3_70", "N73", "CCA",
        "stem_pairs", "seq"]


def acceptor(seq):
    """Acceptor stem 1-7 / 66-72 and discriminator 73 (codon-hunt src/aragorn_trnas.py, unchanged).

    ARAGORN's reported ends vary: CCA encoded or not, plus 0-2 flanking bases at either end. Flanks can pair with
    each other, so pairing alone can pick a frame shifted by one. Anchor the 3' side on the T-loop instead (T54-T55-C56
    sits exactly 12 nt before position 66 in canonical tRNAs), then pick the 5' start that pairs best with 66-72.
    Ties go to the smaller 5' offset, a CCA-like tail, then the shorter trim.
    """
    s = seq.upper()
    best = None
    for o in range(0, 3):
        stem5 = s[o:o + 7]
        for t in range(0, 6):
            e = len(s) - t
            if e - o < 50:
                break
            stem3, tail = s[e - 8:e - 1], s[e:]
            wc = sum((stem5[i], stem3[6 - i]) in PAIRS for i in range(7))
            tloop = s[e - 20:e - 17] == "TTC"
            key = (-(wc + 3 * tloop), o, TAIL_RANK.get(tail, 4), t)
            if best is None or key < best[0]:
                best = (key, stem5, stem3, s[e - 1], tail == "CCA", wc)
    if best is None:  # too short to place an acceptor stem
        return s[:7], "", "?:?", "?", False, 0
    _, stem5, stem3, n73, cca, wc = best
    return stem5, stem3, f"{stem5[2]}:{stem3[4]}", n73, cca, wc


def call(contigs):
    """-> list of row tuples in COLS order."""
    rows = []
    with _LOCK:  # everything read from ARAGORN's results stays inside the lock
        finder = pyaragorn.RNAFinder(11, tmrna=False, linear=True)
        for name, s in contigs:
            for g in finder.find_rna(s):
                if g.type != "tRNA":
                    continue
                seq = g.sequence() if callable(g.sequence) else g.sequence
                rows.append((name, g.begin, g.end, g.strand, g.amino_acid, g.anticodon.upper(), round(g.energy, 1),
                             seq))
    return [r[:7] + acceptor(r[7]) + (r[7],) for r in rows]
