"""Synthetic sense-to-stop genomes, through the production path (Arm A recall and false-positive check).

For a real genome: call genes (pyrodigal, single mode, table 11, as GTDB's Prodigal does); then
  1. in-frame interior uses of codon X (not the start codon) are replaced, each with probability q, by a synonymous
     codon drawn from the genome's own synonymous usage (Trp TGG: by Phe/Tyr codons);
  2. the stop codon of each complete gene is rewritten to X with probability s (downstream DNA untouched, so the
     re-called gene runs on into the genome's real sequence to the next TAA/TAG/TGA);
  codons that overlap another gene's start or stop codon are left alone in both steps (audit 2026-11-10: an ATGA
  overlap or a neighbour's stop would otherwise change too);
  3. optionally fragmented into MAG-like contigs (log-normal lengths, median 15 kb) with a share of them dropped;
  4. genes re-called from scratch (pyrodigal trained on the modified genome) and profiled (codonpos.profile).
s = 0, q = 0 is the null run (same re-calling path, no recoding).

usage: python src/inject.py <fna.gz> <out.npz> X:s:q:frag:seed [X:s:q:frag:seed ...]
       frag = 0 (as is) or the kept fraction of fragments (e.g. 0.75)
"""
import gzip
import random
import sys

import numpy as np
import pyrodigal

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import codonpos as cp

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(fn):
    name, buf = None, []
    with gzip.open(fn, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf).upper()


def call(contigs):
    """As GTDB: Prodigal single mode, table 11 unless table 4 raises coding density by > 5 points to >= 70%."""
    total = sum(len(s) for _, s in contigs)
    best = None
    for table in (11, 4):
        finder = pyrodigal.GeneFinder(meta=total < 100000)
        if total >= 100000:
            finder.train(*[s.encode() for _, s in contigs], translation_table=table)
        out, cod = [], 0
        for name, seq in contigs:
            for g in finder.find_genes(seq.encode()):
                out.append((name, g.begin, g.end, g.strand, g.partial_begin, g.partial_end))
                cod += g.end - g.begin + 1
        dens = cod / max(1, total)
        if best is None or (table == 4 and dens > best[1] + 0.05 and dens >= 0.70):
            best = (out, dens, table)
    return best[0]


def genes_nt(contigs, calls):
    seqs = dict(contigs)
    for name, b, e, st, pb, pe in calls:
        nt = seqs[name][b - 1:e]
        if st == -1:
            nt = rc(nt)
        yield f"{name}_{b}_{e}_{st}", nt, not (pb or pe)


def synonyms(contigs, calls):
    use = np.zeros(64)
    for _, nt, _ in genes_nt(contigs, calls):
        ci = cp.codon_index(nt[: len(nt) // 3 * 3].encode())
        use += np.bincount(ci, minlength=65)[:64]
    syn = {}
    for c in cp.SENSE:
        aa = cp.AA[c]
        sib = [k for k in cp.SENSE if cp.AA[k] == aa and k != c]
        if aa == "W":
            sib = [k for k in cp.SENSE if cp.AA[k] in "FY"]
        if aa == "M":
            sib = [k for k in cp.SENSE if cp.AA[k] in "IL"]
        w = np.array([use[k] + 1 for k in sib])
        syn[cp.CODONS[c]] = ([cp.CODONS[k] for k in sib], w / w.sum())
    return syn


def ends_cover(calls):
    """-> {contig: {position: number of gene start/stop codons covering it}} (0-based positions)."""
    cov = {}
    for name, b, e, st, pb, pe in calls:
        c = cov.setdefault(name, {})
        left, right = (b - 1, e - 3)  # codons at the two ends of the call
        lp, rp = (pb, pe)
        for p0, part in ((left, lp), (right, rp)):
            if not part:
                for x in range(p0, p0 + 3):
                    c[x] = c.get(x, 0) + 1
    return cov


def recode(contigs, calls, X, s, q, rng, syn):
    seqs = {n: bytearray(sq.encode()) for n, sq in contigs}
    sibs, w = syn[X]
    cov = ends_cover(calls)
    n_int = n_stop = 0
    for name, b, e, st, pb, pe in calls:
        arr = seqs[name]
        # gene coordinates in its own orientation: codon k covers genome positions
        L = e - b + 1
        ncod = L // 3
        last = ncod - 1 if not (pe if st == 1 else pb) else ncod  # exclude the stop codon if complete
        c = cov.get(name, {})
        for k in range(1, last):
            if st == 1:
                p0 = b - 1 + 3 * k
                cod = arr[p0:p0 + 3].decode()
            else:
                p0 = e - 3 * (k + 1)
                cod = rc(arr[p0:p0 + 3].decode())
            if cod == X and not any(x in c for x in range(p0, p0 + 3)) and rng.random() < q:
                new = sibs[rng.choices(range(len(sibs)), weights=w)[0]]
                arr[p0:p0 + 3] = (new if st == 1 else rc(new)).encode()
                n_int += 1
    for name, b, e, st, pb, pe in calls:
        if (pe if st == 1 else pb) or rng.random() >= s:
            continue
        arr = seqs[name]
        c = cov.get(name, {})
        if st == 1:
            p0 = e - 3
            if any(c.get(x, 0) > 1 for x in range(p0, p0 + 3)):
                continue
            if arr[p0:p0 + 3].decode() in cp.STOPS:
                arr[p0:p0 + 3] = X.encode()
                n_stop += 1
        else:
            p0 = b - 1
            if any(c.get(x, 0) > 1 for x in range(p0, p0 + 3)):
                continue
            if rc(arr[p0:p0 + 3].decode()) in cp.STOPS:
                arr[p0:p0 + 3] = rc(X).encode()
                n_stop += 1
    return [(n, seqs[n].decode()) for n, _ in contigs], n_int, n_stop


def fragment(contigs, keep, rng):
    out = []
    for n, sq in contigs:
        i = 0
        while i < len(sq):
            L = int(min(200000, max(2000, rng.lognormvariate(np.log(15000), 0.8))))
            if rng.random() < keep:
                out.append((f"{n}_{i}", sq[i:i + L]))
            i += L
    return out


def main(fn, out, specs):
    contigs = list(read_fasta(fn))
    base_calls = call(contigs)
    syn = synonyms(contigs, base_calls)
    res = {}
    for spec in specs:
        X, s, q, frag, seed = spec.split(":")
        s, q, frag, seed = float(s), float(q), float(frag), int(seed)
        rng = random.Random(seed)
        if s > 0 or q > 0:
            mod, n_int, n_stop = recode(contigs, base_calls, X, s, q, rng, syn)
        else:
            mod, n_int, n_stop = contigs, 0, 0
        if frag > 0:
            mod = fragment(mod, frag, rng)
        prof = cp.profile(genes_nt(mod, call(mod)))
        res[spec] = dict(pos=prof["pos"], gc=prof["gc"], ncomplete=prof["ncomplete"], n_int=n_int, n_stop=n_stop,
                         g100=prof["g100"], stops=prof["stops"])
        print(fn, spec, prof["ncomplete"], n_int, n_stop, flush=True)
    np.savez_compressed(out, specs=np.array(list(res)), **{f"{k}|{i}": v[k] for i, (_, v) in enumerate(res.items())
                                                           for k in ("pos", "g100", "stops")},
                        meta=np.array([[v["gc"], v["ncomplete"], v["n_int"], v["n_stop"]] for v in res.values()]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
