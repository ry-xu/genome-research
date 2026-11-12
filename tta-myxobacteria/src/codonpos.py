"""Arm A core: where in called genes does each codon sit?

If sense codon X became a stop, a gene caller using the standard code calls every gene that ends at X as running on
past X to the next TAA/TAG/TGA. X then sits near the 3' end of called genes (at the true end, and in the read-on
tail) far more than in gene interiors. This module turns a set of called genes (nucleotide, stop included) into
per-codon counts by distance from the gene end.

profile(genes) -> dict with
  pos   int64 [65, NB]  codon (0-63, 64 = contains a non-ACGT base) x position bin, sense part of complete genes only
                        bins 0..NB-2 by d = codons from the end (last sense codon d=1), for codons at index >= HEAD;
                        bin NB-1 = first HEAD codons (start region, kept apart)
  g30   int64 [65]      genes with the codon at least once at d <= 30 (index >= HEAD)
  g100  int64 [65]      same, d <= 100
  gany  int64 [65]      genes with the codon anywhere (index >= HEAD)
  stops int64 [3]       complete genes ending TAA, TAG, TGA
  ngenes, ncomplete, gc  (gc over the sense part of complete genes)
  rf    list of (gene_id, protein) for proteins 280-480 aa containing GGQ (release-factor prefilter)
"""
import numpy as np

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CIDX = {c: i for i, c in enumerate(CODONS)}
STOPS = ("TAA", "TAG", "TGA")
STOP_IDX = [CIDX[s] for s in STOPS]
SENSE = [i for i, c in enumerate(CODONS) if c not in STOPS]
HEAD = 20
# bin edges on d (codons from the end, 1-based); last edge open
EDGES = [1, 11, 21, 31, 41, 51, 61, 71, 81, 91, 101, 126, 151, 201, 301]
NB = len(EDGES) + 1  # + start-region bin
D_LUT_MAX = 100000
_dbin = np.searchsorted(np.array(EDGES), np.arange(D_LUT_MAX + 1), side="right") - 1
_dbin[0] = 0
_B = np.full(256, 4, dtype=np.int64)
for i, b in enumerate(BASES):
    _B[ord(b)] = i
    _B[ord(b.lower())] = i
AA_LUT = np.frombuffer((AA + "X").encode(), dtype=np.uint8)


def codon_index(nt_bytes):
    """nt bytes (length multiple of 3) -> codon index array (64 = non-ACGT)."""
    b = _B[np.frombuffer(nt_bytes, dtype=np.uint8)]
    b = b.reshape(-1, 3)
    bad = (b == 4).any(axis=1)
    idx = b[:, 0] * 16 + b[:, 1] * 4 + b[:, 2]
    idx[bad] = 64
    return idx


def profile(genes, rf_min=280, rf_max=480, stops_used=STOPS):
    """genes: iterable of (gene_id, nt_str, complete_bool). nt includes the stop codon. stops_used: the codons that
    end genes under the table the genes were called with (the `stops` count array stays TAA, TAG, TGA)."""
    seqs, ids = [], []
    ngenes = ncomp = 0
    stops = np.zeros(3, dtype=np.int64)
    for gid, nt, complete in genes:
        ngenes += 1
        if not complete or len(nt) % 3 or len(nt) < 3 * (HEAD + 31):
            continue
        last = nt[-3:].upper()
        if last not in stops_used:
            continue
        ncomp += 1
        if last in STOPS:
            stops[STOPS.index(last)] += 1
        seqs.append(nt[:-3])
        ids.append(gid)
    out = {"pos": np.zeros((65, NB), dtype=np.int64), "g30": np.zeros(65, dtype=np.int64),
           "g100": np.zeros(65, dtype=np.int64), "gany": np.zeros(65, dtype=np.int64), "stops": stops,
           "ngenes": ngenes, "ncomplete": ncomp, "gc": float("nan"), "rf": []}
    if not seqs:
        return out
    cat = "".join(seqs).encode("ascii", "replace")
    cidx = codon_index(cat)
    lens = np.array([len(s) // 3 for s in seqs], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(lens)[:-1]])
    gene_of = np.repeat(np.arange(len(lens)), lens)
    i_in = np.arange(len(cidx)) - starts[gene_of]
    d = lens[gene_of] - i_in
    pbin = _dbin[np.minimum(d, D_LUT_MAX)]
    pbin[i_in < HEAD] = NB - 1
    out["pos"] = np.bincount(cidx * NB + pbin, minlength=65 * NB).reshape(65, NB)
    body = i_in >= HEAD
    for key, lim in (("g30", 30), ("g100", 100), ("gany", None)):
        m = body if lim is None else body & (d <= lim)
        pair = np.unique(gene_of[m] * 65 + cidx[m])
        out[key] = np.bincount(pair % 65, minlength=65)
    arr = np.frombuffer(cat, dtype=np.uint8)
    gcn = np.isin(arr, np.frombuffer(b"GCgc", dtype=np.uint8)).sum()
    out["gc"] = float(gcn) / max(1, len(arr))
    # release-factor prefilter: translate only length-matched genes
    aa = AA_LUT[cidx]
    for k in np.nonzero((lens >= rf_min) & (lens <= rf_max))[0]:
        p = aa[starts[k]:starts[k] + lens[k]].tobytes()
        if b"GGQ" in p:
            out["rf"].append((ids[k], p.decode()))
    return out


def parse_fasta_text(text):
    name, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            name, buf = line[1:], []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def prodigal_genes(fasta_text):
    """Prodigal-format headers (GTDB): '>id # start # end # strand # ID=..;partial=00;...'."""
    for h, s in parse_fasta_text(fasta_text):
        part = h.split("partial=")[1][:2] if "partial=" in h else "00"
        yield h.split()[0], s, part == "00"


TABLE_STOPS = {11: ("TAA", "TAG", "TGA"), 4: ("TAA", "TAG"), 25: ("TAA", "TAG"), 16: ("TAA", "TGA"),
               15: ("TAA", "TGA"), 22: ("TCA", "TAA", "TGA"), 23: ("TTA", "TAA", "TAG", "TGA")}
