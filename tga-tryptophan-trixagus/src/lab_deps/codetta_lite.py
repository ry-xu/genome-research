"""Codetta-lite: genetic-code inference from GTDB marker genes.

Model (Shulgina & Eddy 2021, eLife 10:e71402): each codon Z is associated with the set of profile-HMM consensus
columns C_i it aligns to (posterior >= 0.95). Under decoding A, P(C|A) = P(A|C) P(C) / P(A), where P(A|C) is the
column's match emission and P(A) the mean emission over all columns aligned in that genome. The nonspecific model
has P(C|NS) = P(C). So LL(A) = sum_i log P(A|C_i) - N log P(A), and LL(NS) = 0. Uniform prior over 21 models.

Differences from Codetta (for speed): we align only GTDB's single-copy marker genes (bac120/ar53, Prodigal gene
calls) to their own marker HMM (rebuilt --enone from seeds), instead of a six-frame translation against all Pfam.

Per genome we keep sufficient statistics only:
  S[g, z, a]  = sum over columns aligned to codon z of log P(a|C)        (float32)
  N[g, z]     = number of columns aligned to codon z
  P[g, a]     = sum over all aligned columns of P(a|C)                    (for P(A))
  C[g]        = total aligned columns
plus a per-genome, per-stop-codon truncation table (see `stop_stats`).
"""
import io, sys, tarfile, time
from pathlib import Path
import numpy as np
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, TextSequenceBlock
from pyhmmer.plan7 import HMMFile, Pipeline, Background

ROOT = Path(__file__).resolve().parent.parent
AA = "ACDEFGHIKLMNPQRSTVWY"  # Easel amino order
BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
STD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"  # NCBI table 1/11, TCAG order
STOP_IDX = [i for i, a in enumerate(STD) if a == "*"]  # TAA, TAG, TGA

_lut = np.full(256, 4, dtype=np.int8)
for i, b in enumerate(BASES):
    _lut[ord(b)] = i
    _lut[ord(b.lower())] = i


def codon_indices(nt: bytes) -> np.ndarray:
    """Array of codon indices (0..63, or -1 if ambiguous) for a coding sequence."""
    a = _lut[np.frombuffer(nt, dtype=np.uint8)]
    n = len(a) // 3
    a = a[: 3 * n].reshape(n, 3).astype(np.int16)
    idx = a[:, 0] * 16 + a[:, 1] * 4 + a[:, 2]
    idx[(a >= 4).any(axis=1)] = -1
    return idx


_STD_ARR = np.frombuffer(STD.encode(), dtype=np.uint8)


def translate(cidx: np.ndarray) -> str:
    """Standard-code translation with stops (and ambiguous codons) as X."""
    out = _STD_ARR[np.clip(cidx, 0, 63)].copy()
    out[(cidx < 0) | (out == ord("*"))] = ord("X")
    return out.tobytes().decode()


def load_hmms(path=ROOT / "data" / "markers_enone.hmm"):
    hmms = {}
    with HMMFile(str(path)) as f:
        for h in f:
            name = h.name.decode() if isinstance(h.name, bytes) else h.name
            em = np.asarray(h.match_emissions)[1:, :20].astype(np.float64)  # rows k=1..M
            em = em / em.sum(1, keepdims=True)
            hmms[name] = (h, em, np.log(em))
    return hmms


def align_marker(hmm, em, records, min_bits=30.0):
    """records: list of (genome_index, nt_bytes). Returns dict of numpy arrays:
    g, z, k  (one row per confidently aligned codon) and per-record stop/truncation info."""
    abc = Alphabet.amino()
    seqs, cidx_list, gidx = [], [], []
    for j, (g, nt) in enumerate(records):
        c = codon_indices(nt)
        if len(c) and c[-1] >= 0 and STD[c[-1]] == "*":
            term = int(c[-1])
            c = c[:-1]
        else:
            term = -1  # no terminal stop (contig edge / partial gene)
        if len(c) < 30:
            continue
        cidx_list.append((c, term))
        gidx.append(g)
        seqs.append(TextSequence(name=str(len(seqs)).encode(), sequence=translate(c)))
    if not seqs:
        return None
    block = TextSequenceBlock(seqs).digitize(abc)
    pipe = Pipeline(abc, background=Background(abc), T=min_bits, domT=min_bits, Z=len(seqs), domZ=len(seqs))
    hits = pipe.search_hmm(hmm, block)
    G, Z, K = [], [], []
    # per record: terminal stop codon, hmm_to of the alignment reaching furthest, M
    term_codon, hmm_end, hmm_start, seq_cov = [], [], [], []
    rec_g = []
    for hit in hits:
        if not hit.included:
            continue
        j = int(hit.name.decode() if isinstance(hit.name, bytes) else hit.name)
        c, term = cidx_list[j]
        best_end, best_start = 0, 10**9
        for dom in hit.domains.included:
            ali = dom.alignment
            hs = np.frombuffer(ali.hmm_sequence.encode(), dtype=np.uint8)
            ts = np.frombuffer(ali.target_sequence.encode(), dtype=np.uint8)
            pp = np.frombuffer(ali.posterior_probabilities.encode(), dtype=np.uint8)
            is_m = hs != ord(".")
            is_r = ts != ord("-")
            kk = ali.hmm_from + np.cumsum(is_m) - 1  # 1-based match state index
            ii = ali.target_from - 1 + np.cumsum(is_r) - 1  # 0-based residue index
            sel = is_m & is_r & (pp == ord("*"))
            kk = kk[sel]
            ii = ii[sel]
            zz = c[ii]
            ok = zz >= 0
            G.append(np.full(ok.sum(), gidx[j], dtype=np.int32))
            Z.append(zz[ok].astype(np.int8))
            K.append((kk[ok] - 1).astype(np.int32))  # 0-based row into em
            best_end = max(best_end, ali.hmm_to)
            best_start = min(best_start, ali.hmm_from)
        rec_g.append(gidx[j])
        term_codon.append(term)
        hmm_end.append(best_end)
        hmm_start.append(best_start)
    if not G:
        return None
    return dict(
        g=np.concatenate(G), z=np.concatenate(Z), k=np.concatenate(K),
        rec_g=np.array(rec_g, dtype=np.int32), term=np.array(term_codon, dtype=np.int16),
        hmm_end=np.array(hmm_end, dtype=np.int32), hmm_start=np.array(hmm_start, dtype=np.int32),
        M=hmm.M,
    )


class Accumulator:
    def __init__(self, n_genomes):
        self.S = np.zeros((n_genomes, 64, 20), dtype=np.float32)
        self.N = np.zeros((n_genomes, 64), dtype=np.int32)
        self.P = np.zeros((n_genomes, 20), dtype=np.float64)
        self.C = np.zeros(n_genomes, dtype=np.int64)
        # stop stats: for each genome and terminal stop type (TAA,TAG,TGA,none): [n_genes, n_truncated]
        self.T = np.zeros((n_genomes, 4, 2), dtype=np.int32)
        self.markers_seen = np.zeros(n_genomes, dtype=np.int16)

    def add(self, res, em, logem):
        g, z, k = res["g"], res["z"].astype(np.int64), res["k"]
        np.add.at(self.N, (g, z), 1)
        # S[g,z,:] += logem[k,:]
        flat = g.astype(np.int64) * 64 + z
        Sf = self.S.reshape(-1, 20)
        order = np.argsort(flat, kind="stable")
        fs = flat[order]
        lk = logem[k[order]]
        uniq, start = np.unique(fs, return_index=True)
        sums = np.add.reduceat(lk, start, axis=0)
        Sf[uniq] += sums.astype(np.float32)
        Pk = em[k]
        uniq_g, start_g = np.unique(g, return_index=True)  # g sorted? not necessarily
        og = np.argsort(g, kind="stable")
        gs = g[og]
        ug, sg = np.unique(gs, return_index=True)
        self.P[ug] += np.add.reduceat(Pk[og], sg, axis=0)
        np.add.at(self.C, g, 1)
        # truncation: alignment ends > 15 states before the HMM end
        rg, term, he = res["rec_g"], res["term"], res["hmm_end"]
        tcat = np.full(len(term), 3)
        for t, s in enumerate(STOP_IDX):
            tcat[term == s] = t
        trunc = (he < res["M"] - 15).astype(np.int32)
        np.add.at(self.T, (rg, tcat, 0), 1)
        np.add.at(self.T, (rg, tcat, 1), trunc)
        np.add.at(self.markers_seen, rg, 1)

    def save(self, path, accessions):
        np.savez_compressed(path, S=self.S, N=self.N, P=self.P, C=self.C, T=self.T,
                            markers_seen=self.markers_seen, acc=np.array(accessions))


def infer(S, N, P, C, threshold=0.9999):
    """Infer decodings for one genome (or pooled stats). Returns list of (codon, call, posterior, N, LL_best-LL_std)."""
    Pbar = P / max(C, 1)
    LL = S.astype(np.float64) - N[:, None] * np.log(Pbar)[None, :]
    LLall = np.concatenate([LL, np.zeros((64, 1))], axis=1)  # last = nonspecific
    m = LLall.max(1, keepdims=True)
    post = np.exp(LLall - m)
    post /= post.sum(1, keepdims=True)
    best = post.argmax(1)
    out = []
    for z in range(64):
        b = best[z]
        p = post[z, b]
        call = (AA + "?")[b] if p >= threshold else "?"
        std = STD[z]
        d_std = LL[z, AA.index(std)] if std != "*" else 0.0
        out.append((CODONS[z], call, float(p), int(N[z]), float(LLall[z, b] - (d_std if std != "*" else 0.0))))
    return out


class GrowingFile(io.RawIOBase):
    """File reader that waits for more bytes while a download is still running (until `done_flag` exists)."""
    def __init__(self, path, done_flag):
        self.f = open(path, "rb"); self.done_flag = Path(done_flag)
    def readable(self):
        return True
    def readinto(self, b):
        while True:
            n = self.f.readinto(b)
            if n:
                return n
            if self.done_flag.exists():
                n = self.f.readinto(b)
                return n or 0
            time.sleep(2)


def iter_tar_markers(tar_path, kind="fna", done_flag=None):
    """Stream a GTDB marker tarball; yield (marker_id, seq_bytes). Stops after the `kind` section."""
    fobj = io.BufferedReader(GrowingFile(tar_path, done_flag), 1 << 20) if done_flag else open(tar_path, "rb")
    seen = False
    with tarfile.open(fileobj=fobj, mode="r|gz") as tf:
        for m in tf:
            if seen and m.isfile() and f"/{kind}/" not in m.name:
                break
            if not m.isfile() or f"/{kind}/" not in m.name:
                continue
            seen = True
            mid = Path(m.name).name.split(".")[0]
            data = tf.extractfile(m).read()
            yield mid, data


def parse_fasta_filtered(data: bytes, keep: dict):
    """Return list of (genome_index, seq_bytes) for headers in `keep` (accession -> index)."""
    out = []
    for chunk in data.split(b">")[1:]:
        nl = chunk.find(b"\n")
        acc = chunk[:nl].split()[0].decode()
        gi = keep.get(acc)
        if gi is None:
            continue
        out.append((gi, chunk[nl + 1:].replace(b"\n", b"")))
    return out
