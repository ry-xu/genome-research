"""Full-genome Codetta-equivalent confirmation (six-frame translation vs all Pfam --enone profiles).

Follows Codetta v2 (Shulgina & Eddy 2021/2022): sequences split into pieces <100 kb; each piece six-frame
translated with the standard code, stops as X; Pfam-A_enone profiles aligned (here with pyhmmer hmmsearch, with
Z/domZ set so i-Evalues approximate hmmscan's); domains with i-Evalue < 1e-10; aligned consensus columns with
posterior '*' (>= 0.95); excluded Pfam classes (mito, viral, seleno, Pyl, transposon, bad); per-column cap of 1%
of a codon's columns (or 1 per unique column if the codon has < 100 columns); inference as in Codetta.

Also downloads the genome's NCBI annotation (GFF) when available and reports tRNA genes with anticodons.

usage: python3 confirm_genome.py <GCA/GCF accession> [cpus]
writes results/confirm/<acc>.{json,codons.tsv}
"""
import io, json, sys, time, zipfile, re
from pathlib import Path
import numpy as np
import requests
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, TextSequenceBlock
from pyhmmer.plan7 import HMMFile, Pipeline, Background

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codetta_lite as cl

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "tools" / "codetta" / "resources"
OUT = ROOT / "results" / "confirm"
OUT.mkdir(parents=True, exist_ok=True)
GEN = ROOT / "data" / "genomes"
GEN.mkdir(parents=True, exist_ok=True)

COMP = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")


def fetch_genome(acc):
    acc = acc.replace("GB_", "").replace("RS_", "")
    fna = GEN / f"{acc}.fna"
    gff = GEN / f"{acc}.gff"
    if fna.exists():
        return fna, (gff if gff.exists() else None)
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
           f"?include_annotation_type=GENOME_FASTA,GENOME_GFF")
    for attempt in range(4):
        r = requests.get(url, timeout=300)
        if r.status_code == 200:
            break
        time.sleep(3 * (attempt + 1))
    r.raise_for_status()
    z = zipfile.ZipFile(io.BytesIO(r.content))
    names = z.namelist()
    fn = [n for n in names if n.endswith(".fna")]
    gf = [n for n in names if n.endswith(".gff")]
    if not fn:
        raise RuntimeError(f"no fna in package for {acc}: {names}")
    fna.write_bytes(z.read(fn[0]))
    if gf:
        gff.write_bytes(z.read(gf[0]))
    return fna, (gff if gf else None)


def read_fasta(path):
    seqs, name, buf = [], None, []
    for line in open(path, "rb"):
        if line.startswith(b">"):
            if name is not None:
                seqs.append((name, b"".join(buf).upper()))
            name, buf = line[1:].split()[0].decode(), []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs.append((name, b"".join(buf).upper()))
    return seqs


def pieces(seqs, maxlen=100000):
    out = []
    for name, s in seqs:
        parts = [s]
        while any(len(p) > maxlen for p in parts):
            new = []
            for p in parts:
                if len(p) > maxlen:
                    new += [p[: len(p) // 2], p[len(p) // 2:]]
                else:
                    new.append(p)
            parts = new
        out += parts
    return out


def six_frames(piece_list):
    """Yield (name, codon_index_array, translation) for each piece and frame."""
    for pi, p in enumerate(piece_list):
        rc = p.translate(COMP)[::-1]
        for strand, s in (("+", p), ("-", rc)):
            for f in range(3):
                c = cl.codon_indices(s[f:])
                if len(c) < 30:
                    continue
                yield f"p{pi}{strand}{f}", c, cl.translate(c)


def excluded_pfams():
    ex = set()
    for fn in ["mito_pfams.txt", "viral_pfams.txt", "selenocysteine_pfams.txt", "pyrrolysine_pfams.txt",
               "transposon_pfams.txt", "bad_pfams.txt"]:
        p = RES / fn
        if p.exists():
            for l in open(p):
                l = l.strip()
                if l:
                    ex.add(l)
    return ex


def run(acc, cpus=4):
    t0 = time.time()
    fna, gff = fetch_genome(acc)
    seqs = read_fasta(fna)
    gsize = sum(len(s) for _, s in seqs)
    gc = sum(s.count(b"G") + s.count(b"C") for _, s in seqs) / max(gsize, 1)
    pl = pieces(seqs)
    frames = list(six_frames(pl))
    abc = Alphabet.amino()
    block = TextSequenceBlock([TextSequence(name=n.encode(), sequence=t) for n, c, t in frames]).digitize(abc)
    cod = {n: c for n, c, t in frames}
    ex = excluded_pfams()
    print(f"{acc}: {gsize/1e6:.2f} Mb, GC {gc:.3f}, {len(pl)} pieces, {sum(len(t) for _,_,t in frames)/1e6:.1f} Maa",
          flush=True)
    # collect (pfam, column, codon)
    pf_l, col_l, z_l = [], [], []
    emis = {}
    nprof = 0
    with HMMFile(str(RES / "Pfam-A_enone.hmm")) as hf:
        hmms = [h for h in hf]
    nprof = len(hmms)
    print(f"loaded {nprof} profiles in {time.time()-t0:.0f}s", flush=True)
    keep = []
    for h in hmms:
        nm = h.name.decode() if isinstance(h.name, bytes) else h.name
        if nm in ex:
            continue
        keep.append(h)
    print(f"{len(keep)} profiles after exclusions", flush=True)
    n_dom = 0
    for top in pyhmmer.hmmsearch(keep, block, cpus=cpus, Z=nprof, domZ=100, E=1e-5, domE=1e-10, incE=1e-5,
                                  incdomE=1e-10):
        h = top.query
        hname = h.accession.decode() if isinstance(h.accession, bytes) else h.accession
        hname = hname or (h.name.decode() if isinstance(h.name, bytes) else h.name)
        emitted = False
        for hit in top:
            if not hit.included:
                continue
            c = cod[hit.name.decode() if isinstance(hit.name, bytes) else hit.name]
            for dom in hit.domains:
                if not dom.included or dom.i_evalue >= 1e-10:
                    continue
                ali = dom.alignment
                hs = np.frombuffer(ali.hmm_sequence.encode(), dtype=np.uint8)
                ts = np.frombuffer(ali.target_sequence.encode(), dtype=np.uint8)
                pp = np.frombuffer(ali.posterior_probabilities.encode(), dtype=np.uint8)
                is_m = hs != ord(".")
                is_r = ts != ord("-")
                kk = ali.hmm_from + np.cumsum(is_m) - 1
                ii = ali.target_from - 1 + np.cumsum(is_r) - 1
                sel = is_m & is_r & (pp == ord("*"))
                kk, ii = kk[sel], ii[sel]
                zz = c[ii]
                ok = zz >= 0
                if not ok.any():
                    continue
                if not emitted:
                    em = np.asarray(h.match_emissions)[1:, :20].astype(np.float64)
                    emis[hname] = em / em.sum(1, keepdims=True)
                    emitted = True
                n_dom += 1
                pf_l.append(np.full(ok.sum(), len(emis) - 1, dtype=np.int32))
                col_l.append((kk[ok] - 1).astype(np.int32))
                z_l.append(zz[ok].astype(np.int16))
    names = list(emis.keys())
    print(f"{n_dom} domains from {len(names)} Pfams in {time.time()-t0:.0f}s", flush=True)
    pf = np.concatenate(pf_l); col = np.concatenate(col_l); z = np.concatenate(z_l)
    # emission rows for each aligned column
    offs = np.cumsum([0] + [emis[n].shape[0] for n in names])
    E = np.concatenate([emis[n] for n in names])
    gcol = offs[pf] + col
    Pbar = E[gcol].mean(0)
    LL = np.zeros((64, 20)); N = np.zeros(64, dtype=int)
    for zi in range(64):
        g = gcol[z == zi]
        if len(g) == 0:
            continue
        # per-column cap: 1% of total, or 1 per unique column if < 100 columns
        u, cnt = np.unique(g, return_counts=True)
        if len(g) < 100:
            g = u
        else:
            cap = max(1, int(0.01 * len(g)))
            g = np.repeat(u, np.minimum(cnt, cap))
        N[zi] = len(g)
        LL[zi] = np.log(E[g]).sum(0) - len(g) * np.log(Pbar)
    LLall = np.concatenate([LL, np.zeros((64, 1))], 1)
    post = np.exp(LLall - LLall.max(1, keepdims=True)); post /= post.sum(1, keepdims=True)
    best = post.argmax(1)
    rows = []
    for zi in range(64):
        b = best[zi]; p = post[zi, b]
        call = (cl.AA + "?")[b] if p >= 0.9999 else "?"
        std = cl.STD[zi]
        llstd = LL[zi, cl.AA.index(std)] if std != "*" else 0.0
        top2 = np.argsort(-LLall[zi])[:3]
        rows.append(dict(codon=cl.CODONS[zi], std=std, call=call, post=float(p), n=int(N[zi]),
                         margin_vs_std=float(LLall[zi, b] - llstd),
                         top3=",".join(f"{(cl.AA+'?')[t]}:{LLall[zi,t]:.1f}" for t in top2)))
    code = "".join(r["call"] for r in rows)
    # tRNAs from GFF
    trnas = []
    if gff:
        for line in open(gff):
            if "\ttRNA\t" in line or ";gbkey=tRNA" in line and "\ttRNA" in line:
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "tRNA":
                    continue
                attr = f[8]
                prod = re.search(r"product=([^;]+)", attr)
                ac = re.search(r"anticodon=([^;]+)", attr) or re.search(r"seq:([acgtuACGTU]{3})", attr)
                note = re.search(r"Note=([^;]+)", attr)
                trnas.append(dict(contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                                  product=prod.group(1) if prod else "", anticodon=ac.group(1) if ac else "",
                                  note=note.group(1) if note else "", attr=attr[:300]))
    out = dict(acc=acc, size=gsize, gc=gc, n_domains=n_dom, n_pfams=len(names), code=code,
               nonstandard=[r for r in rows if r["call"] not in ("?", r["std"])], trnas=trnas,
               seconds=time.time() - t0)
    (OUT / f"{acc}.json").write_text(json.dumps(out, indent=1))
    import pandas as pd
    pd.DataFrame(rows).to_csv(OUT / f"{acc}.codons.tsv", sep="\t", index=False)
    print("code:", code)
    print("std :", cl.STD)
    for r in out["nonstandard"]:
        print("NONSTANDARD", r)
    print(f"{len(trnas)} tRNAs in GFF; done in {time.time()-t0:.0f}s")
    return out


if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 3)
