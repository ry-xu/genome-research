"""Kill check 1: independent re-implementation of b_read.py site mapping with a self-consistency test.
For every aligned (query residue, subject codon) pair, extract the codon from the genome with Biopython reverse_complement
and translate it with table 4; it must equal the tblastn sseq letter. Then recount ATN at M/I with the same rules.
Writes kill/sites_all.parquet (one row per target codon x reference: best HSP, gene, query position, query length).
usage: uv run python -I checks/chor_wgs/kill/k1_audit.py"""
import os, subprocess, sys
from collections import Counter, defaultdict
import polars as pl
from Bio.Seq import Seq

RD = "scratchpad/chor_work/read"
OUT = "checks/chor_wgs/kill"
REFS = {"Cutleriales": "NC_085327.1", "Desmarestiales": "NC_007684.1", "Dictyotales": "NC_007685.1", "Ectocarpales": "NC_003055.1",
        "Fucales": "NC_007683.1", "Laminariales": "NC_004024.1", "Scytothamnales": "NC_085336.1", "Sporochnales": "NC_085337.1"}
TARGETS = [("cten_WGS", None), ("dmes_WGS", None), ("POS_NC_085326", None), ("NEG_Phaeostrophion", None),
           ("NEG_Protohalopteris", None), ("NEG_Saccharina", "Laminariales"), ("NEG_Petalonia", "Ectocarpales"),
           ("NEG_Sargassum", "Fucales")]
cds = pl.read_parquet("data/cds_mito.parquet", columns=["acc", "gene", "translation"]).filter(pl.col("translation").is_not_null())
tmp = os.path.join(OUT, "tmp"); os.makedirs(tmp, exist_ok=True)
qlen = {}
for acc in REFS.values():
    with open(os.path.join(tmp, f"q_{acc}.fa"), "w") as fh:
        for i, r in enumerate(cds.filter(pl.col("acc") == acc).iter_rows(named=True)):
            s = r["translation"].rstrip("*"); qlen[(acc, f"{i}_{r['gene']}")] = len(s)
            fh.write(f">{i}_{r['gene']}\n{s}\n")
rows = []
for label, drop in TARGETS:
    tf = os.path.join(RD, label + "_target.fa")
    seqs, name = {}, None
    for l in open(tf):
        if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = []
        else: seqs[name].append(l.strip().upper())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    best = {}   # (sid, strand, p, acc) -> row
    nchk = nbad = 0
    for order, acc in REFS.items():
        if order == drop: continue
        out = subprocess.run(["tblastn", "-query", os.path.join(tmp, f"q_{acc}.fa"), "-subject", tf, "-db_gencode", "4",
                              "-seg", "no", "-evalue", "1e-5", "-outfmt", "6 qseqid sseqid qstart qend sstart send bitscore sframe qseq sseq"],
                             check=True, capture_output=True, text=True).stdout
        for l in out.splitlines():
            qid, sid, qs, qe, ss, se, bs, fr, qa, sa = l.split("\t")
            qs, ss, se, bs, fr = int(qs), int(ss), int(se), float(bs), int(fr)
            minus = fr < 0
            assert minus == (se < ss), l[:80]
            g = seqs[sid]
            qp, k = qs - 1, 0
            for a, b in zip(qa, sa):
                if a != "-": qp += 1
                if b == "-": continue
                if minus:
                    hi = ss - 3 * k               # highest base of this codon (1-based)
                    cod = str(Seq(g[hi - 3:hi]).reverse_complement())
                else:
                    lo = ss + 3 * k
                    cod = g[lo - 1:lo + 2]
                p = hi if minus else lo
                k += 1
                tr = str(Seq(cod).translate(table=4))
                nchk += 1
                if tr != b: nbad += 1
                if a == "-": continue
                key = (sid, "-" if minus else "+", p, acc)
                if key not in best or bs > best[key][-1]:
                    best[key] = (label, sid, "-" if minus else "+", p, cod, acc, qid, qp, qlen[(acc, qid)], a, bs)
    print(label, "aligned subject codons checked", nchk, "translation mismatches", nbad)
    rows += list(best.values())
df = pl.DataFrame(rows, schema=["target", "contig", "strand", "pos", "codon", "ref", "qgene", "qpos", "qlen", "res", "bits"], orient="row")
df.write_parquet(os.path.join(OUT, "sites_all.parquet"))
# recount with rules: qpos > 10, >= 2 refs, >= 2/3 agree, top in M/I
f = df.filter(pl.col("qpos") > 10)
for label, _ in TARGETS:
    C = Counter()
    for (sid, st, p), g in f.filter(pl.col("target") == label).group_by(["contig", "strand", "pos"]):
        cod = g["codon"][0]
        if not cod.startswith("AT"): continue
        res = Counter(g["res"].to_list()); n = g.height
        top, c = res.most_common(1)[0]
        if n >= 2 and c / n >= 2 / 3 and top in ("M", "I"): C[(cod, top)] += 1
    print(label, " ".join(f"{c}:M={C[(c, 'M')]}/I={C[(c, 'I')]}" for c in ["ATA", "ATG", "ATT", "ATC"]))
