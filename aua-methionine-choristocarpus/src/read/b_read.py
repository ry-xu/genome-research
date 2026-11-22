"""Read ATN codons at reference Met/Ile positions by tblastn of Ile-code brown-algal mito proteins (RULES.md + amendments).
usage: uv run python -I checks/chor_wgs/b_read.py TARGET.fa LABEL WORKDIR [DROP_ORDER]   -> printed counts, WORKDIR/LABEL_sites.tsv"""
import gzip, os, subprocess, sys
from collections import Counter, defaultdict
import polars as pl

REFS = {"Cutleriales": "NC_085327.1", "Desmarestiales": "NC_007684.1", "Dictyotales": "NC_007685.1", "Ectocarpales": "NC_003055.1",
        "Fucales": "NC_007683.1", "Laminariales": "NC_004024.1", "Scytothamnales": "NC_085336.1", "Sporochnales": "NC_085337.1"}
target, label, wd = sys.argv[1:4]
drop = sys.argv[4] if len(sys.argv) > 4 else None
os.makedirs(wd, exist_ok=True)
seqs, name = {}, None
op = gzip.open if target.endswith(".gz") else open
for l in op(target, "rt"):
    if l.startswith(">"):
        name = l[1:].split()[0]; seqs[name] = []
    else:
        seqs[name].append(l.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
tf = os.path.join(wd, label + "_target.fa")
with open(tf, "w") as fh:
    for k, v in seqs.items():
        fh.write(f">{k}\n{v}\n")
cds = pl.read_parquet("data/cds_mito.parquet", columns=["acc", "gene", "translation"]).filter(pl.col("translation").is_not_null())
COMP = str.maketrans("ACGT", "TGCA")
# per (contig, strand, codon start) -> {ref acc: (bits, residue, query pos)}
best = defaultdict(dict)
for order, acc in REFS.items():
    if order == drop:
        continue
    q = os.path.join(wd, f"q_{acc}.fa")
    with open(q, "w") as fh:
        for i, r in enumerate(cds.filter(pl.col("acc") == acc).iter_rows(named=True)):
            fh.write(f">{i}_{r['gene']}\n{r['translation'].rstrip('*')}\n")
    out = subprocess.run(["tblastn", "-query", q, "-subject", tf, "-db_gencode", "4", "-seg", "no", "-evalue", "1e-5",
                          "-outfmt", "6 qseqid sseqid qstart qend sstart send bitscore qseq sseq"],
                         check=True, capture_output=True, text=True).stdout
    for l in out.splitlines():
        qid, sid, qs, qe, ss, se, bs, qa, sa = l.split("\t")
        qs, ss, se, bs = int(qs), int(ss), int(se), float(bs)
        fwd = se > ss
        qp, k = qs, 0                                   # k = codons consumed in subject
        for a, b in zip(qa, sa):
            if a != "-" and b != "-" and qp > 10:
                p = ss + 3 * k if fwd else ss - 3 * k       # first base (1-based, on its strand: fwd start or rev highest)
                key = (sid, fwd, p)
                old = best[key].get(acc)
                if old is None or bs > old[0]:
                    best[key][acc] = (bs, a, qid)
            if b != "-":
                k += 1
            if a != "-":
                qp += 1
rows, C = [], Counter()
for (sid, fwd, p), d in best.items():
    s = seqs[sid]
    cod = s[p - 1:p + 2] if fwd else s[p - 3:p][::-1].translate(COMP)
    if not cod.startswith("AT"):
        continue
    res = Counter(v[1] for v in d.values())
    n = len(d)
    top, c = res.most_common(1)[0]
    if n >= 2 and c / n >= 2 / 3 and top in "MI":
        C[(cod, top)] += 1
        rows.append((sid, "+" if fwd else "-", p, cod, top, n, c, ";".join(sorted({v[2] for v in d.values()}))))
pl.DataFrame(rows, schema=["contig", "strand", "pos", "codon", "ref_res", "n_refs", "n_agree", "genes"], orient="row") \
  .write_csv(os.path.join(wd, label + "_sites.tsv"), separator="\t")
print(label, "drop", drop, " ".join(f"{c}:M={C[(c, 'M')]}/I={C[(c, 'I')]}" for c in ["ATA", "ATG", "ATT", "ATC"]))
ata_m, ata_i, att_m, att_i = C[("ATA", "M")], C[("ATA", "I")], C[("ATT", "M")], C[("ATT", "I")]
call = ata_m >= 10 and ata_m >= 2 * ata_i and att_m <= 0.1 * (att_m + att_i)
print(label, "CALL ATA=Met" if call else "no call")
