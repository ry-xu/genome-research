"""Kill check 3: reverse view without brown-algal references. Choristocarpus NC_085326 annotated CDS translated with
ATA -> X (TGA -> W), blastp against mito proteins of non-Phaeophyceae stramenopiles (one record per genus, annotated
translations). For each Chor ATA/ATG/ATT/ATC codon (query pos > 10): residue of the best hit per subject genus; site call =
majority residue if >= 3 genera and >= 2/3 agree. Also per lineage group.
usage: uv run python -I checks/chor_wgs/kill/k4_reverse.py"""
import os, re, subprocess
from collections import Counter, defaultdict
import polars as pl
from Bio.Seq import Seq

RD = "scratchpad/chor_work/read"
OUT = "checks/chor_wgs/kill"; tmp = os.path.join(OUT, "tmp"); os.makedirs(tmp, exist_ok=True)
rec = pl.read_parquet("data/records_ranked.parquet").filter((pl.col("kind") == "mito") & pl.col("lineage").str.contains("Stramenopiles")
                                                             & (pl.col("class") != "Phaeophyceae") & pl.col("genus").is_not_null())
rec = rec.sort("acc").group_by("genus", maintain_order=True).first()
grp = {}
for r in rec.iter_rows(named=True):
    L = r["lineage"].split(";")
    grp[r["acc"]] = (r["genus"], "Oomycota" if "Oomycota" in L else "diatoms" if "Bacillariophyta" in L else
                     "Bigyra" if "Bigyra" in L else "other ochrophytes")
print("subject genera", len(grp), "by group", dict(Counter(v[1] for v in grp.values())))
cds = pl.read_parquet("data/cds_mito.parquet", columns=["acc", "gene", "loc", "translation"])
sub = cds.filter(pl.col("acc").is_in(list(grp)) & pl.col("translation").is_not_null())
db = os.path.join(tmp, "nonbrown")
with open(db + ".fa", "w") as fh:
    for i, r in enumerate(sub.iter_rows(named=True)):
        fh.write(f">{r['acc']}|{i}|{r['gene']}\n{r['translation'].rstrip('*')}\n")
print("subject proteins", sub.height)
subprocess.run(["makeblastdb", "-in", db + ".fa", "-dbtype", "prot", "-out", db], check=True, capture_output=True)
g = "".join(l.strip().upper() for l in open(f"{RD}/POS_NC_085326_target.fa") if not l.startswith(">"))
qf = os.path.join(tmp, "chor_X.fa"); qcod = {}
with open(qf, "w") as fh:
    for r in cds.filter(pl.col("acc") == "NC_085326.1").iter_rows(named=True):
        m = re.fullmatch(r"(complement\()?(\d+)\.\.(\d+)\)?", r["loc"])
        nt = g[int(m.group(2)) - 1:int(m.group(3))]
        if m.group(1): nt = str(Seq(nt).reverse_complement())
        cods = [nt[i:i + 3] for i in range(0, len(nt) - 2, 3)]
        if cods and str(Seq(cods[-1]).translate(table=4)) == "*": cods = cods[:-1]
        qcod[r["gene"]] = cods
        fh.write(f">{r['gene']}\n" + "".join("X" if c == "ATA" else str(Seq(c).translate(table=4)) for c in cods) + "\n")
out = subprocess.run(["blastp", "-query", qf, "-db", db, "-evalue", "1e-5", "-num_threads", "1", "-max_target_seqs", "2000",
                      "-outfmt", "6 qseqid sseqid qstart sstart bitscore qseq sseq"], check=True, capture_output=True, text=True).stdout
best = {}   # (qgene, subject acc) -> (bits, line)
for l in out.splitlines():
    q, s, qs, ss, bs, qa, sa = l.split("\t")
    acc = s.split("|")[0]
    if (q, acc) not in best or float(bs) > best[(q, acc)][0]: best[(q, acc)] = (float(bs), l)
votes = defaultdict(dict)  # (gene, qpos) -> {genus: residue}
for (q, acc), (bs, l) in best.items():
    _, s, qs, ss, _, qa, sa = l.split("\t")
    qp = int(qs) - 1
    for a, b in zip(qa, sa):
        if a == "-": continue
        qp += 1
        if b != "-": votes[(q, qp)][acc] = b
def summarise(codon, by=None):
    calls, allv = Counter(), Counter()
    for gene, cods in qcod.items():
        for i, c in enumerate(cods):
            if c != codon or i + 1 <= 10: continue
            v = {a: r for a, r in votes.get((gene, i + 1), {}).items() if by is None or grp[a][1] == by}
            if not v: calls["(no hits)"] += 1; continue
            cnt = Counter(v.values()); allv.update(cnt)
            top, n = cnt.most_common(1)[0]
            calls[top if len(v) >= 3 and n / len(v) >= 2 / 3 else "(mixed/few)"] += 1
    return calls, allv
for by in [None, "Oomycota", "diatoms", "other ochrophytes", "Bigyra"]:
    print(f"\n== subjects: {by or 'all non-brown stramenopiles'}")
    for codon in ["ATA", "ATG", "ATT", "ATC"]:
        calls, allv = summarise(codon, by)
        tot = sum(allv.values())
        print(f"Chor {codon} sites={sum(calls.values())} site calls {dict(calls.most_common(6))} | residue votes M={allv['M']} I={allv['I']} "
              f"V={allv['V']} L={allv['L']} of {tot} (M share {allv['M'] / tot if tot else 0:.3f}, I share {allv['I'] / tot if tot else 0:.3f})")
