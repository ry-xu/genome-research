"""Kill check 1b: tblastn-free recount. Annotated CDS of the target (GenBank loc) translated with ATA -> X (no Ile/Met
assumed), global-aligned (free end gaps, BLOSUM62) to the same-named gene of each reference; codon read at reference
M/I columns with the b_read rules (ref pos > 10, >= 2 refs, >= 2/3 agree). Also hand-prints codons at sites.
usage: uv run python -I checks/chor_wgs/kill/k2_annot.py"""
import re
from collections import Counter, defaultdict
import polars as pl
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

RD = "scratchpad/chor_work/read"
REFS = {"Cutleriales": "NC_085327.1", "Desmarestiales": "NC_007684.1", "Dictyotales": "NC_007685.1", "Ectocarpales": "NC_003055.1",
        "Fucales": "NC_007683.1", "Laminariales": "NC_004024.1", "Scytothamnales": "NC_085336.1", "Sporochnales": "NC_085337.1"}
TARGETS = [("POS_NC_085326", "NC_085326.1", None), ("NEG_Saccharina", "NC_013473.1", "Laminariales"),
           ("NEG_Petalonia", "NC_025227.1", "Ectocarpales"), ("NEG_Sargassum", "NC_024613.1", "Fucales")]
cds = pl.read_parquet("data/cds_mito.parquet", columns=["acc", "gene", "loc", "translation", "pseudo"])
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -10, -0.5
al.target_end_gap_score = al.query_end_gap_score = 0.0
al.mode = "global"

def seq_of(label):
    return "".join(l.strip().upper() for l in open(f"{RD}/{label}_target.fa") if not l.startswith(">"))

def cds_nt(g, loc):
    m = re.fullmatch(r"(complement\()?<?(\d+)\.\.>?(\d+)\)?", loc)
    if not m: return None
    s = g[int(m.group(2)) - 1:int(m.group(3))]
    return str(Seq(s).reverse_complement()) if m.group(1) else s

refprot = defaultdict(dict)
for acc in REFS.values():
    for r in cds.filter(pl.col("acc") == acc).iter_rows(named=True):
        if r["translation"] and r["gene"] and r["gene"] not in refprot[acc]:
            refprot[acc][r["gene"]] = r["translation"].rstrip("*")
rows = []
for label, acc, drop in TARGETS:
    g = seq_of(label)
    for r in cds.filter(pl.col("acc") == acc).iter_rows(named=True):
        nt = cds_nt(g, r["loc"]) if r["loc"] else None
        if nt is None or not r["gene"]:
            continue
        cods = [nt[i:i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]
        prot = "".join("X" if c == "ATA" else str(Seq(c).translate(table=4)) for c in cods).rstrip("*")
        for order, ra in REFS.items():
            if order == drop or r["gene"] not in refprot[ra]: continue
            rp = refprot[ra][r["gene"]]
            a = al.align(prot, rp)[0]
            for (t0, t1), (q0, q1) in zip(*a.aligned):
                for i in range(t1 - t0):
                    rows.append((label, r["gene"], t0 + i, cods[t0 + i], len(cods), ra, q0 + i + 1, rp[q0 + i]))
df = pl.DataFrame(rows, schema=["target", "gene", "cidx", "codon", "ncod", "ref", "rpos", "rres"], orient="row")
df.write_parquet("checks/chor_wgs/kill/annot_sites.parquet")
f = df.filter(pl.col("rpos") > 10)
for label, _, _ in TARGETS:
    C = Counter()
    for _, gg in f.filter(pl.col("target") == label).group_by(["gene", "cidx"]):
        cod = gg["codon"][0]
        if not cod.startswith("AT"): continue
        res = Counter(gg["rres"].to_list()); n = gg.height
        top, c = res.most_common(1)[0]
        if n >= 2 and c / n >= 2 / 3 and top in ("M", "I"): C[(cod, top)] += 1
    print("annot-global", label, " ".join(f"{c}:M={C[(c, 'M')]}/I={C[(c, 'I')]}" for c in ["ATA", "ATG", "ATT", "ATC"]))
