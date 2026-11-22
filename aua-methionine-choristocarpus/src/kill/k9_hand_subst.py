"""Kill check 1c + 6: (a) hand-slice codons listed in b_read's *_sites.tsv straight from the target FASTA (plain string code,
own reverse complement), incl. minus-strand sites and negative-control ATG@M sites, with Desmarestia protein context;
(b) is "ATA = Ile plus Met->Ile substitutions" compatible with Chor codon usage?
usage: uv run python -I checks/chor_wgs/kill/k9_hand_subst.py"""
import random
from collections import Counter
import polars as pl
from scipy.stats import binomtest

RD = "scratchpad/chor_work/read"
cds = pl.read_parquet("data/cds_mito.parquet", columns=["acc", "gene", "translation"]).filter((pl.col("acc") == "NC_007684.1") & pl.col("translation").is_not_null())
desm = {f"{i}_{r['gene']}": r["translation"] for i, r in enumerate(cds.iter_rows(named=True))}
SA = pl.read_parquet("checks/chor_wgs/kill/sites_all.parquet")
RC = {"A": "T", "C": "G", "G": "C", "T": "A"}
def fasta(label):
    return "".join(l.strip().upper() for l in open(f"{RD}/{label}_target.fa") if not l.startswith(">"))
random.seed(1)
for label, codon, res in [("cten_WGS", "ATA", "M"), ("POS_NC_085326", "ATA", "M"), ("NEG_Saccharina", "ATG", "M"), ("NEG_Petalonia", "ATG", "M")]:
    g = fasta(label)
    s = pl.read_csv(f"{RD}/{label}_sites.tsv", separator="\t").filter((pl.col("codon") == codon) & (pl.col("ref_res") == res))
    print(f"\n== {label}: {codon}@{res} rows {s.height}, strands {dict(Counter(s['strand'].to_list()))}")
    rows = list(s.iter_rows(named=True))
    pick = [r for r in rows if r["strand"] == "-"][:3] + random.sample([r for r in rows if r["strand"] == "+"], 3)
    for r in pick:
        p = r["pos"]
        if r["strand"] == "+": c = g[p - 1:p + 2]
        else: c = "".join(RC[b] for b in reversed(g[p - 3:p]))
        d = SA.filter((pl.col("target") == label) & (pl.col("pos") == p) & (pl.col("strand") == r["strand"]) & (pl.col("ref") == "NC_007684.1"))
        ctx = ""
        if d.height:
            q, qp = d["qgene"][0], d["qpos"][0]; prot = desm[q]
            ctx = f"Desmarestia {q} pos {qp}: {prot[max(0, qp - 6):qp - 1]}[{prot[qp - 1]}]{prot[qp:qp + 5]}"
        print(f"  {r['strand']} pos {p}: sliced {c} (tsv {r['codon']}) refs {r['n_agree']}/{r['n_refs']} genes {r['genes'][:40]} | {ctx}")
print("\n== (b) substitution model: if ATA were Ile, ATA@M sites are Met->Ile changes; their Ile codon should follow Chor Ile usage")
C = pl.read_parquet("checks/chor_wgs/kill/consensus.parquet")
for t in ["cten_WGS", "dmes_WGS", "NEG_Phaeostrophion", "NEG_Saccharina"]:
    x = C.filter(pl.col("target") == t)
    i = Counter(x.filter(pl.col("cres") == "I")["codon"].to_list()); m = Counter(x.filter(pl.col("cres") == "M")["codon"].to_list())
    ile_tot = i["ATA"] + i["ATT"] + i["ATC"]; p_ata = i["ATA"] / ile_tot
    n_sub = m["ATA"] + m["ATT"] + m["ATC"]
    bt = binomtest(m["ATA"], n_sub, p_ata, alternative="greater")
    print(f"{t}: Ile-column codons ATA/ATT/ATC = {i['ATA']}/{i['ATT']}/{i['ATC']} (ATA share {p_ata:.3f}); M-column Ile-codons ATA/ATT/ATC = "
          f"{m['ATA']}/{m['ATT']}/{m['ATC']}; expected ATA {n_sub * p_ata:.1f} of {n_sub}; binomial P(>= obs) = {bt.pvalue:.3g}; "
          f"M-column Ile-codon share {n_sub / sum(m.values()):.3f}; ATG+ATA at M {m['ATG'] + m['ATA']} of {sum(m.values())} M cols")
