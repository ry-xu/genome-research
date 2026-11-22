"""Build the frozen tier lists (README, 'Frozen before scoring') from the GTDB R232 QC-failed table.
-> results/tier2.txt, results/tier1.txt (T1 minus T2), results/tier3.txt (500 random others, seed 0), results/tiers.tsv"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REASONS = ["many frameshifted proteins", "abnormal gene to sequence ratio", "annotation fails completeness check",
           "RefSeq annotation failed", "missing tRNA genes", "missing rRNA genes", "low gene count"]

q = pd.read_csv(ROOT / "data" / "qc_failed_r232.tsv.gz", sep="\t")
q["acc"] = q["Accession"].str[3:]
c, n50 = q["No. contigs"], q["N50 contigs"]
good = ((c <= 150) & (n50 >= 20000)) | ((c <= 30) & (n50 >= 5000))
cont = (q["CM2 contamination (%)"] <= 5) & (q["CM1 contamination (%)"] <= 5)
low = (q["CM2 completeness (%)"] < 60) | (q["CM1 completeness (%)"] < 60)
ex = q["NCBI exclude from RefSeq"].fillna("")
t2 = ex.apply(lambda s: any(r in s for r in REASONS))
t1 = good & cont & low & ~t2
rest = q.index[~(t1 | t2)]
t3 = pd.Series(False, index=q.index)
t3[np.random.default_rng(0).choice(rest, 500, replace=False)] = True
q["tier"] = np.select([t2, t1, t3], ["T2", "T1", "T3"], "")
out = ROOT / "results"
out.mkdir(exist_ok=True)
for t in ("T1", "T2", "T3"):
    # T1 ordered by CM2 completeness ascending: the most broken good assemblies first
    sub = q[q.tier == t].sort_values("CM2 completeness (%)")
    (out / f"tier{t[1]}.txt").write_text("".join(f"{a}\n" for a in sub.acc))
    print(t, len(sub))
q[q.tier != ""].drop(columns=["Accession"]).to_csv(out / "tiers.tsv", sep="\t", index=False)
