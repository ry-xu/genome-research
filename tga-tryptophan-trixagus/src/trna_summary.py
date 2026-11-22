"""Summarise results/trna/*.json: genomes whose used codons have no tRNA reader. Lists, per genome, unread codons used
>= 1 per 1000 codons, with assembly size/contigs and tRNA count. Labels come from src/label.py (NCBI lineage).
usage: python src/trna_summary.py -> results/trna_unread.tsv"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
rows = []
for f in sorted((ROOT / "results" / "trna").glob("*.json")):
    d = json.loads(f.read_text())
    big = {c: v for c, v in d["unread"].items() if v["per_1000"] >= 1}
    rows.append((d["acc"], d["bp"], d["contigs"], d["n_trna"], len(d["aa_with_trna"]),
                 " ".join(f"{c}:{v['aa']}:{v['per_1000']}" for c, v in sorted(big.items(), key=lambda x: -x[1]["per_1000"]))))
with open(ROOT / "results" / "trna_unread.tsv", "w") as f:
    f.write("acc\tbp\tcontigs\tn_trna\tn_aa_with_trna\tunread_codons_per1000\n")
    for r in sorted(rows, key=lambda r: (r[2] > 5, -r[3])):
        f.write("\t".join(map(str, r)) + "\n")
        print(*r, sep="\t")
