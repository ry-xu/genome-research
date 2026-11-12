"""Pick NCBI assemblies (not in GTDB R232) to screen: under Polyangiaceae/Polyangiales/Polyangia (NCBI lineage) or with an
uninformative Myxococcota-level name. Writes results/catalogues/ncbi_to_screen.tsv with a priority column:
 1 = released after the R232 cutoff (2025-09-04) and NCBI lineage within Polyangia (Polyangiaceae/Polyangiales/Polyangia, incl. "... bacterium")
 2 = released after the cutoff, uninformative higher-level name (Myxococcota/Myxococcales/Deltaproteobacteria bacterium)
 3 = released on/before the cutoff but absent from GTDB R232 (priority-1/2 rules)
usage: uv run --no-project python -I src/catalogues/ncbi_select.py
"""
import csv, json, time, urllib.request
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
F = ROOT / "results/catalogues/ncbi_myxococcota_assemblies.tsv"
CUT = "2025-09-04"
rows = [r for r in csv.DictReader(open(F), delimiter="\t") if r["in_gtdb_r232"] == "0"]
tids = sorted({r["taxid"] for r in rows})
lin = {}
for k in range(0, len(tids), 100):
    url = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/" + ",".join(tids[k:k + 100])
    d = json.load(urllib.request.urlopen(url, timeout=120))
    for t in d.get("taxonomy_nodes", []):
        tx = t.get("taxonomy", {})
        lin[str(tx.get("tax_id"))] = set(tx.get("lineage", [])) | {tx.get("tax_id")}
    time.sleep(0.5)
POLYANGIA = 3031711
UNINF = {2818507, 2026763, 2026735}  # Myxococcota bacterium, Myxococcales bacterium, Deltaproteobacteria bacterium
UNC = {"uncultured Myxococcales bacterium", "uncultured Myxococcota bacterium", "uncultured Deltaproteobacteria bacterium",
       "uncultured delta proteobacterium"}
# drop RefSeq copies whose GenBank twin is also listed
gca = {r["accession"].split(".")[0] for r in rows if r["accession"].startswith("GCA_")}
out = []
for r in rows:
    a = r["accession"]
    if a.startswith("GCF_") and r["paired"] and r["paired"].split(".")[0] in gca:
        continue
    L = lin.get(r["taxid"], set())
    post = r["release_date"] > CUT
    if POLYANGIA in L:
        p = 1 if post else 3
    elif int(r["taxid"]) in UNINF or r["organism"] in UNC:
        p = 2 if post else 3
    else:
        continue
    out.append(dict(priority=p, **r))
out.sort(key=lambda r: (r["priority"], r["release_date"]))
with open(ROOT / "results/catalogues/ncbi_to_screen.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
from collections import Counter
print("not in GTDB:", len(rows), "selected:", len(out), Counter(r["priority"] for r in out))
print("lineage missing for", sum(1 for r in rows if r["taxid"] not in lin))
