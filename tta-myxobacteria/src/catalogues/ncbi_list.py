"""List NCBI GenBank/RefSeq assemblies for the Polyangiaceae TTA-clade search (NCBI Datasets v2 REST).

Queries every assembly (current versions, GenBank + RefSeq) under:
  Myxococcota (taxid 2818505; includes Polyangia/Polyangiales/Polyangiaceae and the unclassified
  "Myxococcales bacterium", "Myxococcota bacterium", "Polyangia bacterium", "Polyangiales bacterium",
  "Polyangiaceae bacterium") and "Deltaproteobacteria bacterium" (taxid 2026735).
Flags whether each assembly is in GTDB R232 (bac120 metadata, by base accession, GCA or paired GCF).

usage: uv run --no-project python -I src/catalogues/ncbi_list.py
  -> results/catalogues/ncbi_myxococcota_assemblies.tsv
"""
import csv
import gzip
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "catalogues" / "ncbi_myxococcota_assemblies.tsv"
META = ROOT / "data" / "gtdb_meta" / "bac120_metadata.tsv.gz"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{}/dataset_report"
TAXA = {"2818505": "Myxococcota", "2026735": "Deltaproteobacteria bacterium"}


def get(url):
    for i in range(6):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # 429 / transient
            print("retry", i, e, flush=True)
            time.sleep(2 + 3 * i)
    raise RuntimeError(url)


def attrs(bs):
    d = {}
    for a in bs.get("attributes", []) or []:
        d[a.get("name")] = a.get("value")
    return d


def main():
    gtdb = set()
    with gzip.open(META, "rt") as f:
        r = csv.reader(f, delimiter="\t")
        next(r)
        for row in r:
            gtdb.add(row[0][3:].split(".")[0])  # GB_GCA_xxx.1 -> GCA_xxx
    print("GTDB R232 bacterial accessions:", len(gtdb), flush=True)
    rows, seen = [], set()
    for tax, name in TAXA.items():
        tok, n = None, 0
        while True:
            q = {"page_size": 1000, "filters.assembly_source": "all"}
            if tok:
                q["page_token"] = tok
            d = get(API.format(tax) + "?" + urllib.parse.urlencode(q))
            for rep in d.get("reports", []):
                acc = rep["accession"]
                if acc in seen:
                    continue
                seen.add(acc)
                ai, org = rep.get("assembly_info", {}), rep.get("organism", {})
                bs = ai.get("biosample", {}) or {}
                at = attrs(bs)
                ck = rep.get("checkm_info", {}) or {}
                st = rep.get("assembly_stats", {}) or {}
                pair = rep.get("paired_accession", "")
                base = {acc.split(".")[0], pair.split(".")[0] if pair else ""}
                base = {b.replace("GCF_", "GCA_") for b in base if b} | {b for b in base if b}
                rows.append(dict(
                    accession=acc, paired=pair, query_taxon=name, organism=org.get("organism_name", ""),
                    taxid=org.get("tax_id", ""), release_date=ai.get("release_date", ""),
                    assembly_level=ai.get("assembly_level", ""), assembly_name=ai.get("assembly_name", ""),
                    bioproject=ai.get("bioproject_accession", ""), biosample=bs.get("accession", ""),
                    isolation_source=at.get("isolation_source", "") or at.get("env_medium", "") or at.get("host", ""),
                    geo=at.get("geo_loc_name", ""), metagenome_source=at.get("metagenome_source", ""),
                    submitter=ai.get("submitter", ""), total_length=st.get("total_sequence_length", ""),
                    contigs=st.get("number_of_contigs", ""),
                    checkm_completeness=ck.get("completeness", ""), checkm_contamination=ck.get("contamination", ""),
                    in_gtdb_r232=int(bool(base & gtdb))))
                n += 1
            tok = d.get("next_page_token")
            print(name, n, d.get("total_count"), flush=True)
            time.sleep(0.5)
            if not tok:
                break
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT, len(rows))


if __name__ == "__main__":
    main()
