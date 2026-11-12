"""SPIRE v1: list f__Polyangiaceae MAGs and pick download set.

Inputs (all under data/catalogues/spire/):
  spire_v1_genome_metadata.tsv.gz  bulk table (README: SPIRE v01, 2023-09)
  api_family_Polyangiaceae.tsv     live website API (newer GTDB taxonomy)
  microntology_annotations.tsv.gz  sample -> study name / ENA project
Outputs:
  results/catalogues/spire_polyangiaceae.tsv
  data/catalogues/spire/download_list.txt (genome_id <tab> url)
"""
import csv, gzip, sys
from pathlib import Path

ROOT = Path(".")
D = ROOT / "data/catalogues/spire"
OUT = ROOT / "results/catalogues/spire_polyangiaceae.tsv"
TARGET = {"JADLFU01", "JANWXX01", "DASTPP01", "JAXLYQ01", "CAYNEK01"}
R232 = set((ROOT / "results/catalogues/r232_polyangiaceae_genera.txt").read_text().split())

csv.field_size_limit(sys.maxsize)


def genus_of(tax):
    for p in tax.split(";"):
        if p.startswith("g__"):
            return p[3:]
    return ""


bulk = {}
total = 0
with gzip.open(D / "spire_v1_genome_metadata.tsv.gz", "rt") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        total += 1
        if "f__Polyangiaceae;" in row["classification"]:
            bulk[row["genome_id"]] = row

api = {}
with open(D / "api_family_Polyangiaceae.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        api[row["spire_id"]] = row

env = {}
with gzip.open(D / "microntology_annotations.tsv.gz", "rt") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        env[row["ena_ers_sample_id"]] = (row["study_name"], row["ena_project_id"])

ids = sorted(set(bulk) | set(api))
cols = ["genome_id", "taxonomy", "completeness", "contamination", "environment",
        "ncbi_ena_accession", "fasta_url", "taxonomy_website_api", "sample_id",
        "in_bulk_table", "in_api_table", "target_genus", "download_reason"]
dl = []
n_target = 0
with open(OUT, "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for gid in ids:
        b = bulk.get(gid, {})
        a = api.get(gid, {})
        tax_b = b.get("classification", "")
        tax_a = a.get("gtdb_classification", "")
        sample = b.get("derived_from_sample") or a.get("sample_id", "")
        comp = b.get("completeness") or a.get("completeness", "")
        cont = b.get("contamination") or a.get("contamination", "")
        study, proj = env.get(sample, ("", ""))
        envs = f"study={study};ena_project={proj};sample={sample}" if study else f"sample={sample}"
        url = a.get("download_mag") or f"https://spire.embl.de/download_file/{gid}"
        genera = {genus_of(t) for t in (tax_b, tax_a) if t}
        tg = sorted(genera & TARGET)
        reason = ""
        if tg:
            reason = "target_genus"
            n_target += 1
        elif any(g == "" or g not in R232 for g in genera):
            reason = "blank_or_non_r232_genus"
        if reason:
            dl.append((gid, url))
        w.writerow([gid, tax_b, comp, cont, envs, "", url, tax_a, sample,
                    int(bool(b)), int(bool(a)), ",".join(tg), reason])

with open(D / "download_list.txt", "w") as fo:
    for gid, url in dl:
        fo.write(f"{gid}\t{url}\n")

print(f"bulk total genomes: {total}")
print(f"Polyangiaceae: bulk {len(bulk)}, api {len(api)}, union {len(ids)}")
print(f"target-genus genomes: {n_target}; download set: {len(dl)}")
