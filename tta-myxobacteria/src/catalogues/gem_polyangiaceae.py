"""GEM (Nayfach et al. 2021): list f__Polyangiaceae MAGs and pick download set.

Input: data/catalogues/gem/genome_metadata.tsv (https://portal.nersc.gov/GEM/genomes/genome_metadata.tsv,
Last-Modified 2020-11-30; GTDB-Tk taxonomy, GTDB r89).
NOTE: the header has 20 names but column 15 (labelled 'ecosystem') holds the GTDB taxonomy;
columns 16-18 are ecosystem_category/ecosystem_type/habitat, 19-20 longitude/latitude.
Outputs:
  results/catalogues/gem_polyangiaceae.tsv            all f__Polyangiaceae
  results/catalogues/gem_polyangiales_other_families.tsv  other o__Polyangiales MAGs (r89 families
      such as f__NIC37A-2 that are genera inside Polyangiaceae in R232, or blank family)
  data/catalogues/gem/download_list.txt
"""
from pathlib import Path

ROOT = Path(".")
D = ROOT / "data/catalogues/gem"
TARGET = {"JADLFU01", "JANWXX01", "DASTPP01", "JAXLYQ01", "CAYNEK01"}
R232 = set((ROOT / "results/catalogues/r232_polyangiaceae_genera.txt").read_text().split())
URL = "https://portal.nersc.gov/GEM/genomes/fna/{}.fna.gz"
cols = ["genome_id", "taxonomy", "completeness", "contamination", "environment",
        "ncbi_ena_accession", "fasta_url", "img_metagenome_id", "mimag_quality", "otu_id",
        "target_genus", "download_reason"]


def rank(tax, p):
    for x in tax.split(";"):
        if x.startswith(p):
            return x[len(p):]
    return ""


rows_f, rows_o, dl = [], [], []
total = 0
lines = (D / "genome_metadata.tsv").read_text().splitlines()
for line in lines[1:]:
    f = line.split("\t")
    total += 1
    gid, mg, comp, cont, mimag, otu, tax = f[0], f[1], f[9], f[10], f[12], f[13], f[14]
    env = f"{f[15]};{f[16]};{f[17]};lon={f[18]};lat={f[19]}"
    if "o__Polyangiales" not in tax:
        continue
    g = rank(tax, "g__")
    url = URL.format(gid)
    reason = ""
    if g in TARGET:
        reason = "target_genus"
    elif g == "" or g not in R232:
        reason = "blank_or_non_r232_genus"
    rec = [gid, tax, comp, cont, env, "", url, mg, mimag, otu, g if g in TARGET else "", reason]
    if "f__Polyangiaceae" in tax:
        rows_f.append(rec)
        if reason:
            dl.append((gid, url))
    else:
        rec[-1] = "other_polyangiales_family"
        rows_o.append(rec)
        dl.append((gid, url))

for path, rows in [(ROOT / "results/catalogues/gem_polyangiaceae.tsv", rows_f),
                   (ROOT / "results/catalogues/gem_polyangiales_other_families.tsv", rows_o)]:
    with open(path, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            fo.write("\t".join(r) + "\n")
with open(D / "download_list.txt", "w") as fo:
    for gid, url in dl:
        fo.write(f"{gid}\t{url}\n")
print(f"total genomes {total}; f__Polyangiaceae {len(rows_f)}; other Polyangiales {len(rows_o)}")
print(f"target genus: {sum(1 for r in rows_f if r[10])}; download list {len(dl)}")
