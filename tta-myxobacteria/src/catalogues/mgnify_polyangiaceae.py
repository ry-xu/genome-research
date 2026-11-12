"""MGnify genome catalogues: count genomes and f__Polyangiaceae per catalogue from genomes-all_metadata.tsv; list target genera.
usage: uv run --no-project python -I src/catalogues/mgnify_polyangiaceae.py
  reads data/catalogues/mgnify/<catalogue>.genomes-all_metadata.tsv (from <ftp-url>/genomes-all_metadata.tsv)
  -> results/catalogues/mgnify_polyangiaceae.tsv (all f__Polyangiaceae rows), mgnify_counts.tsv
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data/catalogues/mgnify"
R = ROOT / "results/catalogues"
T = ("g__JADLFU01", "g__JANWXX01", "g__DASTPP01", "g__JAXLYQ01", "g__CAYNEK01")
R232 = {l.strip() for l in open(R / "r232_polyangiaceae_genera.txt") if l.strip()}
csv.field_size_limit(10 ** 8)
rows, counts = [], []
for f in sorted(D.glob("*.genomes-all_metadata.tsv")):
    cat = f.name.split(".genomes")[0]
    n = npoly = nt = 0
    with open(f) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n += 1
            lin = r.get("Lineage", "")
            if "f__Polyangiaceae" not in lin:
                continue
            npoly += 1
            g = lin.split(";")[5] if lin.count(";") >= 5 else ""
            tgt = g in T
            nt += tgt
            fna = (r.get("FTP_download", "") or "").strip().replace("ftp://", "https://")  # GFF3 with the contigs in its ##FASTA section (no .fna on the FTP)
            rows.append(dict(catalogue=cat, genome_id=r["Genome"], taxonomy=lin, completeness=r.get("Completeness", ""),
                             contamination=r.get("Contamination", ""), genome_accession=r.get("Genome_accession", ""),
                             species_rep=r.get("Species_rep", ""), sample=r.get("Sample_accession", ""),
                             study=r.get("Study_accession", ""), country=r.get("Country", ""),
                             target_genus=int(tgt), genus_blank_or_not_in_r232=int(g[3:] not in R232),
                             fasta_url=fna))
    counts.append(dict(catalogue=cat, genomes=n, polyangiaceae=npoly, target_genus=nt))
    print(cat, n, npoly, nt)
with open(R / "mgnify_polyangiaceae.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
with open(R / "mgnify_counts.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(counts[0]), delimiter="\t")
    w.writeheader()
    w.writerows(counts)
