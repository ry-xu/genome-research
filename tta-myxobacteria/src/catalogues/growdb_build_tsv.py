"""GROWdb (Borton et al. 2024 Nature 637:103; Supplementary Table MOESM4 Tab1 = 2,093 dereplicated MAGs, GTDB-Tk)
-> results/catalogues/growdb_polyangiaceae.tsv. MAGs are in NCBI BioProject PRJNA946291. Usage: in_tsv out_tsv"""
import csv, sys
TARGET = {"JADLFU01", "JANWXX01", "DASTPP01", "JAXLYQ01", "CAYNEK01"}
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
print("total", len(rows))
n = t = 0
with open(sys.argv[2], "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genome_id", "taxonomy", "completeness", "contamination", "environment", "ncbi_ena_accession", "fasta_url", "genus", "target_genus", "quality", "gtdb_version"])
    for r in rows:
        tax = r["GTDB classification"]
        if "f__Polyangiaceae" not in tax: continue
        n += 1
        g = r["Genus"].replace("g__", ""); t += g in TARGET
        w.writerow([r["user_genome"], tax, r["checkM completion"], r["checkM contamination"],
                    "river surface water/sediment, North America (WHONDRS/GROWdb)", "NCBI BioProject PRJNA946291 (MAG GCA not mapped)",
                    "NCBI PRJNA946291; Zenodo 10.5281/zenodo.8173286", g, "yes" if g in TARGET else "no", r["Quality"], "GTDB-Tk v2.1.1, GTDB r207"])
print("polyangiaceae", n, "target", t)
