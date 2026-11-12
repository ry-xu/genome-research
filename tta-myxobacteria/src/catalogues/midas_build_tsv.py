"""MiDAS / Singleton et al. 2021 HQ MAGs: Supplementary Data 3 (GTDB r89 via GTDB-Tk v1.0.2)
-> results/catalogues/midas_polyangiaceae.tsv, with GCA accessions from NCBI BioProject PRJNA629478.
Usage: data_dir results_tsv"""
import csv, json, sys, os, time, urllib.request, urllib.parse
ddir, out = sys.argv[1], sys.argv[2]
supp = os.path.join(ddir, "singleton2021_SupplData3.tsv")
ncbi_tsv = os.path.join(ddir, "ncbi_PRJNA629478_all.tsv")
TARGET = {"JADLFU01", "JANWXX01", "DASTPP01", "JAXLYQ01", "CAYNEK01"}
if not os.path.exists(ncbi_tsv):
    base = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/bioproject/PRJNA629478/dataset_report"
    tok, rows = None, []
    while True:
        q = {"page_size": "1000", "filters.assembly_source": "all"}
        if tok: q["page_token"] = tok
        with urllib.request.urlopen(base + "?" + urllib.parse.urlencode(q), timeout=120) as r:
            d = json.load(r)
        for rep in d.get("reports", []):
            bs = rep["assembly_info"].get("biosample", {})
            sn = {s["label"]: s["value"] for s in bs.get("sample_ids", []) if "label" in s}
            attrs = {a.get("name"): a.get("value") for a in bs.get("attributes", [])}
            rows.append([rep["accession"], sn.get("Sample name", ""), attrs.get("isolate", ""),
                         rep.get("wgs_info", {}).get("wgs_project_accession", ""),
                         rep.get("assembly_stats", {}).get("total_sequence_length", ""),
                         rep.get("organism", {}).get("organism_name", "")])
        tok = d.get("next_page_token")
        if not tok: break
        time.sleep(0.5)
    with open(ncbi_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "sample_name", "isolate", "wgs_project", "total_sequence_length", "organism"])
        w.writerows(rows)
ncbi = list(csv.DictReader(open(ncbi_tsv), delimiter="\t"))
by_wgs = {r["wgs_project"][:6]: r for r in ncbi if r["wgs_project"]}
by_name = {}
for r in ncbi:
    for k in (r["sample_name"], r["isolate"]):
        if k: by_name[k.replace("_fly", "").replace(".fa", "")] = r
recs = list(csv.DictReader(open(supp), delimiter="\t"))
print("total HQ MAGs", len(recs))
n_poly = n_tgt = 0
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genome_id", "taxonomy", "completeness", "contamination", "environment", "ncbi_ena_accession",
                "fasta_url", "genbank_wgs", "genus_r89", "target_genus", "genome_size", "midas3_7_tax", "gtdb_version"])
    for r in recs:
        tax = r["GTDBTax"]
        if "f__Polyangiaceae" not in tax: continue
        n_poly += 1
        g = tax.split(";g__")[1].split(";")[0]
        tg = g in TARGET; n_tgt += tg
        wgs = r["NCBI_accession_number"].strip()
        mag = r["MAG"].replace(".fa", "")
        hit = by_wgs.get(wgs[:6]) or by_name.get(mag.replace("_cln", "")) or by_name.get(mag)
        gca = hit["accession"] if hit else ""
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{gca}/download?include_annotation_type=GENOME_FASTA" if gca else ""
        wwtp = r["MAG"].split("_")[0]
        w.writerow([r["MAG"], tax, r["Comp"], r["Cont"],
                    f"activated sludge, Danish full-scale WWTP ({wwtp}); Singleton et al. 2021",
                    gca, url, wgs, g, "yes" if tg else "no", r["TotBP"], r["MiDAS3_7Tax"],
                    "GTDB r89 (GTDB-Tk v1.0.2)"])
print("polyangiaceae", n_poly, "target", n_tgt)
