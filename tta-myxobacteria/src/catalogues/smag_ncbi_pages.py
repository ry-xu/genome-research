"""Page NCBI Datasets genome reports for SMAG BioProject PRJNA983538 -> TSV.
Usage: outdir  (writes ncbi_PRJNA983538_all.tsv)"""
import json, sys, time, urllib.request, urllib.parse, csv, os
out = os.path.join(sys.argv[1], "ncbi_PRJNA983538_all.tsv")
base = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/bioproject/PRJNA983538/dataset_report"
tok = None
rows = []
while True:
    q = {"page_size": "1000", "filters.assembly_source": "all"}
    if tok:
        q["page_token"] = tok
    url = base + "?" + urllib.parse.urlencode(q)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                d = json.load(r)
            break
        except Exception as e:
            print("retry", attempt, e, file=sys.stderr); time.sleep(5)
    else:
        sys.exit("failed")
    for rep in d.get("reports", []):
        bs = rep.get("assembly_info", {}).get("biosample", {})
        attrs = {a.get("name"): a.get("value") for a in bs.get("attributes", [])}
        rows.append([rep["accession"], attrs.get("isolate", ""), attrs.get("derived_from", ""),
                     rep.get("assembly_stats", {}).get("total_sequence_length", ""),
                     rep.get("organism", {}).get("organism_name", ""),
                     attrs.get("geo_loc_name", ""), attrs.get("env_local_scale", ""), attrs.get("note", ""),
                     bs.get("accession", "")])
    tok = d.get("next_page_token")
    print(len(rows), file=sys.stderr)
    if not tok:
        break
    time.sleep(0.5)
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "isolate", "derived_from", "total_sequence_length", "organism", "geo_loc_name", "env_local_scale", "note", "biosample"])
    w.writerows(rows)
print("wrote", out, len(rows))
