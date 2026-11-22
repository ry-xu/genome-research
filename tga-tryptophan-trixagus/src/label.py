"""Label flagged genomes: NCBI organism + lineage (Datasets API) and whether the flagged stop matches a known code.
Known genome-wide stop->sense codes in prokaryotes (NCBI code tables 4, 25, 34 + literature + this lab):
  TGA: Mollicutes except Acholeplasmatales/phytoplasmas (Trp; incl. Mycoplasmoidaceae, Entomoplasmatales, Spiroplasma),
       Gracilibacteria/BD1-5 and Absconditabacteria/SR1 (Gly), Hodgkinia, Zinderia, Nasuia, Stammera (Trp),
       Peptacetobacter-like? no. Any other TGA flag = check.
  TAG: code-34 archaea (Pyl): Methanosarcinales, Methanomassiliicoccales, Methermicoccaceae;
       stop-codon-census: Victivallales, Aristaeellaceae (Christensenellales), Lachnospiraceae near UBA1267.
  TAA: none known in any prokaryote.
usage: python src/label.py <acc> [...]   (or no args: every flagged genome in results/scores.tsv) -> results/labels.tsv
"""
import csv
import json
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
KNOWN = {"TGA": ["Mycoplasmatales", "Mycoplasmoidales", "Mycoplasmoidaceae", "Entomoplasmatales", "Spiroplasma",
                 "Mollicutes", "Gracilibacteria", "Absconditabacteria", "Hodgkinia", "Zinderia", "Nasuia", "Stammera"],
         "TAG": ["Methanosarcinales", "Methanomassiliicoccales", "Methermicoccaceae", "Methanosarcinia",
                 "Victivallales", "Christensenellales", "Lachnospiraceae"],
         "TAA": []}
NCBI_CODE_STOP = {"4": "TGA", "25": "TGA"}  # NCBI taxon code that reassigns a stop in prokaryotes (TGA=Trp, TGA=Gly)


def get(url, **kw):
    for i in range(4):
        try:
            r = requests.get(url, timeout=60, **kw)
            if r.status_code == 200:
                return r.json()
        except requests.RequestException:
            pass
        time.sleep(3 * (i + 1))
    return {}


def lineage(taxid):
    d = get(f"{API}/taxonomy/taxon/{taxid}/dataset_report")
    for rep in d.get("reports", []):
        t = rep.get("taxonomy", {})
        c = t.get("classification", {})
        gc = t.get("genetic_code", {}).get("primary", {}).get("id", "")
        return ";".join(c[k]["name"] for k in ("domain", "phylum", "class", "order", "family", "genus", "species")
                        if k in c), str(gc)
    return "", ""


def main():
    if len(sys.argv) > 1:
        todo = {(a, "") for a in sys.argv[1:]}
    else:
        with open(ROOT / "results" / "scores.tsv") as f:
            todo = {(r["acc"], r["stop"]) for r in csv.DictReader(f, delimiter="\t") if r["flag"] == "1"}
    out = []
    accs = sorted({a for a, _ in todo if not a.startswith("SYN")})
    info = {}
    batches = [accs[i:i + 20] for i in range(0, len(accs), 20)]
    while batches:
        b = batches.pop(0)
        d = get(f"{API}/genome/accession/{','.join(b)}/dataset_report")
        got = {x["accession"] for x in d.get("reports", [])}
        if len(b) > 1:
            batches += [[a] for a in b if a not in got]  # retry missing ones singly
        time.sleep(0.4)
        for x in d.get("reports", []):
            o = x.get("organism", {})
            bs = x.get("assembly_info", {}).get("biosample", {})
            attrs = {a.get("name"): a.get("value") for a in bs.get("attributes", [])}
            info[x["accession"]] = dict(org=o.get("organism_name", ""), taxid=o.get("tax_id", ""),
                                        project=x.get("assembly_info", {}).get("bioproject_accession", ""),
                                        source=attrs.get("isolation_source") or attrs.get("host") or
                                        attrs.get("env_medium") or "")
    for a, s in sorted(todo):
        i = info.get(a, {})
        lin, gc = lineage(i["taxid"]) if i.get("taxid") else ("", "")
        known = [k for k in KNOWN.get(s, []) if k in lin or k in i.get("org", "")]
        if not known and gc in NCBI_CODE_STOP and (s in ("", NCBI_CODE_STOP[gc])):
            known = [f"NCBI code {gc}"]
        out.append(dict(acc=a, stop=s, ncbi_gcode=gc, organism=i.get("org", ""), lineage=lin, project=i.get("project", ""),
                        source=i.get("source", ""), known=known[0] if known else ("SYN" if a.startswith("SYN") else "CHECK")))
    cols = ["acc", "stop", "known", "ncbi_gcode", "organism", "lineage", "project", "source"]
    with open(ROOT / "results" / "labels.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in out:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
            print("\t".join(str(r[c]) for c in cols))


if __name__ == "__main__":
    main()
