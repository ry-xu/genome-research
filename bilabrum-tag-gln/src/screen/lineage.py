"""Build data/lineage.tsv (taxid, name, lineage) for every unit's taxid.

Starts from frameshift-codes-euk's data/meta/taxlineage.tsv; fetches the rest from NCBI taxonomy (efetch XML, 200 ids
per request) into data/dl_tax/ (downloaded, untrusted: parsed as XML only). Cache files are named by a hash of their
id list (the first version named them by batch index, so a later run with new taxids reread the old files). Taxids
still missing are filled from the local NCBI taxdump in data/dl_taxdump (nodes.dmp, names.dmp).
usage: lineage.py UNITS.tsv > data/lineage.tsv
"""
import csv
import hashlib
import os
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET

FC = "lab/work/frameshift-codes-euk/data/meta/taxlineage.tsv"
DL = "lab/work/stop-sense-euk/data/dl_tax"
EF = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&retmode=xml&id="

TD = "lab/work/stop-sense-euk/data/dl_taxdump"


def from_taxdump(tids):
    parent, names, merged = {}, {}, {}
    for line in open(os.path.join(TD, "nodes.dmp")):
        x = line.split("\t|\t")
        parent[x[0]] = x[1]
    for line in open(os.path.join(TD, "names.dmp")):
        x = line.split("\t|\t")
        if x[3].startswith("scientific name"):
            names[x[0]] = x[1]
    if os.path.exists(os.path.join(TD, "merged.dmp")):
        for line in open(os.path.join(TD, "merged.dmp")):
            x = line.split("\t|\t")
            merged[x[0]] = x[1].split("\t")[0]
    out = {}
    for t in tids:
        t2 = merged.get(t, t)
        if t2 not in parent:
            continue
        chain, n = [], t2
        while n != "1" and n in parent:
            chain.append(names.get(n, n))
            n = parent[n]
        out[t] = (names.get(t2, ""), "; ".join(reversed(chain)))
    return out


def main(units_tsv):
    want = {r["taxid"] for r in csv.DictReader(open(units_tsv), delimiter="\t") if r["taxid"]}
    lin = {}
    for r in csv.DictReader(open(FC), delimiter="\t"):
        lin[r["taxid"]] = ("", r["lineage"])
    os.makedirs(DL, exist_ok=True)
    todo = sorted(want - set(lin))
    for i in range(0, len(todo), 200):
        ids = ",".join(todo[i:i + 200])
        f = os.path.join(DL, f"ids_{hashlib.sha1(ids.encode()).hexdigest()[:12]}.xml")
        if not os.path.exists(f):
            with urllib.request.urlopen(EF + ids, timeout=120) as r, open(f, "wb") as o:
                o.write(r.read())
            time.sleep(1)
    for f in sorted(os.listdir(DL)):
        for t in ET.parse(os.path.join(DL, f)).getroot().findall("Taxon"):
            tid, name = t.findtext("TaxId"), t.findtext("ScientificName")
            lin[tid] = (name, f"cellular organisms; {t.findtext('Lineage')}; {name}".replace("cellular organisms; cellular organisms", "cellular organisms"))
            for a in t.findall("AkaTaxIds/TaxId"):
                lin.setdefault(a.text, lin[tid])
    missing = want - set(lin)
    if missing:
        lin.update(from_taxdump(missing))
    print("taxid\tname\tlineage")
    for tid in sorted(want, key=int):
        n, l = lin.get(tid, ("", ""))
        print(f"{tid}\t{n}\t{l}")
    print(f"lineages: {sum(1 for t in want if t in lin)} of {len(want)}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1])
