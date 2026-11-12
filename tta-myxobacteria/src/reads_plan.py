"""Raw-read sources for candidate MAGs (check 3 of PREREG.md "From an Arm B pass to a result").

Per candidate: the assembly's BioSample (GTDB metadata) -> its attributes from NCBI (derived_from / sample derived
from / run accessions named there) -> the source samples' runs from the ENA portal (platform, strategy, bases,
fastq URLs). A MAG BioSample with runs of its own is used as is.

usage: uv run --no-project python -I src/reads_plan.py  -> data/reads/plan.tsv
       uv run --no-project python -I src/reads_plan.py members  -> data/reads/plan_members.tsv
       (results/clade_members/candidates.tsv: GTDB members of the TTA clade outside the Arm B labels)
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "reads"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
IDRE = re.compile(r"\b((?:SAMN|SAMEA|SAMD)\d+|(?:SRS|ERS|DRS)\d+|(?:SRR|ERR|DRR)\d+)\b")
FIELDS = "run_accession,sample_accession,instrument_platform,instrument_model,library_strategy,library_layout," \
         "base_count,read_count,fastq_bytes,fastq_ftp"


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))


def biosample_attrs(acc):
    ids = get(EUTILS + "esearch.fcgi?db=biosample&term=" + urllib.parse.quote(acc + "[accn]"))
    uid = re.findall(r"<Id>(\d+)</Id>", ids)
    time.sleep(0.4)
    if not uid:
        return {}
    root = ET.fromstring(get(EUTILS + f"efetch.fcgi?db=biosample&id={uid[0]}"))
    time.sleep(0.4)
    at = {a.get("attribute_name"): (a.text or "") for a in root.iter("Attribute")}
    if not at and acc.startswith("SAMEA"):  # NCBI efetch returns no attributes for these; EBI BioSamples does
        ch = json.loads(get(f"https://www.ebi.ac.uk/biosamples/samples/{acc}.json")).get("characteristics", {})
        at = {k: " ".join(x.get("text", "") for x in v) for k, v in ch.items()}
    return at


def runs(acc):
    q = urllib.parse.urlencode(dict(accession=acc, result="read_run", fields=FIELDS, format="tsv"))
    txt = get(f"{ENA}?{q}")
    time.sleep(0.3)
    return list(csv.DictReader(txt.splitlines(), delimiter="\t")) if txt.strip() else []


def main(members=False):
    OUT.mkdir(parents=True, exist_ok=True)
    src_fn = ROOT / "results" / ("clade_members" if members else "gtdb") / "candidates.tsv"
    cand = list(csv.DictReader(open(src_fn), delimiter="\t"))
    rows = []
    for c in cand:
        bs = c["ncbi_biosample"]
        own = runs(bs)
        src, note = [], "own runs"
        if not own:
            at = biosample_attrs(bs)
            text = " ".join(v for k, v in at.items() if "deriv" in k.lower() or "run" in k.lower() or "sra" in k.lower())
            src = sorted(set(IDRE.findall(text)))
            note = f"derived from {','.join(src) or 'nothing found'}"
            seen = set()
            for s in src:
                for r in runs(s):
                    if r["run_accession"] not in seen:
                        seen.add(r["run_accession"])
                        own.append(r)
        for r in own or [{}]:
            rows.append(dict(gtdb_acc=c["gtdb_acc"], codon=c["codon"], genus=c["genus"], biosample=bs, note=note,
                             **{k: r.get(k, "") for k in FIELDS.split(",")}))
        n = len(own)
        gb = sum(float(r.get("base_count") or 0) for r in own) / 1e9
        print(f"{c['gtdb_acc']} {c['codon']} {c['genus']} {bs}: {note}; {n} runs, {gb:.1f} Gbp, "
              f"{sorted({r.get('instrument_platform', '') for r in own})}", flush=True)
    with open(OUT / ("plan_members.tsv" if members else "plan.tsv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    import sys
    main(sys.argv[1:] == ["members"])
