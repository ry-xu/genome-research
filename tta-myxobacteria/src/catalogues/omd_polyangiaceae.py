"""OMD / OMDB v2.0 (Ocean Microbiomics Database, omdb.microbiomics.io): f__Polyangiaceae genomes.

Taxonomy comes from the website's own table API (POST https://omdb.microbiomics.io/api/genome-cols,
the endpoint the genome-collection page uses; read-only search), gtdb_version = R226.
FASTA URLs come from https://sunagawalab.ethz.ch/share/microbiomics/ocean/db/2.0/data/catalogs/OMDBv2.0_data.tsv.gz
(Last-Modified 2025-05-23; 274,282 genomes).
Outputs: results/catalogues/omd_polyangiaceae.tsv, data/catalogues/omd/download_list.txt
"""
import gzip, json, urllib.request
from pathlib import Path

ROOT = Path(".")
D = ROOT / "data/catalogues/omd"
TARGET = {"JADLFU01", "JANWXX01", "DASTPP01", "JAXLYQ01", "CAYNEK01"}
R232 = set((ROOT / "results/catalogues/r232_polyangiaceae_genera.txt").read_text().split())


def query(term, length=5000):
    body = json.dumps({"repo_id": "ocean", "global_search_value": term, "search_value_pair": {},
                       "start": 0, "length": length, "window": 0, "fetchall": 0,
                       "order": [], "order_cols": []}).encode()
    req = urllib.request.Request("https://omdb.microbiomics.io/api/genome-cols", data=body,
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)


res = query("Polyangiaceae")
(D / "api_search_Polyangiaceae.json").write_text(json.dumps(res))
print("total genomes in OMDB:", res["total_count"], "; search hits:", res["filtered_count"])
rows = [r for r in res["data"] if r["family"] == "Polyangiaceae"]
for g in TARGET:
    print(g, "hits:", query(g, 10)["filtered_count"])

urls = {}
with gzip.open(D / "OMDBv2.0_data.tsv.gz", "rt") as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        urls[f[0]] = f[3]

cols = ["genome_id", "taxonomy", "completeness", "contamination", "environment",
        "ncbi_ena_accession", "fasta_url", "gtdb_version", "study", "sample", "publication",
        "target_genus", "download_reason"]
dl = []
with open(ROOT / "results/catalogues/omd_polyangiaceae.tsv", "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: r["genome"]):
        g = r["genus"]
        g_clean = "" if g.startswith("Unknown") else g
        tax = (f"d__{r['domain']};p__{r['phylum']};c__{r['class_']};o__{r['order']};"
               f"f__{r['family']};g__{g_clean};s__{'' if r['species'].startswith('Unknown') else r['species']}")
        reason = ""
        if g_clean in TARGET:
            reason = "target_genus"
        elif g_clean == "" or g_clean not in R232:
            reason = "blank_or_non_r232_genus"
        acc = ";".join(x for x in [r.get("biosample") or "", r.get("bioproject") or ""] if x)
        url = urls.get(r["genome"], "")
        if reason:
            dl.append((r["genome"], url))
        fo.write("\t".join(str(x) for x in [
            r["genome"], tax, r["completeness"], r["contamination"],
            f"{r.get('environment') or ''};{r.get('public_sample_link') or ''}", acc, url,
            r["gtdb_version"], r["sushistudy"], r["shushisample"], r.get("publication") or "",
            g_clean if g_clean in TARGET else "", reason]) + "\n")
with open(D / "download_list.txt", "w") as fo:
    for gid, url in dl:
        fo.write(f"{gid}\t{url}\n")
print("f__Polyangiaceae:", len(rows), "; download set:", len(dl))
