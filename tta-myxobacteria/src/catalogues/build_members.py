"""Compile results/catalogues/members.tsv from the screen tables.

Call rule (fixed from controls before the final compile; see searched.md):
  candidate member  = n_markers >= 50 and top_hit_clade_frac >= 0.5 and clade_aai > out_aai
  borderline (listed separately in screen_borderline.tsv) = 0.2 <= top_hit_clade_frac < 0.5
Redundant copy = same BioSample as one of the 10 references, or pyskani ANI >= 99.5% to one (redundancy_*.tsv).
usage: uv run --no-project python -I src/catalogues/build_members.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results/catalogues"
D = ROOT / "data/catalogues"
REF_BS = {"SAMN30304787": "GCA_025057345.1", "SAMEA117372376": "GCA_964605685.1", "SAMEA114248066": "GCA_963169165.1",
          "SAMN37623916": "GCA_046351315.1", "SAMEA117295240": "GCA_964396175.1", "SAMN37960669": "GCA_035762175.1",
          "SAMN37334137": "GCA_033578855.1", "SAMN16439656": "GCA_020849655.1", "SAMN27108603": "GCA_025453595.1",
          "SAMEA114248065": "GCA_963169955.1"}


# BioSample attributes looked up by hand for NCBI hits whose assembly report has no isolation source / CheckM
# (NCBI efetch biosample; EBI BioSamples JSON for SAMEA*), 2026-11-11.
ENV_EXTRA = {
    "GCA_982529355.1": "NOT a pure isolate: predator grown in co-culture with its prey and assembled from co-culture reads "
                       "after removing prey reads, deposited at ENA as a MAG (Demina et al., F1000Research 14:1153 v2, "
                       "2026-04-28); 'Predatory bacterium 1-FT3.2 obtained from subarctic soil using Mucilaginibacter "
                       "cryoferens FT3.2 as prey', Finland, collected 2024, Univ. Helsinki (SAMEA122172943, PRJNA1337162)",
    "spire_mag_00879468": "wastewater: SPIRE study Ciesielski_2019_wastewater (SPIRE sample page https://spire.embl.de/spire/v1/sample/mgs763609)",
    "GCA_035762045.1": "= SMAG MAG S400bin.4 (SMAG genus blank): agricultural land soil, China, sample China_S400 "
                       "(SAMN37960678, PRJNA983538); released 2024-01-17, before the R232 cutoff but absent from R232 "
                       "(probably GTDB quality filter: 3,387 contigs, ~50-56% complete)",
    "GCA_057137765.1": "MAG05673 (SRR34280936_bin.33): paddy field soil, wetland ecosystem, China 27.20N 114.73E, collected "
                       "2022-08-31, Zhejiang University (SAMN56415168, PRJNA1430164)",
}
N_MARKERS = 437
# completeness/contamination from a catalogue that holds the same MAG (SMAG Zenodo 8223844 table; smag_polyangiaceae.tsv)
CC_EXTRA = {"GCA_035762045.1": ("56.19", "3.57", "SMAG CheckM (S400bin.4); SMAG CheckM2 49.15/4.1"),
            "GCA_982529355.1": ("81.9", "0.38", "CheckM2, from the ENA assembly description (audit 2026-11-11 12:24; the marker proxy was ~93)")}
# redundancy among catalogue genomes themselves (pyskani, scratch/ani4.py, 2026-11-11 12:3x); the table above is vs the 10 refs
RED_EXTRA = {"GCA_060838125.1": "ANI 99.95 to GCA_060838225.1 (same strain, same project; keep GCA_060838225.1)"}


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t")) if Path(p).exists() else []


def call(r):
    if r.get("n_markers", "NA") == "NA":
        return "NA"
    f, n = float(r["top_hit_clade_frac"]), int(r["n_markers"])
    if n >= 50 and f >= 0.5 and float(r["clade_aai"]) > float(r["out_aai"]):
        return "member"
    if f >= 0.2:
        return "borderline"
    return "no"


red = {}
for p in R.glob("redundancy_*.tsv"):
    for r in tsv(p):
        red[r["genome"]] = r

meta = {}
for r in tsv(R / "ncbi_myxococcota_assemblies.tsv"):
    meta[("ncbi", r["accession"])] = dict(
        source="NCBI GenBank/RefSeq (Datasets v2)", taxonomy=f"NCBI: {r['organism']} (taxid {r['taxid']})",
        url=f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{r['accession']}/download?include_annotation_type=GENOME_FASTA",
        compl=r["checkm_completeness"], cont=r["checkm_contamination"], cc_src="NCBI CheckM",
        env="; ".join(x for x in (r["isolation_source"], r["geo"], r["metagenome_source"], f"released {r['release_date']}") if x),
        biosample=r["biosample"], local=str(D / "ncbi/fna" / f"{r['accession']}.fna.gz"))
SRC = {"spire": ("SPIRE v01 (bulk table + website taxonomy)", D / "spire/mags", ".fa.gz"),
       "gem": ("GEM (Nayfach 2021)", D / "gem/fna", ".fna.gz"),
       "anmagc": ("AnMAGC (anammox MAG catalogue, Wang et al. 2024)", D / "anmagc", ".fa.gz"),
       "brown23": ("Brown et al. 2023 activated-sludge MAGs (Zenodo 10028566)", D / "brown23", ".fna.gz"),
       "mgnify": ("MGnify genome catalogues", D / "mgnify/fna", ".fna.gz")}
for s, (name, ddir, ext) in SRC.items():
    for r in tsv(R / f"{s}_polyangiaceae.tsv"):
        gid = r["genome_id"]
        tax = r["taxonomy"]
        if r.get("taxonomy_website_api"):
            tax += f" | website: {r['taxonomy_website_api']}"
        if r.get("genus_r214"):
            tax += f" | r214 genus: {r['genus_r214']}"
        if s == "mgnify":
            name2 = f"MGnify {r['catalogue']}"
        else:
            name2 = name
        bs = ""
        for k in ("sample", "sample_id", "environment"):
            v = r.get(k, "") or ""
            for tok in v.replace(";", " ").replace("=", " ").split():
                if tok.startswith(("SAMN", "SAMEA", "SAMD")):
                    bs = tok
        meta[(s, gid)] = dict(source=name2, taxonomy=tax, url=r.get("fasta_url", ""), compl=r.get("completeness", ""),
                              cont=r.get("contamination", ""), cc_src="catalogue",
                              env=r.get("environment", "") or f"{r.get('country', '')} {r.get('study', '')} {r.get('sample', '')}".strip(),
                              biosample=bs, local=str(ddir / f"{gid}{ext}"))

out, border = [], []
for s in ["ncbi", "spire", "gem", "anmagc", "brown23", "mgnify"]:
    for r in tsv(R / f"screen_{s}.tsv"):
        c = call(r)
        if c not in ("member", "borderline"):
            continue
        m = meta.get((s, r["genome"]), {})
        rd = red.get(r["genome"], {})
        bsr = REF_BS.get(m.get("biosample", ""), "")
        redundant = []
        if r["genome"] in RED_EXTRA:
            redundant.append(RED_EXTRA[r["genome"]])
        if bsr:
            redundant.append(f"same BioSample as {bsr}")
        if rd and rd.get("ani") not in ("<80", "", None) and float(rd["ani"]) >= 99.5:
            redundant.append(f"ANI {rd['ani']} to {rd['best_ref']}")
        row = dict(
            id=r["genome"], source=m.get("source", s), download_url=m.get("url", ""),
            local_path=m.get("local", "") if Path(m.get("local", "")).exists() else "",
            taxonomy=m.get("taxonomy", ""),
            screen=(f"marker AAI {r['clade_aai']}% to {r['clade_ref']} ({r['clade_ref_genus']}, {r['clade_n']} markers); "
                    f"top-hit clade fraction {r['top_hit_clade_frac']} of {r['n_markers']} markers; "
                    f"nearest outgroup {r['out_ref']} {r['out_aai']}%"
                    + (f"; skani ANI {rd['ani']}% (AF {rd['af_query']}%) to {rd['best_ref']}" if rd and rd.get('best_ref') != 'none' else "; skani ANI <80% to all refs" if rd else "")),
            completeness=m.get("compl", "") or f"~{100 * int(r['n_markers']) / N_MARKERS:.0f} (proxy: {r['n_markers']}/{N_MARKERS} clade core markers found)",
            contamination=m.get("cont", "") or "not computed",
            cc_source=m.get("cc_src", "") if m.get("compl") else "marker proxy (no CheckM in NCBI report)",
            environment=ENV_EXTRA.get(r["genome"], "") + ("; " if r["genome"] in ENV_EXTRA else "") + m.get("env", ""), biosample=m.get("biosample", ""),
            redundant_copy="; ".join(redundant) if redundant else "no", call=c)
        if r["genome"] in CC_EXTRA:
            row["completeness"], row["contamination"], row["cc_source"] = CC_EXTRA[r["genome"]]
        (out if c == "member" else border).append(row)
for p, rows in ((R / "members.tsv", out), (R / "screen_borderline.tsv", border)):
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list((out or border or [{"id": ""}])[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(p.name, len(rows))
