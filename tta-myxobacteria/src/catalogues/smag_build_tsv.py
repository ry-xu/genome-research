"""SMAG (Ma et al. 2023 Nat Commun; Zenodo record 8223844 = doi 10.5281/zenodo.7341719, v1.0, 2023-08-08).
Supplementary Data 2 has GTDB r202 (GTDB-Tk) columns and a GTDB214 column. A block of rows has r202 values
duplicated from a neighbour (same closest-placement ANI/msa/red), so r214 is used as primary and r202 kept.
Polyangiaceae = f__Polyangiaceae in r214 OR r202 (flagged). GCA matched from NCBI PRJNA983538 by exact genome length
(+ SRR if bin name carries one). Usage: data_dir results_tsv want_list"""
import csv, sys, os, re, collections
ddir, out, want = sys.argv[1], sys.argv[2], sys.argv[3]
TARGET = {"JADLFU01", "JANWXX01", "DASTPP01", "JAXLYQ01", "CAYNEK01"}
R232 = set(l.strip() for l in open("data/outside_gtdb/r232_polyangiaceae_genera.txt") if l.strip())
def load(path):
    rows = list(csv.reader(open(path), delimiter="\t"))
    return rows[1], rows[2:]
h, recs = load(os.path.join(ddir, "SupplData2.Supplementary_Data_2.tsv"))
ix = {k: i for i, k in enumerate(h)}
sh, srecs = load(os.path.join(ddir, "SupplData1.Supplementary_Data_1.tsv"))
six = {k: i for i, k in enumerate(sh)}
samp = {}
for s in srecs:
    if len(s) < len(sh): continue
    samp[s[six["user_genome"]]] = s
ncbi = list(csv.DictReader(open(os.path.join(ddir, "ncbi_PRJNA983538_all.tsv")), delimiter="\t"))
by_len = collections.defaultdict(list)
for r in ncbi:
    if r["accession"].startswith("GCA_"):
        by_len[r["total_sequence_length"]].append(r)
print("total MAGs", len(recs))
def genus(t):
    m = re.search(r";g__([^;]*)", t or "")
    return m.group(1) if m else ""
n = collections.Counter(); wl = []
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genome_id", "taxonomy", "completeness", "contamination", "environment", "ncbi_ena_accession", "fasta_url",
                "taxonomy_r202", "genus_r214", "genus_r202", "poly_in", "r202_r214_phylum_agree", "target_genus", "download_set",
                "checkm2_completeness", "checkm2_contamination", "genome_size", "species_cluster", "ncbi_match_note", "gtdb_version"])
    for r in recs:
        t214 = r[ix["GTDB214"]]; t202 = r[ix["GTDB202_classification"]] if ";" in r[ix["GTDB202_classification"]] else ";".join(r[ix[c]] for c in ["GTDB202_classification","GTDB202_phylum","GTDB202_class","GTDB202_order","GTDB202_family","GTDB202_genus","GTDB202_species"])
        p214 = "f__Polyangiaceae" in t214; p202 = r[ix["GTDB202_family"]] == "f__Polyangiaceae"
        if not (p214 or p202): continue
        poly_in = "r214+r202" if (p214 and p202) else ("r214_only" if p214 else "r202_only")
        n[poly_in] += 1
        g214, g202 = genus(t214), r[ix["GTDB202_genus"]].replace("g__", "")
        ph214 = t214.split(";")[1] if ";" in t214 else ""
        agree = "yes" if ph214.replace("p__", "") in ("Myxococcota", "Myxococcota_A") and r[ix["GTDB202_phylum"]] == "p__Myxococcota" else "no"
        tg = (g214 in TARGET) or (g202 in TARGET)
        gid = r[ix["user_genome"]]
        prefix = re.sub(r"bin\..*$", "bin", gid)
        s = samp.get(prefix)
        env = f"{s[six['Ecosystem']]}; {s[six['Nation']]}; {s[six['Source']]}; sample {s[six['Sample']]}; {s[six['Bioproject']]}" if s else ""
        L = r[ix["length"]].replace(".0", "")
        cands = by_len.get(L, [])
        srr = re.match(r"([SED]RR\d+)bin", gid)
        if srr: cands = [c for c in cands if srr.group(1) in c["derived_from"]] or cands
        gca = cands[0]["accession"] if len(cands) == 1 else ""
        note = "" if len(cands) == 1 else ("no NCBI match by length" if not cands else f"{len(cands)} length matches")
        # download set: (a) target, (b) Polyangiaceae (r214 primary) with blank genus or genus not in R232 list
        gsel = g214 if p214 else g202
        if tg: dset = "a"
        elif (p214 and (gsel == "" or gsel not in R232)) or (poly_in == "r202_only" and agree == "yes" and (g202 == "" or g202 not in R232)):
            dset = "b"
        else: dset = ""
        if dset and gca: dset += "_ncbi"
        if dset: n["dl_" + dset] += 1; wl.append((gid, dset, gca))
        n["target"] += tg
        w.writerow([gid, t214 if p214 else t202, r[ix["completeness"]], r[ix["contamination"]], env, gca,
                    (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{gca}/download?include_annotation_type=GENOME_FASTA" if gca else
                     "only in Zenodo 10.5281/zenodo.7341719 split archive mag.tar.gz.* (~38 GB; no per-MAG URL)"),
                    t202, g214, g202, poly_in, agree, "yes" if tg else "no", dset,
                    r[ix["checkm2_completeness"]], r[ix["checkm2_contamination"]], L, r[ix["secondary_cluster"]], note,
                    "GTDB r214 column (primary); GTDB r202 GTDB-Tk columns"])
with open(want, "w") as fh:
    for x in wl: fh.write("\t".join(x) + "\n")
print(dict(n))
