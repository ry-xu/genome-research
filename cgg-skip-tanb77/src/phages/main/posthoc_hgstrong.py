# Post hoc: T1 phages with GC >= 32% and >= 5 CRISPR connections (strong host call, no skip signal).
import csv, subprocess, io
META = "data/uhgv/votus_metadata_extended.tsv.zst"
sets = {r["id"]: r["set"] for r in csv.DictReader(open("data/groups/run_ids.tsv"), delimiter="\t")}
stats = {r["id"]: r for r in csv.DictReader(open("results/main/target.tsv"), delimiter="\t") if sets.get(r["id"]) in ("T1", "Tf")}
txt = subprocess.run(["zstd", "-dcq", META], capture_output=True, text=True).stdout
for r in csv.DictReader(io.StringIO(txt), delimiter="\t"):
    i = r["uhgv_genome"]
    if i in stats and float(stats[i]["gc"]) >= 0.32 and int(r["crispr_total_connections"]) >= 5:
        s = stats[i]; fam = next((x for x in r["uhgv_taxonomy"].split(";") if x.startswith("vFAM-")), "-")
        print(sets[i], i, r["gc_percent"], "O", s["O"], "E", round(float(s["E"]), 1), "skip", s["skip_sites"], "tRNA_CCG", s["trna_ccg"],
              "conn", r["crispr_total_connections"], "agree", r["crispr_host_agreement"], fam, r["lifestyle"], r["checkv_quality"],
              "host", r["crispr_host_lineage_gtdb_r207"].split(";")[-1] or r["crispr_host_lineage_gtdb_r207"].split(";")[-2],
              "kmer", (r["kmers_host_lineage_gtdb_r207"].split(";") + ["", "", "", ""])[3], sep="\t")
