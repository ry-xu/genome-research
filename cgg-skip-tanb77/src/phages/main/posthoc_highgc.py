# Post hoc (descriptive): T1 phages GC >= 32% (no skip signal) vs GC < 32%: host-call evidence and lifestyle.
import csv, subprocess, io, collections, statistics as st
META = "data/uhgv/votus_metadata_extended.tsv.zst"
sets = {r["id"]: r["set"] for r in csv.DictReader(open("data/groups/run_ids.tsv"), delimiter="\t")}
stats = {r["id"]: r for r in csv.DictReader(open("results/main/target.tsv"), delimiter="\t") if sets.get(r["id"]) == "T1"}
txt = subprocess.run(["zstd", "-dcq", META], capture_output=True, text=True).stdout
m = {r["uhgv_genome"]: r for r in csv.DictReader(io.StringIO(txt), delimiter="\t") if r["uhgv_genome"] in stats}
def order(l): return next((x for x in l.split(";") if x.startswith("o__")), "o__?")
for name, ids in (("GC>=32", [i for i in stats if float(stats[i]["gc"]) >= 0.32]),
                  ("GC<32", [i for i in stats if float(stats[i]["gc"]) < 0.32])):
    rows = [m[i] for i in ids]
    O = sum(int(stats[i]["O"]) for i in ids); S = sum(int(stats[i]["skip_sites"]) for i in ids)
    E = sum(float(stats[i]["E"]) for i in ids)
    print(f"== {name}: n={len(ids)} O={O} E={E:.0f} O/E={O/E:.2f} skip={S} skip/O={S/O:.3f}")
    print("  crispr connections median:", st.median(int(r["crispr_total_connections"]) for r in rows),
          " agreement median:", st.median(float(r["crispr_host_agreement"]) for r in rows))
    print("  kmer host order:", collections.Counter(order(r["kmers_host_lineage_gtdb_r207"]) if r["kmers_host_lineage_gtdb_r207"] not in ("", "NULL") else "none" for r in rows).most_common(5))
    print("  lifestyle:", collections.Counter(r["lifestyle"] for r in rows).most_common(), " checkv:", collections.Counter(r["checkv_quality"] for r in rows).most_common())
    print("  host markers>0:", sum(1 for r in rows if r["checkv_host_markers"] not in ("", "0", "NULL")), " trimmed:", collections.Counter(r["checkv_trimmed"] for r in rows).most_common())
    print("  crispr_host_range_red median:", st.median(float(r["crispr_host_range_red"]) for r in rows if r["crispr_host_range_red"] not in ("", "NULL")))
