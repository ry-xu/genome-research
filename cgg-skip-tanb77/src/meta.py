"""GTDB R232 representative table + which genomes are already in the genome cache.
-> data/reps.tsv: acc gca domain phylum class order family genus species compl contam contigs n50 gc size cdens ttable local_path"""
import gzip, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import GENOMES, GTDB_META as META  # noqa: E402
# the genome cache (sharded: genomes/<shard>/<acc>.fna.gz)
STORES = [GENOMES]
local = {}
for d in STORES:
    files = list(d.glob("*/*.fna.gz")) if d == GENOMES else list(d.glob("*.fna*"))
    for f in files:
        acc = f.name.split(".fna")[0]
        local.setdefault(acc.rsplit(".", 1)[0], str(f))
out = open(Path(__file__).resolve().parents[1] / "data/reps.tsv", "w")
cols = ["acc", "gca", "domain", "phylum", "class", "order", "family", "genus", "species", "compl", "contam", "contigs",
        "n50", "gc", "size", "cdens", "ttable", "level", "local"]
out.write("\t".join(cols) + "\n")
n = nl = 0
for fn in ["bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"]:
    with gzip.open(META / fn, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if r[ix["gtdb_representative"]] != "t":
                continue
            tax = [t[3:] for t in r[ix["gtdb_taxonomy"]].split(";")]
            acc = r[ix["accession"]]
            gca = acc[3:]
            gb = r[ix["ncbi_genbank_assembly_accession"]]
            p = local.get(gca.rsplit(".", 1)[0]) or local.get(gb.rsplit(".", 1)[0], "")
            n += 1; nl += bool(p)
            out.write("\t".join([acc, gca] + tax + [r[ix[c]] for c in ["checkm2_completeness", "checkm2_contamination",
                      "contig_count", "n50_contigs", "gc_percentage", "genome_size", "coding_density",
                      "ncbi_translation_table", "ncbi_assembly_level"]] + [p]) + "\n")
print(n, "reps;", nl, "on disk")
