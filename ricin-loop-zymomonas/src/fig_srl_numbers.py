"""Numbers for the SRL GAAA figure (outbox/2026-11-14-srl-gaaa). usage: python -I src/fig_srl_numbers.py OUTDIR
Carriers = the kill agent's final recount (results/lsu/bac/kill/work/recount.json, 'GAAA tetraloop (2661A list)',
every full-length 23S copy reads G2659 A2660 A2661 A2662). Loop strings are E. coli 2652-2667 from copies_kill.tsv
(full copies only). Writes phyla.tsv (GTDB phylum, carriers, metagenome-derived, cultured isolates) and loops.tsv (group, modal string)."""
import os
import collections, json, sys

import duckdb

W = os.environ.get("FROZEN_CORE", "frozen-core")  # lab working tree of the census
out = sys.argv[1]
full = json.load(open(f"{W}/results/lsu/bac/kill/work/recount.json"))["GAAA tetraloop (2661A list)"]["full"]
meta = {a: (t, c) for a, t, c in duckdb.sql(
    f"select acc, gtdb_taxonomy, ncbi_genome_category from '{W}/results/meta.parquet' "
    f"where acc in (select unnest(?::varchar[]))", params=[full]).fetchall()}
assert len(meta) == len(full) == 85, (len(meta), len(full))
# cultured isolates, by hand from NCBI BioSample records (fact-check 14 Nov): the other three genomes that are not
# "derived from metagenome" are Ca. Nesciobacter (a MAG by its BioSample), Ca. Gromoviella and the Pipizella
# endosymbiont (sequenced with their hosts)
ISOLATES = {"RS_GCF_000175255.2", "RS_GCF_000218875.1", "RS_GCF_943192945.1"}
assert ISOLATES <= set(full)
ph = collections.defaultdict(lambda: [0, 0, 0])
for a in full:
    p = meta[a][0].split(";")[1].removeprefix("p__")
    ph[p][0] += 1; ph[p][1] += meta[a][1] == "derived from metagenome"; ph[p][2] += a in ISOLATES
with open(f"{out}/phyla.tsv", "w") as f:
    f.write("phylum\tcarriers\tmetagenome_derived\tcultured_isolates\n")
    for p, (n, m, i) in sorted(ph.items(), key=lambda x: -x[1][0]):
        f.write(f"{p}\t{n}\t{m}\t{i}\n")
srl = collections.defaultdict(collections.Counter)
for a, s in duckdb.sql(f"select acc, C_SRL from read_csv('{W}/results/lsu/bac/kill/copies_kill.tsv', delim='\t', "
                       f"header=true, all_varchar=true) where fullC = 'True'").fetchall():
    srl[a][s] += 1
groups = [("E. coli (reference)", ["RS_GCF_000005845.2"])]
for lab, o in (("Arenicellales", "o__Arenicellales"), ("Christensenellales", "o__Christensenellales"),
               ("Aminicenantales", "o__Aminicenantales"), ("Zymomonas mobilis", "o__Sphingomonadales")):
    groups.append((lab, [a for a in full if f";{o};" in meta[a][0]]))
with open(f"{out}/loops.tsv", "w") as f:
    f.write("group\tgenomes\tmetagenome_derived\tmodal_2652_2667\tgenomes_with_modal\n")
    for lab, accs in groups:
        c = collections.Counter(srl[a].most_common(1)[0][0] for a in accs)
        s, n = c.most_common(1)[0]
        mg = sum(meta.get(a, ("", ""))[1] == "derived from metagenome" for a in accs)
        f.write(f"{lab}\t{len(accs)}\t{mg}\t{s}\t{n}\n")
print(open(f"{out}/phyla.tsv").read(), open(f"{out}/loops.tsv").read(), sep="\n")
