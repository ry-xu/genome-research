"""Prokaryote + virus subset of the 50% Swiss-Prot clusters, for speed. A cluster is kept if any member is from
Bacteria, Archaea or Viruses; its subject is the cluster rep if that is such a member, else the first such member.
-> data/sprot_p50.faa (then diamond makedb). Inputs: data/sprot_c50.tsv (rep, member), Swiss-Prot and NCBI
rankedlineage.dmp (src/paths.py; scripts/setup_reference_data.sh)."""
import gzip
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import RANKEDLINEAGE, ROOT, SPROT_FASTA  # noqa: E402

dom = {}
with open(RANKEDLINEAGE) as fh:
    for line in fh:
        f = [x.strip() for x in line.split("|")]
        dom[int(f[0])] = f[9]
keep_dom = {"Bacteria", "Archaea", "Viruses"}
tax, seqs = {}, {}
with gzip.open(SPROT_FASTA, "rt") as fh:
    sid = None
    for line in fh:
        if line.startswith(">"):
            sid = line[1:].split()[0]
            m = re.search(r" OX=(\d+)", line)
            tax[sid] = dom.get(int(m.group(1)), "") if m else ""
            seqs[sid] = []
        else:
            seqs[sid].append(line.strip())
clusters = {}
for line in open(ROOT / "data/sprot_c50.tsv"):
    rep, mem = line.rstrip("\n").split("\t")
    clusters.setdefault(rep, []).append(mem)
n = 0
with open(ROOT / "data/sprot_p50.faa", "w") as out:
    for rep, mems in clusters.items():
        ok = [m for m in mems if tax.get(m) in keep_dom]
        if not ok:
            continue
        s = rep if rep in ok else ok[0]
        out.write(f">{s}\n{''.join(seqs[s])}\n")
        n += 1
print(len(clusters), "clusters;", n, "kept;", sum(1 for v in tax.values() if v == ""), "entries with unknown domain")
