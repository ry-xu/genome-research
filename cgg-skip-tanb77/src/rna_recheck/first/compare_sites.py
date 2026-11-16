"""Re-check, step 1b: compare the rebuilt sites with the original results/rna_cgg/sites.tsv (and sites_v1.tsv).
  .venv/bin/python -I compare_sites.py      (run after unique_probes.py)
Rebuilt sets: candidates.tsv minus the sites my own Amendment 3 check drops (status ncrna in sites_recheck.tsv);
that drop is applied to A and C (the sets counted); A0 is compared unfiltered and is not used.
Key: (set, genome, contig, 0-based C position, strand); sites.tsv cpos/run_start are 0-based, run_end 0-based
half-open; candidates.tsv cpos1 is 1-based, run_start0/run_end0 0-based half-open. Compares codon, run and CDS.
"""
import csv
import os
from collections import Counter

R = "results/rna_cgg/"
nc = {(r["set"], r["site_id"]) for r in csv.DictReader(open(R + "recheck/sites_recheck.tsv"), delimiter="\t")
      if r["status"] == "ncrna"}
mine_all, mine = {}, {}
for r in csv.DictReader(open(R + "recheck/candidates.tsv"), delimiter="\t"):
    k = (r["set"], r["genome"], r["contig"], int(r["cpos1"]) - 1, r["strand"])
    v = (r["codon"], int(r["run_start0"]), int(r["run_end0"]), r["cds"])
    mine_all[k] = v
    if (r["set"], f'{r["genome"]}:{r["contig"]}:{r["cpos1"]}:{r["strand"]}') not in nc:
        mine[k] = v
print("dropped by my Amendment 3 check:", Counter(s for s, _ in nc))


def load(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        out[(r["set"], r["genome"], r["contig"], int(r["cpos"]), r["strand"])] = (
            r["codon"], int(r["run_start"]), int(r["run_end"]), r["cds"])
    return out


for label, path, m in (("sites.tsv (Amendment 3)", R + "sites.tsv", mine),
                       ("sites_v1.tsv (before Amendment 3)", R + "sites_v1.tsv", mine_all)):
    if not os.path.exists(path):
        print(label, "missing")
        continue
    orig = load(path)
    print("==", label, "vs", "rebuilt after my ncRNA drop" if m is mine else "rebuilt before it")
    for s in ("A", "C", "A0"):
        a = {k for k in m if k[0] == s}
        o = {k for k in orig if k[0] == s}
        both = a & o
        diff = Counter(f for k in both
                       for f, x, y in zip(("codon", "run_start", "run_end", "cds"), m[k], orig[k]) if x != y)
        print(f"set {s}: rebuilt {len(a)}, original {len(o)}, shared {len(both)}, only rebuilt {len(a - o)}, "
              f"only original {len(o - a)}, field differences {dict(diff) or 0}")
