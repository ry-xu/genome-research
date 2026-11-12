"""Natural controls for Arm A: call a real genome's genes with a fixed translation table and profile them.
Calling an ordinary genome with table 4 (TGA = Trp) or 16 (TAG = Leu) makes a real stop look like a sense codon to
the method, with real gene ends and real read-on tails: the method must flag it. Calling a genome that really reads
TAG (or TGA) as sense with table 16 (or 4) must not flag it.

usage: python src/controls.py <out.npz> <fna.gz>:<table>[,<table>...] ...   (table 'g' = GTDB's 11-or-4 rule)
"""
import sys

import numpy as np
import pyrodigal

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import codonpos as cp
from inject import call, genes_nt, read_fasta


def call_table(contigs, table):
    if table == "g":
        return call(contigs), None
    table = int(table)
    # single mode only: meta mode picks its own models and cannot be held to one table
    assert sum(len(s) for _, s in contigs) >= 100000, "genome too small for a fixed-table call"
    finder = pyrodigal.GeneFinder()
    finder.train(*[s.encode() for _, s in contigs], translation_table=table)
    found = [(n, finder.find_genes(s.encode())) for n, s in contigs]
    out = [(n, g.begin, g.end, g.strand, g.partial_begin, g.partial_end) for n, gs in found for g in gs]
    return out, table


def main(out, args):
    res = {}
    for a in args:
        fn, tabs = a.rsplit(":", 1)
        contigs = list(read_fasta(fn))
        for t in tabs.split(","):
            calls, table = call_table(contigs, t)
            stops = cp.TABLE_STOPS[table] if table else cp.STOPS
            prof = cp.profile(genes_nt(contigs, calls), stops_used=stops)
            key = f"{fn.rsplit('/', 1)[-1].replace('.fna.gz', '')}|{t}"
            res[key] = prof
            print(key, prof["ncomplete"], prof["stops"], flush=True)
    np.savez_compressed(out, keys=np.array(list(res)), meta=np.array([[p["gc"], p["ncomplete"]] for p in res.values()]),
                        **{f"pos|{i}": p["pos"] for i, p in enumerate(res.values())})


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
