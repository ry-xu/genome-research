"""Re-check, step 1: candidate Set A and Set C sites from the rule's Sites section, built independently.
Run with the atlas python (proteomics_db imports codon_atlas):
  python -I build_sites.py
Every catalogue CDS of the 24 genomes is walked from its start with walk(s, "skip") of src/proteomics_db.py.
Set A: the C of every in-frame CGG met by a walk (primary genomes). Set C: the C of every in-frame CGA/CGC/CGT met by
the same walks (the walk's own reading frame, i.e. the +1 frame after each skip). Positions deduplicated per genome;
sites with the C within 30 nt of a contig end dropped (index < 30 from either end).
Site run: maximal run of identical plus-strand bases containing the C position.
Probes (plus strand): keep = left20 + run + right20; lack = left20 + run minus one base + right20.
Also writes the non-primary genomes' CGG sites (set A0) only so the probe-uniqueness step sees all 24 genomes;
they are not counted.
Output: candidates.tsv (one row per site, before the probe-uniqueness filter)
Variants recorded for the comparison with sites.tsv:
  walk_has_cgg  1 if at least one walk meeting this position has a CGG (Set C restricted to walks with a skip)
  best_walk     1 if the position is met by the longest walk among walks sharing an end (proteomics_db's dedup)
"""
import csv
import sys
from pathlib import Path

ROOT = Path(".")
sys.path.insert(0, str(ROOT / "src"))
from proteomics_db import gff, rc, read_fasta, walk  # noqa: E402

MDIR = ROOT / "data/proteomics/mgnify_pig"
OUT = ROOT / "results/rna_cgg/recheck"
FLANK, EDGE = 20, 30
CTRL = {"CGA", "CGC", "CGT"}


def main():
    rows = list(csv.DictReader(open(ROOT / "results/proteomics_pig/genomes.tsv"), delimiter="\t"))
    primary = {r["genome"]: int(r["primary"]) for r in rows}
    assert len(primary) == 24, len(primary)
    print("primary genomes:", sum(primary.values()), file=sys.stderr)
    out = []
    stats = {}
    for g in sorted(primary):
        seqs = read_fasta(MDIR / f"{g}.fna")
        cds = gff(MDIR / f"{g}.gff")
        found = {}  # (contig, pos1, strand) -> dict
        best = {}  # walk end -> (len, cid)
        walks = []
        for c, a, b, st, cid in cds:
            x = seqs[c]
            s = x[a - 1:] if st == "+" else rc(x[:b])
            prot, used, sites, complete, pos = walk(s, "skip")
            end = (c, st, a - 1 + used if st == "+" else b - used) + (() if complete else (cid,))
            if end not in best or len(prot) > best[end][0]:
                best[end] = (len(prot), cid)
            hits = []
            for _, i in sites:
                assert s[i:i + 3] == "CGG"
                hits.append((i, "A", "CGG"))
            for r, p in enumerate(pos):
                cod = s[p:p + 3]
                if cod in CTRL:
                    hits.append((p, "C", cod))
            walks.append((c, a, b, st, cid, end, bool(sites), hits))
        best_cids = {v[1] for v in best.values()}
        for c, a, b, st, cid, end, has_cgg, hits in walks:
            for i, sset, cod in hits:
                p1 = a + i if st == "+" else b - i
                key = (c, p1, st)
                d = found.setdefault(key, {"set": sset, "codon": cod, "cds": cid, "has_cgg": 0, "best": 0})
                assert d["set"] == sset and d["codon"] == cod, (g, key)
                d["has_cgg"] |= has_cgg
                d["best"] |= cid in best_cids
        n_edge = 0
        for (c, p1, st), d in sorted(found.items()):
            x = seqs[c]
            p0 = p1 - 1
            base = x[p0]
            assert base == ("C" if st == "+" else "G"), (g, c, p1, st, base)
            if p0 < EDGE or len(x) - 1 - p0 < EDGE:
                n_edge += 1
                continue
            rs, re = p0, p0 + 1
            while rs > 0 and x[rs - 1] == base:
                rs -= 1
            while re < len(x) and x[re] == base:
                re += 1
            if rs < FLANK or re + FLANK > len(x):
                n_edge += 1
                continue
            left, run, right = x[rs - FLANK:rs], x[rs:re], x[re:re + FLANK]
            sset = d["set"] if primary[g] else ("A0" if d["set"] == "A" else "C0")
            out.append([sset, g, c, p1, st, d["codon"], rs, re, re - rs, d["cds"], d["has_cgg"], d["best"],
                        left + run + right, left + run[:-1] + right])
        stats[g] = (len(found), n_edge)
        print(g, primary[g], len(cds), "positions", len(found), "edge-dropped", n_edge, file=sys.stderr, flush=True)
    with open(OUT / "candidates.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["set", "genome", "contig", "cpos1", "strand", "codon", "run_start0", "run_end0", "run_len", "cds",
                    "walk_has_cgg", "best_walk", "keep_probe", "lack_probe"])
        w.writerows(out)
    from collections import Counter
    print(Counter(r[0] for r in out), file=sys.stderr)


if __name__ == "__main__":
    main()
