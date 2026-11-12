"""Which GTDB genomes go through Arm B (PREREG.md, "From Arm A to Arm B (budget)", Arm C extra entry, and the
Coriobacteriales lead).

  flag   Arm A flags: all if <= 2,500 genomes; else (1) codons flagged in >= 2 species of a genus or >= 3 of a family,
         then (2) every flagged genome not yet chosen, by most negative log10p over any flagged codon; at most 20
         genomes per genus at each step
  strat  one genome per GTDB family (most complete genes, ties by accession), 1,500 families drawn with seed 20261110
  armc   for a release-factor class the genome carries, no copy has the canonical motif (copies aligned over < 60% of
         the reference, or with a gap in the motif, are missing data); one genome per (genus, class, non-canonical
         motif residues), most complete genes, ties by accession (Changes 2026-11-10 in PREREG.md)
  actino every Actinomycetota family's stratified-rule genome (Coriobacteriales contrast; all Coriobacteriales
         representatives are run separately, data/corio/)

usage: uv run --no-project python -I src/armB_select.py
  -> results/gtdb/armB_select.tsv  acc (NCBI, no RS_/GB_), reasons, Arm B table (11 for a suspect table-4 call,
     else the Arm A calling table), taxonomy
"""
import collections
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
G = ROOT / "results" / "gtdb"
BUDGET, PER_GENUS, NFAM, SEED, COVER = 2500, 20, 1500, 20261110, 0.6


def rank(tax, r):
    for x in tax.split(";"):
        if x.startswith(r + "__"):
            return x
    return ""


def main():
    genomes = {}
    for line in open(G / "armA_genomes.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "acc":
            continue
        genomes[f[0]] = dict(table=11 if f[3] == "1" else int(f[2]), ncomplete=int(f[5]), flags=f[7], tax=f[8])
    reasons = collections.defaultdict(set)

    # flag
    flags = [line.rstrip("\n").split("\t") for line in open(G / "armA_flags.tsv")][1:]
    flagged = {f[0] for f in flags}
    if len(flagged) <= BUDGET:
        chosen = flagged
    else:
        best = {}
        for f in flags:  # acc, ..., codon (col 6), ..., log10p (col 13)
            best[f[0]] = min(best.get(f[0], 0.0), float(f[13]))
        by_codon = collections.defaultdict(set)
        for f in flags:
            by_codon[f[6]].add(f[0])
        recur = set()
        for c, accs in by_codon.items():
            gen = collections.Counter(rank(genomes[a]["tax"], "g") for a in accs)
            fam = collections.Counter(rank(genomes[a]["tax"], "f") for a in accs)
            recur |= {a for a in accs if gen[rank(genomes[a]["tax"], "g")] >= 2 or fam[rank(genomes[a]["tax"], "f")] >= 3}
        chosen = set()
        for step in (recur, flagged):
            per = collections.Counter()
            for a in sorted(step - chosen, key=lambda a: (best[a], a)):
                if len(chosen) >= BUDGET:
                    break
                g = rank(genomes[a]["tax"], "g")
                if per[g] < PER_GENUS:
                    chosen.add(a)
                    per[g] += 1
    for a in chosen:
        reasons[a].add("flag")

    # strat
    fam_pick = {}
    for a, v in genomes.items():
        fam = rank(v["tax"], "f")
        if fam in ("", "f__"):
            continue
        b = fam_pick.get(fam)
        if b is None or (-v["ncomplete"], a) < (-genomes[b]["ncomplete"], b):
            fam_pick[fam] = a
    fams = sorted(fam_pick)
    for fam in random.Random(SEED).sample(fams, min(NFAM, len(fams))):
        reasons[fam_pick[fam]].add("strat")

    # armc
    copies = collections.defaultdict(lambda: collections.defaultdict(dict))  # acc -> class -> gene -> state
    noncanon = collections.defaultdict(set)  # (acc, class, gene) -> {motif=residues} of its non-canonical motifs
    for line in open(G / "rf_motifs.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "acc":
            continue
        acc, gene, cls, motif, res, canon, cover = f[0], f[1], f[2], f[5], f[6], f[7] == "1", float(f[9])
        missing = cover < COVER or "-" in res
        st = copies[acc][cls].get(gene, "canon")
        if missing:
            st = "missing"
        elif not canon and st != "missing":
            st = "noncanon"
        if not canon and not missing:
            noncanon[acc, cls, gene].add(f"{motif}={res}")
        copies[acc][cls][gene] = st
    armc_pick = {}
    for acc, cl in copies.items():
        for cls, genes in cl.items():
            st = set(genes.values())
            if "canon" not in st and "noncanon" in st and acc in genomes:
                res = set().union(*(noncanon[acc, cls, g] for g, s in genes.items() if s == "noncanon"))
                key = (rank(genomes[acc]["tax"], "g"), cls, ",".join(sorted(res)))
                b = armc_pick.get(key)
                if b is None or (-genomes[acc]["ncomplete"], acc) < (-genomes[b]["ncomplete"], b):
                    armc_pick[key] = acc
    for (g, cls, _), acc in armc_pick.items():
        reasons[acc].add(f"armc:{cls}")

    # actino
    for fam, a in fam_pick.items():
        if rank(genomes[a]["tax"], "p") == "p__Actinomycetota":
            reasons[a].add("actino")

    with open(G / "armB_select.tsv", "w") as out:
        out.write("acc\tgtdb_acc\treasons\ttable\tncomplete\ttaxonomy\n")
        for a in sorted(reasons):
            v = genomes[a]
            out.write(f"{a[3:]}\t{a}\t{','.join(sorted(reasons[a]))}\t{v['table']}\t{v['ncomplete']}\t{v['tax']}\n")
    n = collections.Counter(r.split(":")[0] for v in reasons.values() for r in v)
    print(f"{len(flagged)} flagged genomes; selected {len(reasons)}: {dict(n)}; {len(fams)} families; "
          f"Arm C entries {len(armc_pick)}")


if __name__ == "__main__":
    main()
