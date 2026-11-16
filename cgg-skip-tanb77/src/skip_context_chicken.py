"""Replication of the CUU-before-skip count (rule: results/skip_context_chicken/README.md, fixed 15 Nov before any
upstream codon was read). usage (from the folder root): uv run python -I src/skip_context_chicken.py
inputs: results/rna_chicken/sites.tsv and site_class.tsv (shipped as data/cuu_context/chicken_sites_setA.tsv.gz and
chicken_site_class.tsv.gz; gunzip them into results/rna_chicken/), genomes and pyrodigal GFFs in data/rna_chicken/genomes/"""
import csv
import random
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "rna_chicken"
GDIR = ROOT / "data" / "rna_chicken" / "genomes"
OUT = ROOT / "results" / "skip_context_chicken"
COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"TAA", "TAG", "TGA"}
ORIG_UP = ["CTT", "TTG", "TCA", "CCA", "GAC"]
ORIG_DOWN = ["TTT", "AAA", "GAA", "TAT"]


def rc(s):
    return s.translate(COMP)[::-1]


def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def background(acc, seqs):
    """codons immediately preceding an in-frame GGN codon in the genome's CDSs"""
    c = Counter()
    for line in open(GDIR / f"{acc}.gff"):
        f = line.rstrip("\n").split("\t")
        if line[0] == "#" or len(f) < 9 or f[2] != "CDS":
            continue
        s = seqs[f[0]][int(f[3]) - 1:int(f[4])]
        if f[6] == "-":
            s = rc(s)
        cod = [s[i:i + 3] for i in range(0, len(s) - 2, 3)]
        for a, b in zip(cod, cod[1:]):
            if b.startswith("GG") and a not in STOPS and "N" not in a + b:
                c[a] += 1
    return c


def main():
    excl = {r["gca"] for r in csv.DictReader(open(ROOT / "data" / "frozen_test" / "cgg_test_T.tsv"), delimiter="\t")}
    cls = {r["site_id"]: r["class"] for r in csv.DictReader(open(RES / "site_class.tsv"), delimiter="\t")}
    sites = [r for r in csv.DictReader(open(RES / "sites.tsv"), delimiter="\t") if r["set"] == "A"]
    genomes = sorted({r["genome"] for r in sites})
    kept = [g for g in genomes if g not in excl]
    lines = [f"Set A genomes {len(genomes)}; excluded (in the original count's set) {len(genomes) - len(kept)}; "
             f"kept {len(kept)}"]
    obs = {k: defaultdict(Counter) for k in ("skip", "inframe")}
    nsite = {k: Counter() for k in ("skip", "inframe")}
    bg = {}
    conv_ok = Counter()
    for g in kept:
        seqs = fasta(GDIR / f"{g}.fna")
        bg[g] = background(g, seqs)
        for r in sites:
            if r["genome"] != g or cls.get(r["site_id"]) not in ("skip", "inframe"):
                continue
            s, p = seqs[r["contig"]], int(r["cpos"])
            # cpos is the 0-based C (checked on 3,000 sites before this run: plus codon s[p:p+3]; minus codon
            # rc(s[p-2:p+1]), so the codon 5' of it on the minus strand is rc(s[p+1:p+4]))
            if r["strand"] == "+":
                codon, up = s[p:p + 3], s[p - 3:p]
            else:
                codon, up = rc(s[p - 2:p + 1]), rc(s[p + 1:p + 4])
            conv_ok[codon == r["codon"]] += 1
            if codon != r["codon"] or len(up) != 3:
                continue
            k = cls[r["site_id"]]
            obs[k][g][up] += 1
            nsite[k][g] += 1
    lines.append(f"codon convention check: {dict(conv_ok)} (True = site codon reproduced)")

    def ratio(k, gs, cod):
        o = sum(obs[k][g][cod] for g in gs)
        e = sum(nsite[k][g] * bg[g][cod] / sum(bg[g].values()) for g in gs if sum(bg[g].values()))
        return o, e, (o / e if e else float("nan"))
    for k in ("skip", "inframe"):
        n = sum(nsite[k].values())
        lines.append(f"{k}-like sites: {n} in {sum(1 for g in kept if nsite[k][g])} genomes")
    tab = {}
    allc = sorted({c for g in kept for c in bg[g]})
    with open(OUT / "codon_ratios.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["codon", "skip_obs", "skip_exp", "skip_ratio", "inframe_obs", "inframe_exp", "inframe_ratio"])
        for c in allc:
            a, b = ratio("skip", kept, c), ratio("inframe", kept, c)
            tab[c] = (a, b)
            w.writerow([c, a[0], f"{a[1]:.1f}", f"{a[2]:.3f}", b[0], f"{b[1]:.1f}", f"{b[2]:.3f}"])
    for c in ORIG_UP + ORIG_DOWN:
        (o, e, r), (o2, e2, r2) = tab[c]
        lines.append(f"  {c}: skip-like {o} vs {e:.1f} expected (ratio {r:.2f}); in-frame-like {o2} vs {e2:.1f} ({r2:.2f})")
    top = sorted(((tab[c][0][2], c, tab[c][0][0], tab[c][0][1]) for c in allc if tab[c][0][1] >= 5), reverse=True)
    lines.append("top skip-like ratios (expected >= 5): " + ", ".join(f"{c} {r:.2f} ({o}/{e:.1f})" for r, c, o, e in top[:8]))
    lines.append("bottom: " + ", ".join(f"{c} {r:.2f} ({o}/{e:.1f})" for r, c, o, e in top[-6:]))
    rng = random.Random(1)
    gs = [g for g in kept if nsite["skip"][g]]
    boots = []
    for _ in range(1000):
        smp = [rng.choice(gs) for _ in gs]
        o = sum(obs["skip"][g]["CTT"] for g in smp)
        e = sum(nsite["skip"][g] * bg[g]["CTT"] / sum(bg[g].values()) for g in smp)
        boots.append(o / e if e else float("nan"))
    boots.sort()
    lines.append(f"CTT skip-like ratio bootstrap over genomes (1,000, seed 1): 95% {boots[25]:.2f}-{boots[974]:.2f}")
    r_skip, r_in = tab["CTT"][0][2], tab["CTT"][1][2]
    p1, p2 = r_skip >= 2, r_skip > r_in
    p3 = sum(tab[c][0][2] > 1 for c in ORIG_UP) >= 4
    lines.append(f"P1 (CTT >= 2): {p1}; P2 (skip > in-frame): {p2}; P3 (>= 4 of 5 > 1): {p3}")
    lines.append("decision: " + ("replicated" if p1 and p2 and p3 else ("not replicated" if not p1 else "partial")))
    (OUT / "result.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
