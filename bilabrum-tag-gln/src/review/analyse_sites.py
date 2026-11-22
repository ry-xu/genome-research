"""Summaries from stops_in_hsps.py outputs.

usage: python3 -I analyse_sites.py PREFIX
"""
import math
import sys
from collections import Counter, defaultdict

B = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in B for b in B for c in B]
CODE = dict(zip(CODONS, AA))
STD = "ACDEFGHIKLMNPQRSTVWY"
STOPS = ("TAA", "TAG", "TGA")


def rows(path):
    with open(path) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(hdr, line.rstrip("\n").split("\t")))


def main(pre):
    sites = list(rows(pre + ".sites.tsv"))
    for s in sites:
        s["edge_dist"], s["win_ident"], s["qgaps6"] = int(s["edge_dist"]), int(s["win_ident"]), int(s["qgaps6"])
        s["clean"] = s["qres"] in STD and s["edge_dist"] >= 10 and s["qgaps6"] == 0 and s["win_ident"] >= 8
    print("== stop sites by class (distinct genomic sites, best-scoring HSP) ==")
    for c in STOPS:
        S = [s for s in sites if s["codon"] == c]
        al = [s for s in S if s["qres"] in STD]
        inner = [s for s in al if s["edge_dist"] >= 10]
        clean = [s for s in S if s["clean"]]
        print(f"{c}: all {len(S)} | aligned to a query residue {len(al)} | aligned & >=10 from HSP ends {len(inner)} "
              f"| clean (also no query gap within 6 cols, >=8/25 identical in window) {len(clean)} "
              f"in {len({s['query'] for s in clean})} queries, {len({s['contig'] for s in clean})} contigs")
        print(f"     clean residues: {Counter(s['qres'] for s in clean).most_common()}")
    # codon x residue matrix
    pair = defaultdict(Counter)
    for r in rows(pre + ".pairs.tsv"):
        if r["residue"] in STD:
            pair[r["codon"]][r["residue"]] += int(r["n"])
    tot = Counter()
    for c in pair:
        if c not in STOPS:
            tot.update(pair[c])
    N = sum(tot.values())
    comp = {a: (tot[a] + 0.5) / (N + 10) for a in STD}
    prof = {}
    for a in STD:
        v = Counter()
        for c in CODONS:
            if CODE[c] == a:
                v.update(pair[c])
        n = sum(v.values())
        prof[a] = {b: (v[b] + 0.5) / (n + 10) for b in STD}
    own = {a: sum(prof[a][a] for _ in [0]) for a in STD}
    for c in STOPS:
        res = Counter(s["qres"] for s in sites if s["codon"] == c and s["clean"])
        n = sum(res.values())
        if n == 0:
            continue
        ll = {a: sum(k * math.log(prof[a][b] / comp[b]) for b, k in res.items()) / n for a in STD}
        best = sorted(ll, key=lambda a: -ll[a])
        print(f"{c} clean n={n}: best-fitting codon profile {best[0]} (llr/col {ll[best[0]]:.3f}), next {best[1]} "
              f"({ll[best[1]]:.3f}), {best[2]} ({ll[best[2]]:.3f}); share Q {res['Q'] / n:.3f} vs Q codons' Q share "
              f"{prof['Q']['Q']:.3f}, background Q {comp['Q']:.3f}")
    for c in ("CAG", "CAA", "AAG", "GAG", "TGG"):
        v = pair[c]
        n = sum(v.values())
        print(f"   {c} aligned n={n}: " + ", ".join(f"{a} {v[a] / n:.3f}" for a, _ in v.most_common(5)))
    # single-substitution error model: if TAG sites were point errors at rate r per base, TAA/TGA should appear too.
    cod_n = {c: sum(pair[c].values()) for c in CODONS}
    def neighbours(t):
        out = []
        for i in range(3):
            for b in "ACGT":
                if b != t[i]:
                    out.append(t[:i] + b + t[i + 1:])
        return [x for x in out if x not in STOPS]
    src = {t: sum(cod_n[x] for x in neighbours(t)) for t in STOPS}
    obs = {t: sum(pair[t].values()) for t in STOPS}
    r = obs["TAG"] / src["TAG"]
    print("== uniform point-error model (aligned columns, all HSP positions) ==")
    for t in STOPS:
        print(f"   {t}: sense codons one step away {src[t]}, observed aligned {obs[t]}, expected at TAG's rate {r * src[t]:.1f}")
    # C->T at first position only (deamination-like): CAG->TAG, CAA->TAA, CGA->TGA
    r2 = obs["TAG"] / cod_n["CAG"]
    print(f"   C->T first-position model: rate {r2:.4f}; expected TAA from CAA {r2 * cod_n['CAA']:.1f} (obs {obs['TAA']}), "
          f"TGA from CGA {r2 * cod_n['CGA']:.1f} (obs {obs['TGA']}); CAG n={cod_n['CAG']} CAA n={cod_n['CAA']} CGA n={cod_n['CGA']}")
    # HSPs: GC3 and chromosomes with vs without TAG
    H = list(rows(pre + ".hsps.tsv"))
    clean_by_hsp = Counter((s["contig"], s["query"]) for s in sites if s["codon"] == "TAG" and s["clean"])
    w, wo = [], []
    for h in H:
        if int(h["ncod"]) < 60:
            continue
        (w if clean_by_hsp[(h["contig"], h["query"])] else wo).append(h)
    def mean(v):
        return sum(v) / len(v) if v else float("nan")
    print("== HSPs (>=60 codons, owning >=1 codon) with vs without a clean TAG ==")
    print(f"   with: n={len(w)} GC3 mean {mean([float(h['gc3']) for h in w]):.3f}, pid mean {mean([float(h['pid']) for h in w]):.1f}")
    print(f"   without: n={len(wo)} GC3 mean {mean([float(h['gc3']) for h in wo]):.3f}, pid mean {mean([float(h['pid']) for h in wo]):.1f}")
    cw = Counter(h["contig"] for h in w)
    cwo = Counter(h["contig"] for h in wo)
    both = [c for c in cw if cwo[c]]
    print(f"   contigs with TAG HSPs: {len(cw)}; of these also carrying TAG-free HSPs: {len(both)}; "
          f"unplaced (CBCXTV) contigs among TAG HSPs: {sum(c.startswith('CBCXTV') for c in cw)}; "
          f"among TAG-free: {sum(c.startswith('CBCXTV') for c in cwo)}")
    # gene ends
    E = list(rows(pre + ".ends.tsv"))
    near = [e for e in E if e["first_stop_dist"] not in ("None",) and int(e["first_stop_dist"]) <= 20]
    print("== gene ends: HSPs reaching query C-terminus (qend >= qlen-3), first in-frame stop within 20 codons ==")
    print(f"   loci {len(E)}; with a stop within 20 codons {len(near)}: {dict(Counter(e['first_stop'] for e in near))}")
    near2 = [e for e in E if e["first_TAA_TGA_dist"] != "None" and int(e["first_TAA_TGA_dist"]) <= 20]
    print(f"   first TAA/TGA within 20 codons: {len(near2)}; TAGs passed before it: "
          f"{sum(int(e['tags_before']) for e in near2)}")


if __name__ == "__main__":
    main(sys.argv[1])
