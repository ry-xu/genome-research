"""Compile the per-genome CGG table: union of in-frame CGG loci (03, *.cgg_loci.tsv) and CGG at frameshift events
(04, *.fs_cgg_sites.tsv), keyed by (contig, strand, codon start), with repeated copies collapsed by +-30 nt context.
Frame behaviour per unique CGG: +1 / -1 if it is a 04 site of that direction or a 03 locus with an event within 9 nt,
else same. Residue: 04 consensus opposite the CGG+N unit for shifted sites, 03 consensus opposite CGG otherwise.
Null codons use the 03 statistics from *.summary.json (same definitions as the CGG in-frame loci)."""
import sys, json, ast
from collections import Counter
sys.path.insert(0, "src/recheck")
from mp import read_fasta, rc
W = "recheck_work"

def rows(path):
    L = open(path).read().splitlines()
    h = L[0].split("\t")
    return [dict(zip(h, l.split("\t"))) for l in L[1:]]

out = []
for lab in sys.argv[1:]:
    g = read_fasta(f"{W}/genomes/{lab}.fna")
    def ctx(c, st, s):
        x = g[c][max(0, s - 30):s + 33]
        return x if st == "+" else rc(x)
    U = {}
    for r in rows(f"{W}/out/{lab}.cgg_loci.tsv"):
        k = (r["contig"], r["strand"], int(r["start0"]))
        fs = ast.literal_eval(r["fs_near"]) if r["fs_near"] else {}
        fr = "+1" if fs.get("+1") else ("-1" if fs.get("-1") else "same")
        U[k] = dict(frame=fr, res=r["consensus"], clean=r["clean_residues"], src="03")
    for r in rows(f"{W}/out/{lab}.fs_cgg_sites.tsv"):
        k = (r["contig"], r["strand"], int(r["cgg_start0"]))
        U[k] = dict(frame=r["dir"], res=r["consensus"], clean=r["clean_residues"], src="04")
    seen = set(); n = Counter(); res = {"same": Counter(), "+1": Counter(), "-1": Counter()}; cl = {"same": Counter(), "+1": Counter(), "-1": Counter()}
    for k in sorted(U):
        c = ctx(*k)
        if c in seen: continue
        seen.add(c)
        v = U[k]; n[v["frame"]] += 1
        if v["res"]: res[v["frame"]][v["res"]] += 1
        for x in v["clean"]: cl[v["frame"]][x] += 1
    S = json.load(open(f"{W}/out/{lab}.summary.json"))
    def fmt(c, k=4): return " ".join(f"{a}{b}" for a, b in c.most_common(k)) or "none"
    def null(cod):
        d = S["codons"][cod] if "codons" in S else S[cod]
        fs = d["fs"]; p, m = fs.get("+1", 0), fs.get("-1", 0)
        return f"{cod} {d['n_loci']}: {d['n_loci']-p-m}/{p}/{m}; {fmt(Counter(d['cons']),2)}; junc {d['junction_old']}"
    out.append((lab, sum(n.values()), n, res, cl, [null(c) for c in ("CGA", "CGC", "AGG")],
                {c: (S["codons"][c] if "codons" in S else S[c])["junction_old"] for c in ("CGG",)}))
    print(lab, "unique CGG", sum(n.values()), dict(n), "| res same:", fmt(res["same"], 6), "| res +1:", fmt(res["+1"], 6), "| res -1:", fmt(res["-1"], 6),
          "| clean same:", fmt(cl["same"], 6), "| clean +1:", fmt(cl["+1"], 6), "|", " | ".join(out[-1][5]))
