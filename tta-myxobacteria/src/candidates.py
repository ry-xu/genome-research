"""GTDB Arm B candidates (results/gtdb/armB_passes.tsv): check 2 of PREREG.md "From an Arm B pass to a result" and
the facts the other checks need.

Per candidate (genome, X):
  synonym contrast  X's ratio / pooled ratio (sum O / sum E) of X's synonymous codons in the same genome (check 2:
                    >= 3 passes; not applicable to ATG and TGG under table 11)
  offsets           where X sits in the window among the genome's events (e*-1, e*, e*+1, e*+2, e*+3; a stop at the
                    homolog's end is expected at e*+1)
  genus             other GTDB genomes of the genus that went through Arm B, and how many pass X (check 1 counts any
                    other assembly of the genus that passes X)
  metadata          CheckM2 completeness / contamination, contigs, N50, GTDB GC, NCBI category, BioProject,
                    BioSample, submitter, isolation source

usage: uv run --no-project python -I src/candidates.py  -> results/gtdb/candidates.tsv
"""
import collections
import csv
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
G = ROOT / "results" / "gtdb"
META = ROOT / "data" / "gtdb_meta"
NA = ("ATG", "TGG")
MCOLS = ("checkm2_completeness", "checkm2_contamination", "contig_count", "n50_contigs", "gc_percentage",
         "ncbi_genome_category", "ncbi_bioproject", "ncbi_biosample", "ncbi_submitter", "ncbi_isolation_source")


def tsv(fn):
    return list(csv.DictReader(open(fn), delimiter="\t"))


def meta(accs):
    out = {}
    for fn in ("ar53_metadata.tsv.gz", "bac120_metadata.tsv.gz"):
        with gzip.open(META / fn, "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["accession"] in accs:
                    out[r["accession"]] = r
    return out


def main():
    passes = tsv(G / "armB_passes.tsv")
    label_of = {r["gtdb_acc"]: r["label"] for r in tsv(ROOT / "data" / "gtdb_armB" / "labels.tsv")}
    sel = {r["gtdb_acc"]: r for r in tsv(G / "armB_select.tsv")}
    rows = collections.defaultdict(dict)
    events = collections.defaultdict(list)
    for i in (0, 1):
        for r in tsv(G / f"armB_{i}.tsv"):
            rows[r["label"]][r["codon"]] = r
        for r in tsv(G / f"armB_{i}_events.tsv"):
            events[r["label"]].append(r["window"].split(","))
    md = meta({p["gtdb_acc"] for p in passes})
    genus = lambda tax: next(x for x in tax.split(";") if x.startswith("g__"))
    by_genus = collections.defaultdict(list)
    for a, s in sel.items():
        by_genus[genus(s["taxonomy"])].append(a)
    passing = {(p["gtdb_acc"], p["codon"]) for p in passes}
    out = []
    for p in passes:
        a, x, lab = p["gtdb_acc"], p["codon"], p["label"]
        d = rows[lab]
        sib = [k for k in d if k != x and d[k]["aa"] == d[x]["aa"]]
        so, se = sum(int(d[k]["O"]) for k in sib), sum(float(d[k]["E"]) for k in sib)
        sr = so / se if se else float("nan")
        contrast = "n/a" if (x in NA and p["table"] == "11") else f"{float(p['ratio']) / sr:.2f}" if sr else "inf"
        off = [sum(w[i] == x for w in events[lab]) for i in range(5)]
        g = genus(p["taxonomy"])
        others = [b for b in by_genus[g] if b != a]
        m = md.get(a, {})
        out.append(dict(gtdb_acc=a, label=lab, codon=x, aa=p["aa"], O=p["O"], O_fam=p["O_fam"], ratio=p["ratio"],
                        log10p=p["log10p"], syn_O=so, syn_E=f"{se:.2f}", syn_ratio=f"{sr:.3f}", syn_contrast=contrast,
                        offsets="/".join(map(str, off)), genus=g, genus_others_run=len(others),
                        genus_others_pass=sum((b, x) in passing for b in others),
                        **{k: m.get(k, "") for k in MCOLS}, taxonomy=p["taxonomy"]))
    with open(G / "candidates.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    for r in sorted(out, key=lambda r: (r["codon"], r["taxonomy"])):
        print(r["gtdb_acc"], r["codon"], f"ratio {r['ratio']} syn {r['syn_ratio']} contrast {r['syn_contrast']}",
              f"offsets {r['offsets']}", r["genus"], f"genus {r['genus_others_pass']}/{r['genus_others_run']}",
              f"cm2 {r['checkm2_completeness']}/{r['checkm2_contamination']} contigs {r['contig_count']}",
              f"N50 {r['n50_contigs']} GC {r['gc_percentage'][:4]}", r["ncbi_genome_category"][:12],
              r["ncbi_bioproject"], r["ncbi_submitter"][:30], "|", r["ncbi_isolation_source"][:40], sep="  ")


if __name__ == "__main__":
    main()
