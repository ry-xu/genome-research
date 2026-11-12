"""GTDB members of the TTA clade outside the frozen Arm B labels: the frozen Arm B rule and check 2 on them.

All GTDB R232 genomes of the clade's five genera (g__JADLFU01, g__JANWXX01, g__DASTPP01, g__JAXLYQ01, g__CAYNEK01) are
the 8 species representatives plus two non-representatives, GCA_963169955.1 (JADLFU01 sp020849655) and
GCA_033578855.1 (JANWXX01 sp025057345). CAYNEK01 GCA_964396175.1 was a representative that never got Arm B labels.
Arm B on these three: results/clade_members/armB*.tsv (src/term_codon.py, frozen, table "g").

Same rules as src/check1_eval.py (Arm B pass) and src/candidates.py (synonym contrast, check 2 >= 3).

usage: uv run --no-project python -I src/clade_members.py
  -> results/clade_members/members.tsv (all three), candidates.tsv (Arm B passes; read by reads_check.py)
"""
import csv
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
M = ROOT / "results" / "clade_members"
META = ROOT / "data" / "gtdb_meta"
X = "TTA"
MCOLS = ("checkm2_completeness", "checkm2_contamination", "contig_count", "ncbi_bioproject", "ncbi_biosample",
         "ncbi_isolation_source", "gtdb_representative", "gtdb_taxonomy")


def tsv(fn):
    return list(csv.DictReader(open(fn), delimiter="\t"))


def ok(r):  # src/check1_eval.py
    return int(r["O_fam"]) >= 10 and float(r["ratio"]) >= 3 and float(r["log10p"]) <= -10


def main():
    rows = {}
    for r in tsv(M / "armB.tsv"):
        rows.setdefault(r["label"], {})[r["codon"]] = r
    accs = {"GB_" + lab for lab in rows}
    md = {}
    with gzip.open(META / "bac120_metadata.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["accession"] in accs:
                md[r["accession"]] = r
    out = []
    for lab, d in rows.items():
        a, r = "GB_" + lab, d[X]
        sib = [k for k in d if k != X and d[k]["aa"] == d[X]["aa"]]
        so, se = sum(int(d[k]["O"]) for k in sib), sum(float(d[k]["E"]) for k in sib)
        sr = so / se
        m = md[a]
        out.append(dict(gtdb_acc=a, label=lab, codon=X, aa=r["aa"], O=r["O"], O_fam=r["O_fam"], ratio=r["ratio"],
                        log10p=r["log10p"], rank=r["rank"], armB_pass="yes" if ok(r) else "no", syn_O=so,
                        syn_E=f"{se:.2f}", syn_ratio=f"{sr:.3f}", syn_contrast=f"{float(r['ratio']) / sr:.2f}",
                        genus=next(x for x in m["gtdb_taxonomy"].split(";") if x.startswith("g__")),
                        **{k: m[k] for k in MCOLS}))
    for fn, keep in (("members.tsv", out), ("candidates.tsv", [r for r in out if r["armB_pass"] == "yes"])):
        with open(M / fn, "w") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
            w.writeheader()
            w.writerows(keep)
    for r in out:
        print(r["gtdb_acc"], r["genus"], f"O_fam {r['O_fam']} ratio {r['ratio']} log10p {r['log10p']} rank {r['rank']}",
              f"ArmB {r['armB_pass']}", f"syn {r['syn_ratio']} contrast {r['syn_contrast']}",
              f"cm2 {r['checkm2_completeness']}/{r['checkm2_contamination']}", f"rep {r['gtdb_representative']}",
              r["ncbi_biosample"], r["ncbi_isolation_source"][:50], sep="  ")


if __name__ == "__main__":
    main()
