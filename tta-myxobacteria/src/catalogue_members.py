"""Catalogue candidates of the TTA clade (outside GTDB R232): the frozen Arm B rule and check 2 on them.

Candidates come from the marker-AAI screen of NCBI and MAG catalogues (results/catalogues/members.tsv, 2026-11-11).
Arm B on them: results/clade_members/armB_cat*.tsv (src/term_codon.py, frozen, table "g"); the redundant copy
MGYG000501200 (= GCA_035762175.1) was not rerun.

Same rules as src/clade_members.py (Arm B pass as src/check1_eval.py; synonym contrast as src/candidates.py, check 2 >= 3).

usage: uv run --no-project python -I src/catalogue_members.py
  -> results/clade_members/cat_members.tsv (all), cat_candidates.tsv (Arm B passes; read by reads_check.py)
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from clade_members import M, X, ok, tsv  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CAT = ROOT / "results" / "catalogues" / "members.tsv"
MCOLS = ("source", "local_path", "completeness", "contamination", "cc_source", "environment", "biosample")


def label_of(path):  # as src/term_codon.py labels a genome
    return Path(path).name.split(".fna")[0].split(".fa")[0]


def main():
    rows, evf = {}, {}
    for p in sorted(M.glob("armB_cat*.tsv")):  # armB_cat.tsv (10:4x), armB_cat2.tsv (NCBI resume, 12:34)
        if p.name.endswith("_events.tsv"):
            continue
        for r in tsv(p):
            assert evf.setdefault(r["label"], p.stem) == p.stem, f"{r['label']} in two Arm B runs"
            rows.setdefault(r["label"], {})[r["codon"]] = r
    cat = tsv(CAT)
    labs, paths = [label_of(r["local_path"]) for r in cat], [r["local_path"] for r in cat]
    # outputs and data/reads/<label>/ are named by label alone (audit 2026-11-11 12:24)
    assert len(set(labs)) == len(labs) and len(set(paths)) == len(paths), "duplicate catalogue label or path"
    md = dict(zip(labs, cat))
    out = []
    for lab, d in rows.items():
        r, m = d[X], md[lab]
        sib = [k for k in d if k != X and d[k]["aa"] == d[X]["aa"]]
        so, se = sum(int(d[k]["O"]) for k in sib), sum(float(d[k]["E"]) for k in sib)
        sr = so / se
        near = m["screen"].split(" (")[1].split(",")[0] if " (" in m["screen"] else ""
        out.append(dict(gtdb_acc=lab, label=lab, codon=X, aa=r["aa"], O=r["O"], O_fam=r["O_fam"], ratio=r["ratio"],
                        log10p=r["log10p"], rank=r["rank"], armB_pass="yes" if ok(r) else "no", syn_O=so,
                        syn_E=f"{se:.2f}", syn_ratio=f"{sr:.3f}", syn_contrast=f"{float(r['ratio']) / sr:.2f}",
                        nearest_genus=near, events=f"{evf[lab]}_events.tsv", **{k: m[k] for k in MCOLS}))
    for fn, keep in (("cat_members.tsv", out), ("cat_candidates.tsv", [r for r in out if r["armB_pass"] == "yes"])):
        with open(M / fn, "w") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
            w.writeheader()
            w.writerows(keep)
    for r in out:
        print(r["label"], r["nearest_genus"], f"O_fam {r['O_fam']} ratio {r['ratio']} log10p {r['log10p']} rank {r['rank']}",
              f"ArmB {r['armB_pass']}", f"syn {r['syn_ratio']} contrast {r['syn_contrast']}",
              f"cc {r['completeness']}/{r['contamination']}", r["source"][:30], r["environment"][:40], sep="  ")


if __name__ == "__main__":
    main()
