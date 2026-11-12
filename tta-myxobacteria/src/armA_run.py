"""Arm A final scoring of every scanned GTDB genome with the frozen rule (PREREG.md).

Null: every genome with calling table 11 and >= 300 complete genes (score.null_table adds the >= 1,000 sense codons
per region filter). Scored: every genome with >= 300 complete genes, under its own calling table.

usage: uv run --no-project --with numpy python -I src/armA_run.py
  -> results/gtdb/armA_flags.tsv    one row per flagged (genome, codon), with taxonomy
     results/gtdb/armA_genomes.tsv  one row per genome: table, suspect table-4 call, GC, complete genes, flags
     results/gtdb/armA_rows.tsv.gz  every scored (genome, codon) with lr >= 0.5, and TGG always (clade contrasts)
     results/gtdb/armA_missing.txt  R232 representatives with no profile
Refuses to run unless scan.log has its "done:" line and the profiles hold exactly the R232 representatives
(ALLOW_PARTIAL=1 overrides; Changes 2026-11-10 in PREREG.md).
"""
import gzip
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import armA
import codonpos as cp

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "gtdb"
META = [ROOT / "data" / "gtdb_meta" / f for f in ("bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz")]
MIN_COMPLETE = 300
TGA, TGG = cp.CIDX["TGA"], cp.CIDX["TGG"]


def flag(r):
    """Frozen rule: log10p <= -10, excess >= 10, and T1 (lr >= 1.5) or T2 (lr >= 1.0, depl <= 0.2)."""
    t1 = r["lr"] >= 1.5
    t2 = r["lr"] >= 1.0 and r["depl"] <= 0.2  # depl is nan for TGA under table 4: T1 only
    return r["log10p"] <= -10 and r["excess"] >= 10 and (t1 or t2)


def taxonomy():
    """GTDB taxonomy of every genome, and the set of species representatives (gtdb_representative == t)."""
    tax, reps = {}, set()
    for fn in META:
        with gzip.open(fn, "rt") as f:
            assert next(f).split("\t")[18:20] == ["gtdb_representative", "gtdb_taxonomy"], fn
            for line in f:
                x = line.split("\t", 20)
                tax[x[0]] = x[19]
                if x[18] == "t":
                    reps.add(x[0])
    return tax, reps


def fmt(v):
    return f"{v:.6g}" if isinstance(v, float) else str(v)


def main():
    files = sorted(OUT.glob("prof_*.npz"))
    acc, dom, pos, stops, meta = armA.load_gtdb(files)
    ncomp, gc = meta[:, 1], meta[:, 2]
    assert len(set(acc.tolist())) == len(acc), "duplicate accessions"
    assert (ncomp == stops.sum(axis=1)).all(), "meta column 1 is not the number of complete genes"
    tax, reps = taxonomy()
    done = any(line.startswith("done:") for line in open(OUT / "scan.log"))
    missing, extra = sorted(reps - set(acc.tolist())), set(acc.tolist()) - reps
    (OUT / "armA_missing.txt").write_text("".join(a + "\n" for a in missing))
    print(f"scan.log done line: {done}; {len(reps)} representatives, {len(acc)} profiles, {len(missing)} missing, "
          f"{len(extra)} not representatives")
    if not (done and not missing and not extra) and os.environ.get("ALLOW_PARTIAL") != "1":
        sys.exit("scan incomplete (ALLOW_PARTIAL=1 to score anyway)")
    tabs = np.array([armA.gtdb_table(s) for s in stops])
    keep = (tabs == 11) & (ncomp >= MIN_COMPLETE)
    rho, irate = armA.null_from(gc[keep], pos[keep])
    nflag = 0
    with open(OUT / "armA_flags.tsv", "w") as ff, open(OUT / "armA_genomes.tsv", "w") as fg, \
            gzip.open(OUT / "armA_rows.tsv.gz", "wt") as fr:
        cols = ["codon", "aa", "n3", "nI", "exp3", "excess", "lr", "log10p", "depl"]
        ff.write("acc\tdomain\ttable\tsuspect_t4\tgc\tncomplete\t" + "\t".join(cols) + "\trule\ttaxonomy\n")
        fr.write("acc\ttable\t" + "\t".join(cols) + "\tflag\n")
        fg.write("acc\tdomain\ttable\tsuspect_t4\tgc\tncomplete\tin_null\tflags\ttaxonomy\n")
        for i in range(len(acc)):
            a, t = acc[i], int(tabs[i])
            trp = pos[i, TGA].sum() + pos[i, TGG].sum()
            suspect = int(t == 4 and pos[i, TGA].sum() < 0.5 * trp)
            fl = []
            if ncomp[i] >= MIN_COMPLETE:
                for r in armA.score(pos[i], gc[i], t, rho, irate):
                    f = flag(r)
                    if r["lr"] >= 0.5 or r["codon"] == "TGG":
                        fr.write(f"{a}\t{t}\t" + "\t".join(fmt(r[c]) for c in cols) + f"\t{int(f)}\n")
                    if f:
                        fl.append(r["codon"])
                        rule = "T1" if r["lr"] >= 1.5 else "T2"
                        ff.write(f"{a}\t{dom[i]}\t{t}\t{suspect}\t{gc[i]:.4f}\t{int(ncomp[i])}\t" +
                                 "\t".join(fmt(r[c]) for c in cols) +
                                 f"\t{rule}\t{tax.get(a, '')}\n")
            nflag += bool(fl)
            fg.write(f"{a}\t{dom[i]}\t{t}\t{suspect}\t{gc[i]:.4f}\t{int(ncomp[i])}\t{int(keep[i])}\t"
                     f"{','.join(fl)}\t{tax.get(a, '')}\n")
    print(f"{len(acc)} genomes, {int((ncomp >= MIN_COMPLETE).sum())} scored, null {int(keep.sum())}, "
          f"{nflag} with a flag")


if __name__ == "__main__":
    main()
