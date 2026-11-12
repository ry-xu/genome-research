"""Arm A injection calibration under the final null and the frozen flag rule (PREREG.md: "The final null (all GTDB)
changes these slightly; the calibration is recomputed with it and reported, and the rule is not changed").

Null as in src/armA_run.py (every scanned GTDB genome with calling table 11 and >= 300 complete genes); each
injection run is scored under the calling table of its own (injected) profile, runs with < 300 complete genes are
skipped, and a run counts as recalled when the injected codon X is flagged by armA_run.flag (T1 or T2).
Reported: recall by (s, q); null runs flagging any codon; injected runs also flagging a codon other than X.

usage: uv run --no-project --with numpy python -I src/inject_eval.py results/inject2 [results/inject ...]
"""
import collections
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import armA
from armA_run import MIN_COMPLETE, OUT, flag


def final_null():
    acc, dom, pos, stops, meta = armA.load_gtdb(sorted(OUT.glob("prof_*.npz")))
    tabs = np.array([armA.gtdb_table(s) for s in stops])
    keep = (tabs == 11) & (meta[:, 1] >= MIN_COMPLETE)
    print(f"final null: {keep.sum()} of {len(acc)} scanned genomes")
    return armA.null_from(meta[keep, 2], pos[keep])


def evaluate(d, rho, irate):
    rec, tot = collections.Counter(), collections.Counter()
    nnull = ninj = 0
    null_fl, inj_fl = [], []
    for f in sorted(Path(d).glob("*.npz")):
        z = np.load(f)
        a = f.name[:-4]
        for i, spec in enumerate(z["specs"]):
            X, s, q, frag, seed = spec.split(":")
            gc, nc = z["meta"][i][:2]
            if nc < MIN_COMPLETE:
                continue
            sc = armA.score(z[f"pos|{i}"], gc, armA.gtdb_table(z[f"stops|{i}"]), rho, irate)
            other = [r["codon"] for r in sc if r["codon"] != X and flag(r)]
            if X == "null":
                nnull += 1
                if other:
                    null_fl.append((a, frag, other))
            else:
                ninj += 1
                k = (float(s), float(q))
                tot[k] += 1
                rec[k] += any(r["codon"] == X and flag(r) for r in sc)
                if other:
                    inj_fl.append((a, X, other))
    print(f"\n{d}: {nnull} null runs, {ninj} injected runs (>= {MIN_COMPLETE} complete genes)")
    print(f"null runs flagging any codon: {len(null_fl)}/{nnull}", null_fl)
    print(f"injected runs also flagging a non-target codon: {len(inj_fl)}/{ninj}")
    for qq in sorted({k[1] for k in tot}, reverse=True):
        ks = sorted(k for k in tot if k[1] == qq)
        print(f"q = {qq}: recall " + " / ".join(f"{rec[k] / tot[k]:.2f}" for k in ks) + " at s = " +
              " / ".join(f"{k[0]}" for k in ks) + "  (n " + " / ".join(str(tot[k]) for k in ks) + ")")
    s_all = sorted({k[0] for k in tot})
    print("pooled over q = 1: " + ", ".join(f"s={s}: {rec[(s, 1.0)]}/{tot[(s, 1.0)]}" for s in s_all))


if __name__ == "__main__":
    rho, irate = final_null()
    for d in sys.argv[1:]:
        evaluate(d, rho, irate)
