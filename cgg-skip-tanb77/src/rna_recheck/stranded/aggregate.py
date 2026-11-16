"""Re-check, step 4: pool the per-run counts, bootstrap D over sites, apply the rule's Decision.
Follow-up S copy (recheck_S/): runs from ../followS_runs.txt, counts from recheck_S/counts/, outputs to recheck_S/;
sites still from recheck/sites_recheck.tsv. Adds D per study (per_study). --plus-old additionally pools the old
re-check's SRR37726778 units (recheck/counts/SRR37726778.*, spots from ../main_runs.txt), with the same exact-spot
check, and writes *_plus_SRR37726778.* outputs (reported only).
  .venv/bin/python -I aggregate.py --plus-old               the 8 Follow-up S runs + old SRR37726778 units
  .venv/bin/python -I aggregate.py                          the real run: counts/, runs from ../main_runs.txt
  .venv/bin/python -I aggregate.py --counts <dir> --partial  test: runs = the summaries in <dir>, no check against
                                                            main_runs.txt, outputs written to <dir>
  .venv/bin/python -I aggregate.py --totals-only            observation totals only (keeps + lacks, sites observed,
                                                            discordant; per set and per run), no keeps/lacks split,
                                                            no D; writes totals.json only. For the powered check
                                                            (Amendment 2: 10 more runs if Set A < 200 observations).
Inputs: sites_recheck.tsv (status ok, sets A and C), and for every run in ../main_runs.txt ("dir 1 SPOTS 2 RUN") its
count units in counts/: either one whole-run unit <run>.tsv, or ranges <run>.<A>-<B>.tsv (count_run.py --from/--to).
Exact-spot check (Amendment 4: spots 1..SPOTS): every unit exited 0; the ranges start at 1, are contiguous with no
gap or overlap, and end at SPOTS; each range holds exactly B - A + 1 spots, its first read name is spot A and its last
spot B; the spots summed over the units equal SPOTS. Any failure stops the aggregate.
D = lacks / (keeps + lacks), pooled over a set's sites. 95% interval: percentile bootstrap over the sites with >= 1
observation (10,000 resamples, seed 1).
Outputs: per_site_counts.tsv (every ok site of A and C, zeros included), summary.json (numbers for RECHECK.md)
"""
import csv
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

R = "results/rna_cgg/recheck/"  # sites_recheck.tsv, old SRR37726778 units
S = "results/rna_cgg/recheck_S/"  # Follow-up S counts and outputs
OLD = "SRR37726778"
STUDY = {"SRR377267": "PRJNA1441222", "SRR296696": "PRJNA1129603"}  # run-accession prefix -> study
NBOOT = 10_000


def boot(k, l, rng):
    """pooled D and its 95% percentile interval, resampling sites with >= 1 observation"""
    n = k + l
    m = n > 0
    k, l, n = k[m], l[m], n[m]
    if not len(n):
        return None, None, None
    d = l.sum() / n.sum()
    idx = rng.integers(0, len(n), size=(NBOOT, len(n)))
    ds = l[idx].sum(1) / n[idx].sum(1)
    return float(d), float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5))


def decide(A, C):
    """the rule's Decision section, on Set A against Set C"""
    powered = A["obs"] >= 200 and A["sites_observed"] >= 30 and C["obs"] >= 1000
    if not powered:
        return False, "not powered (Amendment 2: 10 more runs if Set A < 200 observations, cap 30; else reported as is)"
    if A["hi"] < 0.10:
        return True, "transcripts keep the C"
    if A["D"] >= 0.25 and A["lo"] >= max(5 * C["D"], 0.05):
        return True, "transcripts lack the C"
    return True, "inconclusive"


def main():
    args = sys.argv[1:]
    cdir = args[args.index("--counts") + 1].rstrip("/") + "/" if "--counts" in args else S + "counts/"
    odir = cdir if "--counts" in args else S
    plus = "--plus-old" in args
    suf = f"_plus_{OLD}" if plus else ""
    partial = "--partial" in args
    rng = np.random.default_rng(1)
    sites = [r for r in csv.DictReader(open(R + "sites_recheck.tsv"), delimiter="\t")
             if r["status"] == "ok" and r["set"] in ("A", "C")]
    info = {r["site_id"]: r for r in sites}
    summaries = defaultdict(list)  # run -> [(path stem, meta)]
    for p in sorted(glob.glob(cdir + "*.summary")):
        meta = dict(x.rstrip("\n").split("\t") for x in open(p))
        summaries[meta["run"]].append((p[:-len(".summary")], meta))
    if plus:
        assert OLD not in summaries, "SRR37726778 units already in the Follow-up S counts"
        for p in sorted(glob.glob(R + f"counts/{OLD}.*.summary")):
            meta = dict(x.rstrip("\n").split("\t") for x in open(p))
            assert meta["run"] == OLD, p
            summaries[OLD].append((p[:-len(".summary")], meta))
    if partial:
        runs, expect = sorted(summaries), {}
    else:
        lines = [x.split() for x in open(S + "../followS_runs.txt") if x.strip()]
        if plus:
            old = [x.split() for x in open(R + "../main_runs.txt") if x.strip() and x.split()[4] == OLD]
            assert len(old) == 1, old
            lines += old
        for x in lines:
            assert len(x) == 5 and x[1] == "1" and x[3] == "2", x
        runs, expect = [x[4] for x in lines], {x[4]: int(x[2]) for x in lines}
        assert len(set(runs)) == len(runs), "duplicate run in the run list"
        extra = sorted(set(summaries) - set(runs))
        if extra:
            print("WARNING: counts of runs not in followS_runs.txt are ignored:", extra)
    print("runs", len(runs), runs)
    tot = defaultdict(lambda: np.zeros(5, int))  # keeps lacks disc keeps_m1plus lacks_m1plus
    per_run = {}
    run_meta = {}
    for run in runs:
        parts = summaries[run]
        assert parts, (run, "no count units")
        for stem, meta in parts:
            assert meta["exit"] == "0", (stem, meta)
            assert os.path.exists(stem + ".tsv"), (stem, "no .tsv")
        ranged = [(int(m["from"]), int(m["to"]), st, m) for st, m in parts if m.get("from", "") != ""]
        if not ranged:
            assert len(parts) == 1, (run, "more than one whole-run unit")
            spots = int(parts[0][1]["spots"])
        else:
            assert len(ranged) == len(parts), (run, "whole-run and range units mixed")
            ranged.sort()
            assert ranged[0][0] == 1, (run, "ranges do not start at spot 1", ranged[0][:2])
            for (a1, b1, _, _), (a2, b2, st, _) in zip(ranged, ranged[1:]):
                assert a2 == b1 + 1, (run, "gap or overlap between ranges", (a1, b1), (a2, b2))
            for a, b, st, m in ranged:
                assert int(m["spots"]) == b - a + 1, (st, "spots", m["spots"], "expected", b - a + 1)
                assert int(m["first_spot"]) == a and int(m["last_spot"]) == b, (st, m["first_spot"], m["last_spot"])
            spots = sum(int(m["spots"]) for _, m in parts)
            assert spots == ranged[-1][1] - ranged[0][0] + 1
        assert partial or spots == expect[run], (run, "spots", spots, "expected", expect[run])
        run_meta[run] = {"spots": spots, "reads": sum(int(m["reads"]) for _, m in parts), "units": len(parts),
                         "unit_seconds": sum(int(m["seconds"]) for _, m in parts)}
        pr = defaultdict(lambda: np.zeros(5, int))
        for stem, _ in parts:
            for r in csv.DictReader(open(stem + ".tsv"), delimiter="\t"):
                v = np.array([int(r[k]) for k in ("keeps", "lacks", "discordant", "keeps_m1plus", "lacks_m1plus")])
                tot[r["site_id"]] += v
                pr[r["site_id"]] += v
        per_run[run] = pr
    if "--totals-only" in args:
        ids = {s: [r["site_id"] for r in sites if r["set"] == s] for s in ("A", "C")}

        def totals(idl, counts):
            v = [counts[i] for i in idl if i in counts]
            return {"sites": len(idl), "sites_observed": sum(int(x[0] + x[1] > 0) for x in v),
                    "obs": sum(int(x[0] + x[1]) for x in v), "discordant": sum(int(x[2]) for x in v)}
        t = {"runs": run_meta, "sets": {s: totals(ids[s], tot) for s in ("A", "C")},
             "per_run": {run: {s: totals(ids[s], per_run[run]) for s in ("A", "C")} for run in runs}}
        A, C = t["sets"]["A"], t["sets"]["C"]
        t["powered"] = A["obs"] >= 200 and A["sites_observed"] >= 30 and C["obs"] >= 1000
        json.dump(t, open(odir + f"totals{suf}.json", "w"), indent=1)
        print(json.dumps({k: t[k] for k in ("sets", "powered")}, indent=1))
        for run in runs:
            print(run, run_meta[run]["spots"], "A obs", t["per_run"][run]["A"]["obs"], "C obs",
                  t["per_run"][run]["C"]["obs"])
        return
    with open(odir + f"per_site_counts{suf}.tsv", "w") as fh:
        fh.write("site_id\tset\tgenome\tstrand\trun_len\tcodon\tkeeps\tlacks\tdiscordant\tkeeps_m1plus\tlacks_m1plus\n")
        for r in sites:
            v = tot[r["site_id"]]
            fh.write(f'{r["site_id"]}\t{r["set"]}\t{r["genome"]}\t{r["strand"]}\t{r["run_len"]}\t{r["codon"]}\t'
                     + "\t".join(map(str, v)) + "\n")

    def setstats(ids, counts):
        k = np.array([counts[i][0] if i in counts else 0 for i in ids])
        l = np.array([counts[i][1] if i in counts else 0 for i in ids])
        dc = int(sum(counts[i][2] for i in ids if i in counts))
        d, lo, hi = boot(k, l, rng)
        return {"sites": len(ids), "sites_observed": int(((k + l) > 0).sum()), "keeps": int(k.sum()),
                "lacks": int(l.sum()), "obs": int((k + l).sum()), "discordant": dc, "D": d, "lo": lo, "hi": hi}

    ids = {s: [r["site_id"] for r in sites if r["set"] == s] for s in ("A", "C")}
    out = {"runs": run_meta, "sets": {s: setstats(ids[s], tot) for s in ("A", "C")}}
    A, C = out["sets"]["A"], out["sets"]["C"]
    out["powered"], out["decision"] = decide(A, C)
    # extras (reported, not part of the decision)
    out["C_walks_with_cgg"] = setstats([i for i in ids["C"] if info[i]["walk_has_cgg"] == "1"], tot)
    out["C_walks_without_cgg"] = setstats([i for i in ids["C"] if info[i]["walk_has_cgg"] == "0"], tot)
    out["per_run"] = {run: {s: setstats(ids[s], per_run[run]) for s in ("A", "C")} for run in runs}
    study_of = {run: next((v for k, v in STUDY.items() if run.startswith(k)), "other") for run in runs}
    per_study = defaultdict(lambda: defaultdict(lambda: np.zeros(5, int)))
    for run in runs:
        for i, v in per_run[run].items():
            per_study[study_of[run]][i] += v
    out["study_of_run"] = study_of
    out["per_study"] = {st: {s: setstats(ids[s], per_study[st]) for s in ("A", "C")} for st in sorted(per_study)}
    out["by_run_len"] = {s: {L: setstats([i for i in ids[s] if min(int(info[i]["run_len"]), 3) == L], tot)
                             for L in (1, 2, 3)} for s in ("A", "C")}
    out["by_codon_C"] = {c: setstats([i for i in ids["C"] if info[i]["codon"] == c], tot) for c in ("CGA", "CGC", "CGT")}
    genomes = sorted({r["genome"] for r in sites})
    out["per_genome"] = {g: {s: setstats([i for i in ids[s] if info[i]["genome"] == g], tot) for s in ("A", "C")}
                         for g in genomes}
    out["sites_ge10_A"] = sorted([(i, int(tot[i][0]), int(tot[i][1])) for i in ids["A"] if tot[i][:2].sum() >= 10],
                                 key=lambda x: -(x[1] + x[2]))
    nC10 = [i for i in ids["C"] if tot[i][:2].sum() >= 10]
    out["sites_ge10_C"] = {"n": len(nC10), "with_any_lack": sum(int(tot[i][1] > 0) for i in nC10),
                           "D_values": sorted(float(tot[i][1] / tot[i][:2].sum()) for i in nC10 if tot[i][1] > 0)}

    # strandedness: share of observations whose read 1 is antisense to the gene (sense if dUTP-style read 2)
    def strand(counts, idl):
        n = r1_sense = 0
        for i in idl:
            if i not in counts:
                continue
            o = counts[i][0] + counts[i][1]
            plus = counts[i][3] + counts[i][4]
            n += o
            r1_sense += plus if info[i]["strand"] == "+" else o - plus
        return {"obs": int(n), "read1_sense": int(r1_sense), "share_read1_sense": (r1_sense / n) if n else None}
    out["strandedness"] = {run: {s: strand(per_run[run], ids[s]) for s in ("A", "C")} for run in runs}
    out["strandedness_all"] = {s: strand(tot, ids[s]) for s in ("A", "C")}
    json.dump(out, open(odir + f"summary{suf}.json", "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("sets", "powered", "decision")}, indent=1))
    for run in runs:
        a, c = out["per_run"][run]["A"], out["per_run"][run]["C"]
        print(run, run_meta[run]["spots"], "A", a["obs"], a["D"], "C", c["obs"], c["D"],
              "r1 sense A/C", out["strandedness"][run]["A"]["share_read1_sense"],
              out["strandedness"][run]["C"]["share_read1_sense"])


if __name__ == "__main__":
    main()
