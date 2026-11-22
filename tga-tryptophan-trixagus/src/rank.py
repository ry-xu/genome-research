"""Score read-through results with stop-codon-census rank_readers' formula (unchanged) and the frozen bars.
Stage 1 dirs: results/rt/ and the 662 QC-failed genomes stop-codon-census already tested (its results/rt/).
Stage 2: results/rt_full/ (every gene) replaces stage 1 for genomes retested there.
  flag: obs >= 10, excess >= 8, rate >= 0.1, abs >= 4      (rank_readers' bar)
  lead: obs >= 3, excess >= 4, rate >= 0.1 at stage 1      (-> retest every gene; README amendment 11:31)
  near: stage 2, not flagged, obs >= 10, rate >= 0.25, excess >= 4, abs >= 4  (reported, never called; README 11:52)
-> results/scores.tsv (every genome x stop), results/leads.txt (stage-2 todo), flags printed
usage: python src/rank.py
"""
import gzip
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCC = ROOT.parent / "stop-codon-census" / "results" / "rt"


def score(d):
    rt, te, rates = d["readthrough"], d["tested"], d.get("rates") or {}
    out = []
    for s in ("TAG", "TGA", "TAA"):
        o = [x for x in ("TAA", "TAG", "TGA") if x != s]
        other_rate = sum(rt[x] for x in o) / max(1, sum(te[x] for x in o))
        exp = te[s] * other_rate + 0.5
        est = {x: rt[x] / rates.get(x, 1.0) for x in rt}
        out.append(dict(stop=s, obs=rt[s], tested=te[s], rate=round(rt[s] / max(1, te[s]), 3), exp=round(exp, 2),
                        excess=round(rt[s] / exp, 1), abs=round(est[s] / max(1.0, max(est[x] for x in o)), 1),
                        other=f"{rt[o[0]]}/{te[o[0]]} {rt[o[1]]}/{te[o[1]]}"))
    return out


def main():
    tiers = {}
    with open(ROOT / "results" / "tiers.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tiers[r["acc"]] = r
    for a in (ROOT / "results" / "controls.txt").read_text().split():
        tiers.setdefault(a, {"tier": "ctrl"})
    res = {}
    for f in SCC.glob("*.stats.json"):
        a = f.name[:-11]
        if a in tiers:
            res[a] = ("scc", json.loads(f.read_text()))
    for sub in ("rt", "rt_full"):
        for f in (ROOT / "results" / sub).glob("*.stats.json"):
            res[f.name[:-11]] = (sub, json.loads(f.read_text()))
    rows = []
    for a, (src, d) in res.items():
        if d.get("readthrough") is None:
            continue
        t = tiers.get(a, {})
        for r in score(d):
            r.update(acc=a, src=src, tier=t.get("tier", ""), bp=d["bp"], contigs=d["contigs"],
                     cm2=t.get("CM2 completeness (%)", ""), cm1=t.get("CM1 completeness (%)", ""),
                     ncbi=t.get("NCBI exclude from RefSeq", ""))
            r["flag"] = int(r["obs"] >= 10 and r["excess"] >= 8 and r["rate"] >= 0.1 and r["abs"] >= 4)
            r["lead"] = int(src != "rt_full" and r["obs"] >= 3 and r["excess"] >= 4 and r["rate"] >= 0.1)
            r["near"] = int(src == "rt_full" and not r["flag"] and r["obs"] >= 10 and r["rate"] >= 0.25 and
                            r["excess"] >= 4 and r["abs"] >= 4)
            rows.append(r)
    cols = ["acc", "tier", "src", "stop", "obs", "tested", "rate", "exp", "excess", "abs", "other", "flag", "lead", "near", "bp",
            "contigs", "cm1", "cm2", "ncbi"]
    with open(ROOT / "results" / "scores.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: -r["excess"]):
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    leads = sorted({r["acc"] for r in rows if r["lead"] and r["acc"] not in
                    {x for x, (s, _) in res.items() if s == "rt_full"}})
    (ROOT / "results" / "leads.txt").write_text("".join(f"{a}\tfull\n" for a in leads))
    by_tier = {}
    for a, (src, d) in res.items():
        by_tier.setdefault(tiers.get(a, {}).get("tier", "?"), set()).add(a)
    print("genomes scored by tier:", {k: len(v) for k, v in sorted(by_tier.items())})
    print("leads awaiting stage 2:", len(leads))
    for s, kind in [(x, k) for k in ("flag", "near") for x in ("TAG", "TGA", "TAA")]:
        hits = [r for r in rows if r["stop"] == s and r[kind]]
        print(f"\n== {s}: {len(hits)} {'flagged' if kind == 'flag' else 'near-miss (not called)'}")
        for r in sorted(hits, key=lambda r: -r["excess"])[:80]:
            print(f"  {r['acc']:22s} {r['tier']:4s} {r['src']:7s} {r['obs']:4d}/{r['tested']:<5d} {r['rate']:<6} "
                  f"x{r['excess']:<6} abs {r['abs']:<5} others {r['other']:<14} {r['bp'] / 1e6:.2f} Mb "
                  f"{r['contigs']} ctg CM1 {r['cm1']} CM2 {r['cm2']}  {r['ncbi'][:60]}")


if __name__ == "__main__":
    main()
