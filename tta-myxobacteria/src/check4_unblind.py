"""Check 4 unblinding: join the blinded tblastn re-check (check4/results/codon_table.tsv, genome_summary.tsv) with
the key (results/check4_key.tsv) and the candidates' codons (results/gtdb/candidates.tsv). Frozen rule (PREREG,
check 4): a candidate (genome, X) is confirmed when X passes in that genome, i.e. O >= 5, O/E >= 3, log10p <= -4;
the rule is re-applied here from O, ratio and log10p, and must agree with the subagent's pass column.
Controls are reported with every codon that passes in them (false positives of the method).

usage: uv run --no-project python -I src/check4_unblind.py  -> results/check4_verdict.txt
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "check4" / "results"


def tsv(fn):
    return list(csv.DictReader(open(fn), delimiter="\t"))


def ok(r):
    return int(r["O"]) >= 5 and float(r["ratio"]) >= 3 and float(r["log10p"]) <= -4


def main():
    key = {r["id"]: r for r in tsv(ROOT / "results" / "check4_key.tsv")}
    cand = {r["gtdb_acc"]: r for r in tsv(ROOT / "results" / "gtdb" / "candidates.tsv")}
    tab = {}
    for r in tsv(RES / "codon_table.tsv"):
        assert ok(r) == (r["pass"] in ("1", "True", "true", "PASS", "yes")), r
        tab[(r["id"], r["codon"])] = r
    summ = {r["id"]: r for r in tsv(RES / "genome_summary.tsv")}
    lines, conf = [], 0
    lines.append("candidates (id, genome, genus, X: O n_hits E ratio log10p -> check 4; other passing codons)")
    for gid, k in sorted(key.items(), key=lambda kv: (kv[1]["role"], cand.get(kv[1]["gtdb_acc"], {}).get("codon", ""),
                                                      kv[1]["gtdb_acc"])):
        if k["role"] != "candidate":
            continue
        c = cand[k["gtdb_acc"]]
        x = c["codon"]
        r = tab.get((gid, x))
        others = [cd for (g, cd), rr in tab.items() if g == gid and cd != x and ok(rr)]
        n = summ.get(gid, {}).get("n_hits", "?")
        if r is None:
            lines.append(f"  {gid} {k['gtdb_acc']} {c['genus']} {x}: not scored (n_hits {n}) -> not confirmed")
            continue
        conf += ok(r)
        lines.append(f"  {gid} {k['gtdb_acc']} {c['genus']} {x}: O {r['O']} n_hits {r['n_hits']} E {float(r['E']):.2f} "
                     f"ratio {float(r['ratio']):.2f} log10p {float(r['log10p']):.2f} -> "
                     f"{'CONFIRMED' if ok(r) else 'not confirmed'}; others passing: {others or 'none'}")
    lines.insert(0, f"check 4: {conf} of {sum(k['role'] == 'candidate' for k in key.values())} candidates confirmed")
    lines.append("controls (id, genome, role, n_hits, canonical-stop share at e*+1; passing codons)")
    fp = 0
    for gid, k in sorted(key.items()):
        if k["role"] == "candidate":
            continue
        ps = [f"{cd} (O {rr['O']}, log10p {float(rr['log10p']):.1f})" for (g, cd), rr in sorted(tab.items())
              if g == gid and ok(rr)]
        fp += bool(ps)
        s = summ.get(gid, {})
        lines.append(f"  {gid} {k['gtdb_acc']} {k['role']} n_hits {s.get('n_hits', '?')} "
                     f"stop share {s.get('canonical_stop_share_e1', '?')}: {', '.join(ps) or 'none'}")
    lines.insert(1, f"controls with any passing codon: {fp} of {sum(k['role'] != 'candidate' for k in key.values())}")
    (ROOT / "results" / "check4_verdict.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
