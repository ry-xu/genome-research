"""Follow-up S re-check: print the numbers for RECHECK_S.md from summary.json and summary_plus_SRR37726778.json.
  ../recheck/.venv/bin/python -I report_S.py
Extra, not in the rule: Clopper-Pearson one-sided 97.5% upper bound when a set has 0 lacks (1 - 0.025**(1/n)).
"""
import json
import os

S = "results/rna_cgg/recheck_S/"


def fmt(x):
    if x["D"] is None:
        return f'obs {x["obs"]:,} (keeps {x["keeps"]:,}, lacks {x["lacks"]:,}) at {x["sites_observed"]:,} sites; D n/a'
    return (f'obs {x["obs"]:,} (keeps {x["keeps"]:,}, lacks {x["lacks"]:,}) at {x["sites_observed"]:,} of '
            f'{x["sites"]:,} sites, discordant {x["discordant"]:,}; D {x["D"]:.4f} [{x["lo"]:.4f}, {x["hi"]:.4f}]')


def cp0(n):
    return 1 - 0.025 ** (1 / n) if n else None


def show(path, label):
    d = json.load(open(path))
    print(f"===== {label}: {os.path.basename(path)}")
    print("runs:", list(d["runs"]))
    print("spots:", sum(v["spots"] for v in d["runs"].values()), "reads:", sum(v["reads"] for v in d["runs"].values()))
    print("decision:", d["decision"], "| powered:", d["powered"])
    for s in ("A", "C"):
        x = d["sets"][s]
        print(s, fmt(x))
        if x["lacks"] == 0:
            print(f'   CP 97.5% upper, per obs {cp0(x["obs"]):.5f}; per observed site {cp0(x["sites_observed"]):.4f}')
    print("C walks with CGG   ", fmt(d["C_walks_with_cgg"]))
    print("C walks without CGG", fmt(d["C_walks_without_cgg"]))
    print("\nper run: spots, reads, units, unit-seconds; A; C; read-1 sense share (A, C) and m1plus counts")
    for run, v in d["per_run"].items():
        st, rm = d["strandedness"][run], d["runs"][run]
        print(f'  {run} [{d["study_of_run"][run]}] spots {rm["spots"]:,} reads {rm["reads"]:,} units {rm["units"]} '
              f'sec {rm["unit_seconds"]}')
        print(f'    A {fmt(v["A"])}')
        print(f'    C {fmt(v["C"])}')
        for s in ("A", "C"):
            sh = st[s]["share_read1_sense"]
            print(f'    read-1 sense {s}: {st[s]["read1_sense"]:,} of {st[s]["obs"]:,} = '
                  f'{"n/a" if sh is None else round(sh, 4)}')
    print("\nper study:")
    for stdy, v in d["per_study"].items():
        print(f'  {stdy}: runs {[r for r, x in d["study_of_run"].items() if x == stdy]}')
        print(f'    A {fmt(v["A"])}')
        print(f'    C {fmt(v["C"])}')
    sa = d["strandedness_all"]
    print("strandedness all:", {s: (sa[s]["read1_sense"], sa[s]["obs"], sa[s]["share_read1_sense"]) for s in "AC"})
    print("\nby run length:")
    for s in ("A", "C"):
        for L, v in d["by_run_len"][s].items():
            print(f"  {s} len {L}{'+' if L == '3' else ''}: {fmt(v)}")
    print("by codon C:")
    for c, v in d["by_codon_C"].items():
        print(f"  {c}: {fmt(v)}")
    print("\nSet A sites with >= 10 observations (site, keeps, lacks):", len(d["sites_ge10_A"]))
    for x in d["sites_ge10_A"][:25]:
        print("  ", *x)
    print("Set C sites with >= 10 observations:", d["sites_ge10_C"]["n"], "with any lack:",
          d["sites_ge10_C"]["with_any_lack"], "D values:", d["sites_ge10_C"]["D_values"][:20])
    print("\nper genome, sorted by Set A obs (A obs, sites, lacks, D_A | C obs, lacks, D_C):")
    g = sorted(d["per_genome"].items(), key=lambda kv: -kv[1]["A"]["obs"])
    totA = d["sets"]["A"]["obs"]
    for name, v in g:
        a, c = v["A"], v["C"]
        if a["obs"] or c["obs"]:
            print(f'  {name} A {a["obs"]} ({a["sites_observed"]} sites, {100 * a["obs"] / totA if totA else 0:.1f}%) '
                  f'lacks {a["lacks"]} D {a["D"] if a["D"] is None else round(a["D"], 4)} | '
                  f'C {c["obs"]} lacks {c["lacks"]} D {c["D"] if c["D"] is None else round(c["D"], 4)}')
    print("genomes with Set A obs:", sum(v["A"]["obs"] > 0 for _, v in g), "with Set C obs:",
          sum(v["C"]["obs"] > 0 for _, v in g))
    print()


show(S + "summary.json", "Follow-up S, 8 new runs (the decision)")
if os.path.exists(S + "summary_plus_SRR37726778.json"):
    show(S + "summary_plus_SRR37726778.json", "Follow-up S + old SRR37726778 (reported only)")
