"""Re-check, step 5: print the numbers for RECHECK.md section 6 from summary.json (run after aggregate.py).
  .venv/bin/python -I report.py
"""
import json

R = "results/rna_cgg/recheck/"
d = json.load(open(R + "summary.json"))


def fmt(x):
    if x["D"] is None:
        return f'obs {x["obs"]:,} (keeps {x["keeps"]:,}, lacks {x["lacks"]:,}) at {x["sites_observed"]:,} sites; D n/a'
    return (f'obs {x["obs"]:,} (keeps {x["keeps"]:,}, lacks {x["lacks"]:,}) at {x["sites_observed"]:,} of '
            f'{x["sites"]:,} sites, discordant {x["discordant"]:,}; D {x["D"]:.4f} [{x["lo"]:.4f}, {x["hi"]:.4f}]')


print("decision:", d["decision"], "| powered:", d["powered"])
for s in ("A", "C"):
    print(s, fmt(d["sets"][s]))
print("C walks with CGG   ", fmt(d["C_walks_with_cgg"]))
print("C walks without CGG", fmt(d["C_walks_without_cgg"]))
print("\nper run (spots; A; C; read-1-sense share A, C):")
for run, v in d["per_run"].items():
    st = d["strandedness"][run]
    sh = [st[s]["share_read1_sense"] for s in ("A", "C")]
    print(f'  {run} {d["runs"][run]["spots"]:,}')
    print(f'    A {fmt(v["A"])}')
    print(f'    C {fmt(v["C"])}')
    print(f'    read-1 sense share A {sh[0] if sh[0] is None else round(sh[0], 4)} '
          f'(n {st["A"]["obs"]:,}), C {sh[1] if sh[1] is None else round(sh[1], 4)} (n {st["C"]["obs"]:,})')
sa = d["strandedness_all"]
print("strandedness all:", {s: (sa[s]["share_read1_sense"], sa[s]["obs"]) for s in ("A", "C")})
print("\nby run length:")
for s in ("A", "C"):
    for L, v in d["by_run_len"][s].items():
        print(f"  {s} len {L}{'+' if L == '3' else ''}: {fmt(v)}")
print("by codon C:")
for c, v in d["by_codon_C"].items():
    print(f"  {c}: {fmt(v)}")
print("\nSet A sites with >= 10 observations (site, keeps, lacks):", len(d["sites_ge10_A"]))
for x in d["sites_ge10_A"]:
    print("  ", *x)
print("Set C sites with >= 10 observations:", d["sites_ge10_C"]["n"], "with any lack:", d["sites_ge10_C"]["with_any_lack"])
print("\nper genome (A obs, D_A | C obs, D_C):")
for g, v in d["per_genome"].items():
    a, c = v["A"], v["C"]
    if a["obs"] or c["obs"]:
        print(f'  {g} A {a["obs"]} {a["D"] if a["D"] is None else round(a["D"], 4)} | '
              f'C {c["obs"]} {c["D"] if c["D"] is None else round(c["D"], 4)}')
