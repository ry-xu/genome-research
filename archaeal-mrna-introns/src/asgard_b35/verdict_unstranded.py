"""The words of culture/PREREG_CULTURE.md for runs that strand_identity.py finds not stranded.

  cd culture && uv run --project ../../reading-ledger python src/verdict_unstranded.py RNA_OUTNAME DNA_OUTNAME

Written 13 Nov 23:45 after the count, because verdict_culture.py (registered) implements only the stranded case, and the
registered text says: "If a run is not stranded, all its templates count and I say so." Same thresholds as
verdict_culture.py; J and T are templates of both orientations; the one-removal share is read from
out/<RNA>_mol_gaps.tsv (templates of both orientations with the commonest removal at the site, slid left).
"""
import csv, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def table(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


rna, dna = sys.argv[1], sys.argv[2]
M = {(r["site"], int(r["place"]), r["class"]): r for r in table(os.path.join(HERE, "out", f"{rna}_mol_sites.tsv"))}
G = table(os.path.join(HERE, "out", f"{rna}_mol_gaps.tsv"))
D = {r["site"]: r for r in table(os.path.join(HERE, "out", f"{dna}_sites.tsv"))}
for p in table(os.path.join(HERE, "data", "predicted.tsv")):
    site = p["site"]
    j, t = M.get((site, 0, "joined"), {}), M.get((site, 0, "through"), {})
    J, T = int(j.get("templates") or 0), int(t.get("templates") or 0)
    inf = J + T
    gs, gl = j.get("sense_gap_start"), j.get("sense_gap_len")
    k = sum(int(g["templates"]) for g in G if g["gap_start"] == gs and g["gap_len"] == gl) if J else 0
    one_gap = J > 0 and 5 * k >= 4 * J
    frame = j.get("frame_restored") == "True"
    if inf < 10:
        word = "unread"
    elif J >= 5 and 2 * J >= inf and one_gap and frame:
        word = "removed in RNA"
    elif inf >= 20 and J == 0:
        word = "not removed"
    else:
        word = "mixed"
    dj, dt = int(D[site]["joined"]), int(D[site]["through"])
    if word == "removed in RNA" and dj > 0:
        word = "withdrawn (the stretch is missing from some DNA)"
    ends = (gs == p["pred_start"] and gl == p["pred_len"]) if J else ""
    print(f"{site} {p['fam']} ({p['label']}): all templates joined {J}, through {T}; commonest removal {gs}+{gl} in {k};"
          f" predicted {p['pred_start']}+{p['pred_len']}; ends as predicted: {ends}; DNA joined {dj}, through {dt}")
    print(f"   WORD BY THE RULE (run not stranded, all templates): {word}")
