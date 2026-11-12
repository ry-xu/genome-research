"""Per-genome marker-AAI profile of candidates: AAI to each of the 10 clade genomes and to the 5 nearest outgroups,
plus the same profile for the 10 clade references (leave-self-out) for comparison. Same diamond filters/markers as screen.py.
Purpose: is a candidate distinctly closer to one clade genus (inside the clade) or about equally far from all (could be sister)?
usage: uv run --no-project python -I src/catalogues/candidate_profile.py OUT.tsv cand.faa...
"""
import sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen

ROOT = Path(__file__).resolve().parents[2]
DB = ROOT / "data/catalogues/markerdb/markers"
files = sys.argv[2:] + [str(ROOT / f"data/catalogues/ref_faa/{a}.faa") for a in screen.CLADE]
q = ROOT / "data/catalogues/markerdb/query_profile.faa"
with open(q, "w") as o:
    for f in files:
        o.write(open(f).read())
qbest, phit = defaultdict(dict), defaultdict(dict)
for qq, s, pid, bs in screen.diamond(q, DB, 500, threads=1):
    sg, _, m = s.split("|")
    qg = qq.split("|")[0]
    if sg == qg:
        continue
    cur = qbest[qg].get(m)
    if cur is None or bs > cur[0]:
        qbest[qg][m] = (bs, qq)
    d = phit[qq]
    if sg not in d or bs > d[sg][0]:
        d[sg] = (bs, pid)
q.unlink()
with open(sys.argv[1], "w") as o:
    o.write("genome\tis_reference\t" + "\t".join(f"{a}({g})" for a, g in screen.CLADE.items()) + "\tnearest_outgroups\n")
    for f in files:
        qg = screen.gname(f)
        aai = defaultdict(list)
        for m, (bs, qp) in qbest[qg].items():
            for sg, (b2, pid) in phit[qp].items():
                aai[sg].append(pid)
        st = {g: sum(v) / len(v) for g, v in aai.items() if len(v) >= 20}
        og = sorted(((a, g) for g, a in st.items() if g not in screen.CLADE), reverse=True)[:5]
        o.write(f"{qg}\t{int(qg in screen.CLADE)}\t" + "\t".join(f"{st[a]:.1f}" if a in st else "NA" for a in screen.CLADE)
                + "\t" + "; ".join(f"{g} {a:.1f}" for a, g in og) + "\n")
