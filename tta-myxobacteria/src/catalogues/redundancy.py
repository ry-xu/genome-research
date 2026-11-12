"""pyskani ANI of candidate genomes to the 10 clade references (redundancy check: ANI >= 99.5% = same genome/strain).
usage: uv run --no-project --with pyskani python -I src/catalogues/redundancy.py OUT.tsv genome.fna.gz...
"""
import gzip, sys
from pathlib import Path
import pyskani
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen
def seqs(p):
    return [s.upper().encode() for _, s in screen.read_fa(p)]
db = pyskani.Database()
for a in screen.CLADE:
    db.sketch(a, *seqs(ROOT / f"data/genomes/{a}.fna.gz"))
with open(sys.argv[1], "w") as o:
    o.write("path\tgenome\tbest_ref\tani\taf_query\taf_ref\n")
    for p in sys.argv[2:]:
        hs = db.query(screen.gname(p), *seqs(p))
        if hs:
            h = max(hs, key=lambda h: h.identity)
            o.write(f"{p}\t{screen.gname(p)}\t{h.reference_name}\t{h.identity*100:.2f}\t{h.query_fraction*100:.1f}\t{h.reference_fraction*100:.1f}\n")
        else:
            o.write(f"{p}\t{screen.gname(p)}\tnone\t<80\t\t\n")
# also: pairwise ANI among the candidates themselves (same-species / duplicate MAGs across catalogues)
db2 = pyskani.Database()
S = {screen.gname(p): seqs(p) for p in sys.argv[2:]}
for g, s in S.items():
    db2.sketch(g, *s)
with open(sys.argv[1].replace(".tsv", "_pairwise.tsv"), "w") as o:
    o.write("genome\tother\tani\taf_query\taf_ref\n")
    for g, s in S.items():
        for h in db2.query(g, *s):
            if h.reference_name != g:
                o.write(f"{g}\t{h.reference_name}\t{h.identity*100:.2f}\t{h.query_fraction*100:.1f}\t{h.reference_fraction*100:.1f}\n")
