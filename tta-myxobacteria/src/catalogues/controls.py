"""Positive/negative controls for the marker-AAI screen (same FAST gene-calling mode as candidates).
 positive: the 10 clade references, (a) leave-self-out, (b) leave-own-genus-out (simulates a new genus of the clade)
 negative: (c) 84 outgroup panel genomes leave-self-out (= leave-own-genus-out; one rep per non-clade Polyangiaceae genus),
           (d) held-out reps of other Myxococcota families (results/catalogues/neg_controls_heldout.tsv), not in the DB
usage: uv run --no-project --with pyrodigal python -I src/catalogues/controls.py
"""
import os, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen
screen.diamond.__defaults__ = (1,)  # 1 thread: the NCBI driver is running alongside
ROOT = Path(__file__).resolve().parents[2]
C = ROOT / "data/catalogues/controls"
R = ROOT / "results/catalogues"
def main():
    neg = [l.split("\t")[0] for l in open(R / "neg_controls_heldout.tsv") if l.strip()]
    (C / "neg_fna").mkdir(parents=True, exist_ok=True)
    (C / "accs.txt").write_text("\n".join(neg) + "\n")
    subprocess.run(["uv", "run", "--no-project", "python", "-I", str(ROOT / "src/catalogues/fetch_ncbi.py"), str(C / "accs.txt"), str(C / "neg_fna"), "30"])
    screen.predict(str(C / "neg_faa"), "FAST", [str(p) for p in (C / "neg_fna").glob("*.fna.gz")], procs=1)
    refs = [ROOT / f"data/genomes/{a}.fna.gz" for a in screen.CLADE]
    screen.predict(str(C / "ref_fast_faa"), "FAST", [str(p) for p in refs], procs=1)
    for f in [R / "controls_neg_heldout.tsv", R / "controls_pos_self.tsv", R / "controls_pos_genus.tsv"]:
        f.unlink(missing_ok=True)
    screen.run(str(ROOT / "data/catalogues/markerdb"), str(C / "neg_faa"), str(R / "controls_neg_heldout.tsv"))
    screen.run(str(ROOT / "data/catalogues/markerdb"), str(C / "ref_fast_faa"), str(R / "controls_pos_self.tsv"))
    for g in sorted(set(screen.CLADE.values())):
        d = C / f"logo_{g}"
        d.mkdir(exist_ok=True)
        mem = [a for a, x in screen.CLADE.items() if x == g]
        for a in mem:
            l = d / f"{a}.faa"
            if not l.exists():
                os.symlink(C / "ref_fast_faa" / f"{a}.faa", l)
        screen.run(str(ROOT / "data/catalogues/markerdb"), str(d), str(R / "controls_pos_genus.tsv"), ",".join(mem))
    print("done")


if __name__ == "__main__":
    main()
