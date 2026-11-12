"""Run the injection grid on the reference panel (data/panel.tsv) -> results/inject/<acc>.npz.
Per genome: 2 null runs (as is; fragmented to 75%) and N injections with X uniform over sense codons (not ATG),
s (share of stops rewritten to X), q (share of interior X removed), frag drawn from fixed lists. Deterministic seeds.
usage: [INJECT_OUT=<dir>] python src/calibrate.py <n_inject> <workers>   (default results/inject; the rerun with the
       fixed inject.recode, PREREG Changes 19:59, is results/inject2)"""
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(os.environ.get("INJECT_OUT", ROOT / "results" / "inject"))
sys.path.insert(0, str(ROOT / "src"))
import codonpos as cp
from fetch import fetch

S = [0.005, 0.01, 0.02, 0.03, 0.05, 0.1]
Q = [1.0, 0.9, 0.5, 0.0]
FRAG = [0, 0.75]
XS = [cp.CODONS[c] for c in cp.SENSE if cp.CODONS[c] != "ATG"]


def specs(acc, n):
    rng = random.Random(acc)
    out = ["null:0:0:0:1", "null:0:0:0.75:2"]
    for i in range(n):
        out.append(f"{rng.choice(XS)}:{rng.choice(S)}:{rng.choice(Q)}:{rng.choice(FRAG)}:{100 + i}")
    return out


def run(acc, n):
    o = OUT / f"{acc}.npz"
    if o.exists():
        return
    fn = fetch(acc)
    subprocess.run([sys.executable, str(ROOT / "src" / "inject.py"), str(fn), str(o)] + specs(acc, n), check=True,
                   stdout=open(OUT / f"{acc}.log", "w"), stderr=subprocess.STDOUT)


if __name__ == "__main__":
    n, w = int(sys.argv[1]), int(sys.argv[2])
    OUT.mkdir(parents=True, exist_ok=True)
    accs = [l.split("\t")[0] for l in open(ROOT / "data" / "panel.tsv") if l.strip()]
    with ThreadPoolExecutor(w) as ex:
        for f in [ex.submit(run, a, n) for a in accs]:
            f.result()
    print("done")
