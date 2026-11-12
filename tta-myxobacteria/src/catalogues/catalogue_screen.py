"""Marker-AAI screen (src/catalogues/screen.py, FAST mode) of genomes downloaded from non-NCBI catalogues.
usage: uv run --no-project --with pyrodigal python -I src/catalogues/catalogue_screen.py SOURCE FILE...
  -> data/catalogues/<SOURCE>/faa/, results/catalogues/screen_<SOURCE>.tsv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen
ROOT = Path(__file__).resolve().parents[2]
def main():
    src, files = sys.argv[1], sys.argv[2:]
    faa = ROOT / "data/catalogues" / src / "faa"
    screen.predict(str(faa), "FAST", files, procs=3)
    out = ROOT / f"results/catalogues/screen_{src}.tsv"
    out.unlink(missing_ok=True)
    screen.run(str(ROOT / "data/catalogues/markerdb"), str(faa), str(out))
if __name__ == "__main__":
    main()
