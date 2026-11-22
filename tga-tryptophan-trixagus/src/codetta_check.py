"""Independent method for flagged genomes: full-genome Codetta (codon-hunt's confirm_genome.py, imported unchanged).
Six-frame translation with stops as X vs all Pfam-A_enone profiles (pyhmmer); each codon's meaning is inferred from
the profile columns it aligns to. This shares nothing with the read-through test except the genome itself: no gene
calls, no extension past a stop, no DIAMOND, no Swiss-Prot.
Pfam-A_enone (Codetta's Pfam 35 build, eddylab.org/publications/Shulgina21) lives in ~/dev/tools/codetta-pfam.
usage: uv run --no-project --with pyhmmer,numpy,pandas,requests python src/codetta_check.py <acc> [cpus]
writes results/codetta/<acc>.{json,codons.tsv}
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "codon-hunt" / "src"))
import confirm_genome as cg  # noqa: E402

cg.RES = Path.home() / "dev" / "tools" / "codetta-pfam"
cg.OUT = HERE.parent / "results" / "codetta"
cg.GEN = HERE.parent / "data" / "codetta_genomes"
cg.OUT.mkdir(parents=True, exist_ok=True)
cg.GEN.mkdir(parents=True, exist_ok=True)

if __name__ == "__main__":
    cg.run(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 3)
