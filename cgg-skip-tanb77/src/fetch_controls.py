"""Fetch the control genomes into the shared codon atlas cache (no copies in this project).
usage: uv run python -I src/fetch_controls.py"""
from codon_atlas.cache import fetch

CONTROLS = {
    "ecoli": "GCF_000005845.2",     # prfB +1 CUU-UGA-C, dnaX -1
    "bsub": "GCF_000009045.1",      # prfB +1
    "mtb": "GCF_000195955.2",       # Actinomycetota: no prfB shift
    "mlep": "GCF_000195855.1",      # pseudogene decay (~1,100 pseudogenes)
    "sodalis": "GCF_000010085.1",   # pseudogene-rich endosymbiont
    "efoc": "GCA_001880345.2",      # Euplotes focardii, code 10, natural +1 positive
}
for k, a in CONTROLS.items():
    print(k, a, fetch(a), flush=True)
