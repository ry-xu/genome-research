# Code

Run from the folder root. Usage is in each script's docstring.

**Runs directly on the included `data/` tables:**
- `uv run python -I src/check_numbers.py`: recounts and asserts every headline number in `README.md` from `data/`.
- `uv run python -I src/fig_l02_public.py`: rebuilds `figures/asgard_b35_l02.png` from `data/asgard_b35/`.

**Full analysis pipelines (require public genomes from GTDB r232 / NCBI and SRA runs not included here):**
- `census/`: 3,666-genome archaeal DNA scan, BHB motif caller, decoy generator, supergroup comparison, and RefSeq vs.
  MAG stratification (`work/archaeal-intron-census/src/`).
- `asgard_b35/`: B35 culture read mapping, molecule collapsing, junction parsing, PCR-skip test (`c1.py`), and
  DNA-hairpin scan (`c2.py`, `c2_read.py`).
- `tack_c0/`: TACK metatranscriptome screening, metagenome DNA checks, 22-nt RNA hairpin background scan, and
  25-genome BHB fold check (`c0_homologs2.py`).

Scripts copied from the lab under `census/`, `asgard_b35/`, and `tack_c0/` retain their original relative/lab paths as
a faithful record of how the tables in `data/` were produced. External tools used: minimap2, SRA toolkit, DIAMOND,
ViennaRNA, primer3, Sage.
