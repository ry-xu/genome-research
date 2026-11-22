# Code

Run from the folder root. Usage is in each file's first lines. The scripts were copied from the lab; the only edits
are path edits (absolute lab paths became repo paths, `data/external/...` placeholders or bare tool names).

**Runs on what is here** (`uv run python -I <script> OUT.png`):
- `fig_div_ptc_rows.py`: the first README figure, from `data/silent_operon/origin/aln23.sto`.
- `fig_div_rrna_bars.py`: the DNA-vs-RNA bar figure. It recounts the DNA run and the first three RNA runs from
  `data/expression_23s/` with `frozen_core/fig_div23.py`; the two later runs (0 : 74 and 0 : 43) are constants, taken
  from `data/silent_operon/rna/`, because their PAFs are too large to include.
- `frozen_core/fig_div23.py`: the technical figure (`figures/ptc_bases_and_reads.png`).

**Needs inputs that are not included** (genomes from GTDB release 232 / NCBI, SRA runs, Rfam models, PDB files):
- `frozen_core/`: the multi-copy 23S scan (`div23_scan.py`, `div23_class.py`, `div23_check.py`, `div23_chim*.py`,
  `div23_post.py`), the PTC base reads (`aq_ptc_coords.py`, `ptc_bases_pair.py`), pairing (`aq_pairing.py`), read
  pileups (`contig_pileup.py`) and expression counts (`expr_count.py`).
- `checks/kill_23s/`, `checks/expr_audit/`: the kill agent's and the expression audit's own code.
- `orthogonal_ribosome/`: the 16S tests (T1 tail channel, T4 census, T5 16S counts, T6 contacts; `t6_contacts.py`
  needs the PDB 6OG7 file, not included).
- `uba6777_reads/`: the UBA6777 read tests.
- The code of the 22 Nov checks (`data/silent_operon/`) is not included yet; their results, rules and tables are.

External tools: minimap2, SRA toolkit, HMMER, Infernal, mafft, BLAST+, wgsim, SPAdes.

## Code for the 22 Nov checks
`silent_operon/A_processing`, `B_origin` and `C_rna` hold the scripts behind the three checks summarised in the README
(processing stems and promoter; origin, neighbourhood and selection; extra RNA). They were written for the lab's local
layout, so paths start with `lab/` or `~/`. Their inputs are public genomes and SRA runs listed in each check's RESULT
and RULES files under `data/silent_operon/`.
