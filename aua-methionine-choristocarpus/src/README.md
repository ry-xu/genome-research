# Code

Run everything from the folder root with `uv run python -I <script>`.

**Runs on what is here:**
- `check_numbers.py`: recounts the README's numbers (88 of 94, AUU 245/3, the strict 68/0, 20 of 21 genes, the 67%
  and 2-7%, the controls, AUG per 1,000 codons, the 70% excess, the RNA and screen counts) from `../data/`.
- `fig_aua_bars.py OUT.png`: the figure in the README, from `../data/codons/b_read.out` and `b_read_context.out`.

**Needs files not included** (kept as run; paths refer to the lab's working tree, `lab/work/organelle-reading/`):
- `read/`: `a_contigs.py` (find the mitochondrial contigs in the Genoscope assemblies; BLAST+), `b_read.py` (the
  reading: tblastn of the 8 reference proteomes against a genome, then the codon at each conserved Met/Ile position;
  needs the lab's `data/cds_mito.parquet`), `c_rna.py` (the pre-registered RNA test; reads SAM from minimap2).
- `kill/`: the kill agent's scripts. `k1_audit.py` re-implements the reading and wrote `sites_all.parquet`;
  `k2_annot.py` is the pipeline without tblastn; `k3`-`k5`, `k9`-`k11` test column bias, the reverse view, the two
  genome sequences and the starts; `k6`-`k8` the tRNAs; `k12_rnaseq.py`, `k13_rna_vs_dna.py` the RNA.
- `dbcheck/`: the reference-database check (`fetch*.sh` downloads, `d*_*.py` read them; `d4_run_trna.sh`,
  `d5_run_trnaB.sh` run tRNAscan-SE and ARAGORN).
- `composition/`: Altair's alignment-free codon-usage check (`comp.py`, `gene.py`, `cons.py`, `fish.py`).
- `screen/`: the screen over GenBank and RefSeq organelle records (`a_anchors.py`, `b_codons.py`, `c_survey.py`,
  `d_test.py`); needs the organelle-reading ledger runs and its `read_k1` module.
