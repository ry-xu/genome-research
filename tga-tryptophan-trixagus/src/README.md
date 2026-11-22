# Code

Run from the folder root with `uv run python -I <script>`.

**Runs on what is here:**
- `check_numbers.py`: recounts the README's numbers (codon usage, orthologs, reads, Codetta calls, TGA sites, GTDB
  markers) from `../data/`.
- `throscid_figure.py [treefile] [out.png]`: the figure in the README, from `../data/throscidae/concat.faa.treefile`
  and `../data/symbionts/symbiont_profile.tsv`.
- `rf2_readscan/rf2_scan_summary.py HITS OWN ../data/rf2_readscan/rfdb.faa`: summarises a run's read hits (the hit
  tables are not included; `rf2_readscan.sh` and the other shell scripts made them with fastq-dump and DIAMOND, and
  expect a working folder in `$FACTCHECK_DIR`).

**Kept as run** (they expect the lab's working tree: `results/`, the streamed genomes, and the sibling projects
`stop-codon-census` and `codon-hunt`, whose modules are copied to `lab_deps/` for reading):
- The screen: `tiers.py` (tiers), `synth.py` (synthetic controls), `screen.py` (read-through test, with
  `lab_deps/readthrough.py`, `batch.py`, `batch_pooled.py`), `rank.py` (flags), `label.py` (NCBI lineage and known
  codes), `place.py` (GTDB placement), `tga_meaning.py`, `hit_taxa.py`, `qc_bias.py`, `trna_missing.py`,
  `trna_summary.py`.
- The symbionts: `symbiont_profile.py` (codon usage, coding fraction, release factors), `ortholog_tga.py`,
  `read_support_tga.py` (with `lab_deps/read_support.py`; needs minimap2 and fastq-dump), `codetta_check.py` (with
  `lab_deps/confirm_genome.py` and `codetta_lite.py`; needs Pfam-A profiles), `gene_content.py`, `throscid_trp.py`
  (needs tRNAscan-SE), `throscidae_tree.py` (needs DIAMOND, MAFFT and IQ-TREE 3).
- `rf2_readscan/`: the second fact-check agent's scripts for the RF2 read scan and the *T. dermestoides* junction
  (`junction_*`, `depth_window.py`, `holes.py`, `clip_tails.py`, `edge_reads.py`).
- `recheck_mintaka_tga.py`: the later independent recount of codon usage, run on prodigal `-g 4` gene calls.
