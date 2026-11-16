# Code

Run everything from the folder root with `uv run python -I <script>`. Large inputs (genomes, Swiss-Prot, reads, UHGV)
are not included: `../scripts/setup_reference_data.sh` fetches the reference data, and `paths.py` says where each input
is expected (each can be moved with an environment variable). Scripts write into working folders that are not
included (`results/`, `scratch/`, `recheck_work/`, `phages_work/`); the tables they made that support the result are
in `../data/` (see `../data/README.md`).

- **Frameshift scan:** `fs_scan.py` (DIAMOND blastx `-F 15`, one genome set per call), `summarize.py` (break classes
  and the per-genome flag rules), `run_set.py` (a genome set in batches), `select_set.py` and `meta.py` (one genome
  per GTDB family), `prok_db.py` (the Swiss-Prot database), `candidate.py`, `pool.py` (pooled test),
  `inject.py` and `score_inj.py` (planted shifts and recall). `fetch_set.py`, `fetch_controls.py`, `trna_products.py`
  and `codon_atlas/` fetch genomes and make the ARAGORN tRNA calls.
- **The CGG test and its follow-ups** (`analysis/`): `cgg_test.py` (calls per genome), `shift_residue.py` (residue
  opposite the skip), `inframe_cgg.py` (CGG read in frame), `skip_context.py` (codon before the skip; N after it),
  `cgg_reader_census.py` (genomes with no CGG tRNA), `cgg_figure.py`, `allele_kmers.py` (raw reads),
  `site_turnover.py`, `turnover_prfb_windows.py`, `turnover_ggn_pool.py` (turnover), `reader_contig.py` (the two
  reader genomes).
- **RNA test:** `rna_cgg.py` (sites, known-answer test, counting, decision), `rna_stream.py` (streams SRA runs into
  minimap2), `rna_strand.py` (strandedness probe), `rna_skiplike.py` (Amendment 5 split). `proteomics_db.py` is
  included only for the gene walk that `rna_cgg.py` imports. `rna_recheck/`: the independent aligner-free recount
  (`first/` for the first 10 runs, `stranded/` for the stranded follow-up).
- **CUU context:** `skip_context_chicken.py`; `context_recheck/` is the independent re-derivation.
- **Blinded re-check:** `recheck/` (miniprot; numbered in run order). It expects its working folder in `recheck_work/`
  (genomes G1-G9, references, outputs).
- **tRNA check in the pig genomes:** `trna_pig/` (tRNAscan-SE and ARAGORN, the planted 8-nt-loop control).
- **Phages:** `phages/main/` (gene calls and codon statistics per phage, groups, scores, post hoc counts),
  `phages/audit/` (the independent audit), `phages/out_of_sample/` (the referee check), `phages/pfam_junction/` (the
  Pfam check). They expect UHGV in `data/uhgv/` and their working folder in `phages_work/`.
- **Figures:** `figures/site.py` (the site's letters from the genome FASTA, into `site.json`), `figures/image.py`
  (`../figures/rpoC_skip_site.png`), `figures/anim.py` (frames of `../figures/ribosome.gif`; the frames were joined into
  the GIF separately).
