# Code

Run everything from the folder root with `uv run python -I <script>`.

**Runs on what is here:**
- `check_numbers.py`: recounts the README's numbers (166 / 1 / 4 aligned stops, glutamine at 44 of 91 clean TAG
  sites, the gene ends, TAG at 41 of 60 RNA-covered sites, 12.3 / 14.4 / 0.8 per 1,000) from `../data/`.
- `review/analyse_sites.py data/review/bil`: reproduces `../data/review/bil.analysis.txt`.
- `review/terminal_codons.py data/ends/GCA_*.ends.tsv`: reproduces `../data/ends/terminal_codons_relatives.tsv`.
- `rna/site_table.py data/review/bil.sites.tsv data/rna/ERR*.sites.tsv`: reproduces `../data/rna/all3.table.txt`.
- `rna/strain_compare.py data/rna/{LLXJ,ALZF,JIWJ}.codons.tsv.gz`: reproduces `../data/rna/strain_compare.txt`.

**Needs the genomes, which are not included.** Download GCA_965178125.1 and GCA_008330645.1 from NCBI and put the
FASTA files in `data/dl/` (`GCA_965178125.1_puBilCfla1.1_genomic.fna`, `GCA_008330645.1_CrBVI_genomic.fna`):
- `fig_bilabrum_enoa.py OUT.png`: the figure in the README.
- `review/stops_in_hsps.py GENOME.fna HITS.tsv OUTPREFIX`: the per-site tables in `../data/review/` from the tblastn
  hits (`tblastn -db_gencode 1` of `human76.fa`, outfmt 6 with the extra columns listed in the script).
- `review/codons.py`: per-position codons of one query in one genome.
- `rna/run_1kp.sh`: streams the three 1KP runs from the SRA AWS mirror (fasterq-dump) into minimap2, then
  `rna/site_reads.py` (codons at the stop sites); `rna/codon_pairs.py GENOME HITS OUTPREFIX ALN.paf ...`: the
  genome-against-read codon tables (`*.codons.tsv.gz`, `*.pairs.txt`). Also needs fasterq-dump and minimap2. Paths in
  `run_1kp.sh` refer to the lab's working tree.

**The screen, kept as run** (`screen/`; reads the stored DIAMOND scans of the lab's genome-scan project, which are
not included; paths refer to the lab's working tree): `units.py` (the unit list), `stopsense.py` (re-align each
unit's hit bodies and tally codon against reference residue), `score.py` and `call.py` (the frozen call rules),
`lineage.py` (each unit's NCBI lineage) and `label.py` (known / lineage / new labels against
`../data/screen/known_codes.tsv`), `stage2.py` and
`stage2_table.py` (contamination and dispersion tests against a control codon), `summary.py`
(`../data/screen/summary_all.txt`).

Not kept: the script behind the eRF1 count over UniProt reference proteomes (see `../data/README.md`).
