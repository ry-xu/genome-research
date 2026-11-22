# Code

Run everything from the folder root with `uv run python -I <script>`.

**Runs on what is here:**
- `check_numbers.py`: recounts the README's numbers from `../data/` (the census tables, the reads, the Fisher test)
  and recomputes the isolate pairing from the two GenBank records in `../data/genomes/`.
- `newleads/codon_pair.py A.gb B.gb CODON`: the isolate pairing (`../data/isolates/isolate_pair.txt`), e.g.
  `uv run python -I src/newleads/codon_pair.py data/genomes/NC_062302.1.gb data/genomes/MT483995.1.gb AGG`.
- `agg_met/fig_agg.py OUT.png`: the figure in the README (counts written into the script, from `isolate_pair.txt`).
- `audit/fisher.py`: the Fisher tests the audit compared. `audit/pair_nt.py` and `audit/nt_identity.py` (the audit's
  nucleotide-level pairing and identity) also run on `../data/genomes/`, but need MAFFT.
- `codes_census/org_trnas.py OUT.tsv genomes.fa`: the tRNA calls in `../data/trna/` (genomes as FASTA from GenBank).

**Needs the census, which is not included** (about 22,700 genomes; three files, about 215 MB, built by the lab's
genome-fetch and profile-building steps, which are not included either):
- `codes_census/census_org.py` (with `../scan_org.py` and `../codes.py`): the census: translate each genome with its
  declared code, align each protein to the profile models, record the codon at every match column.
- `codes_census/sitescan_org.py`: the census loader and the column consensus (strict = at least 90% of the
  family-weighted residues agree).
- `newleads/codon_by_clade.py` (`../data/census/rhodo_agg.tsv`, run with `Rhodophyta AGG,ATG,AGA,CGN,ATA M,R,I`),
  `newleads/same_column.py` (`../data/census/porph_agg_columns.tsv`).
- `codes_census/reads_codons.py` (`../data/reads/`; also needs minimap2, samtools and the read runs, with
  `MIN_MAPQ=0 MAX_DIV=0.05`).
- `audit/census_check.py`, `audit/red_by_file.py`: the audit's own recount of the census numbers.

**Kept as run, by the checking agent** (paths refer to its scratch folder):
- `blind_recheck/`: the blinded re-check's own pipeline: Codetta runs (`run_codetta.sh`, `parse_codetta.py`), its
  reference-protein likelihood (`select_refs.py`, `m2.py`; `m2_v1.py`, `m2_v2.py` are earlier versions), the frozen
  calling rule (`tiers.py`) and diagnostics.
