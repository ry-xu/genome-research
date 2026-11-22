# Data

Small tables only. The genomes, the 1KP read runs and the read alignments (PAF, about 425 MB) are not included.
Genomes: *Bilabrum* cf. *latius* GCA_965178125.1 (Wellcome Sanger Institute, Aquatic Symbiosis Genomics; chromosomes
OZ237735-OZ237777), *Cafeteria roenbergensis* GCA_008330645.1, *Caecitellus* GCA_965644315.1, *Halocafeteria*
GCA_057816405.1, uncultured Cafeteriaceae GCA_964026335.1. Times are US Eastern, November 2026.

## `review/`
The independent re-check (a separate agent): tblastn of 76 conserved human proteins (`human76.fa`, accessions in
`human76_acc.txt`) against each genome with the standard code (`bil_human76.tsv`, `caf_human76.tsv`; outfmt 6 with
qlen, sframe, slen, qseq, sseq). From these, `src/review/stops_in_hsps.py` wrote, per genome (`bil.*`, `caf.*`):
- `*.sites.tsv`: every distinct in-frame stop site, with the human residue opposite it, the distance from the
  alignment ends, the query gaps nearby and the identity in a 25-column window (the "clean" filter).
- `*.ends.tsv`, `*.hsps.tsv`, `*.pairs.tsv`: gene ends after full-length hits, per-locus best hits (GC3, stop counts),
  and codon-residue pairs.
- `*.analysis.txt` (`src/review/analyse_sites.py`) and `*.summary.txt`: the counts in the README (166 aligned TAG,
  1 TAA, 4 TGA; 91 clean TAG, glutamine at 44; *Cafeteria* TAG 0).

## `ends/`
Gene ends from the screen's stored DIAMOND scans (first version; `GCA_*.ends.tsv`: in-frame body and downstream sequence of every
best hit that reaches its reference's C-terminus). `terminal_codons_relatives.tsv` (`src/review/terminal_codons.py`,
duplicate loci removed): *Bilabrum* TAG 0 of 149, *Cafeteria* 11 of 145, the other three 35-39 of 149-172.

## `rna/`
Reads of three 1KP cultures that carry a related *Bilabrum* (ENA ERR2041128 = LLXJ *Chroodactylon*, ERR2041020 = ALZF
*Halochlorococcum*, ERR3487419 = JIWJ *Acrosiphonia*), mapped to the genome with minimap2 (`src/rna/run_1kp.sh`;
`ERR*.mm2.log`).
- `ERR*.sites.tsv`, `ERR*.summary.txt`: read codons at each stop site of `review/bil.sites.tsv`
  (`src/rna/site_reads.py`). `ERR*.table.txt` per run, `all3.table.txt` pooled (`src/rna/site_table.py`): clean TAG
  sites covered 60, majority TAG 41, CAG 18.
- `*.codons.tsv.gz`, `*.pairs.txt`: genome codon against read codon at every codon of the 76-protein panel
  (`src/rna/codon_pairs.py`); `ALL3.pairs.txt` is the pooled table behind 12.3 (CAG to TAG), 14.4 (CTG to TTG) and
  0.8 per 1,000 (His to Tyr, Arg to Trp and Pro to Ser together).
- `strain_compare.txt` (`src/rna/strain_compare.py`): the three cultures carry one strain (97.9-99.0% of shared
  codons agree), so the runs are pooled, not counted as replicates.

## `trna/`
tRNAscan-SE 2.0.13 on the genome (`trna.out`, `trna.ss`, `trna.stats`), and the isotype models (cmscan,
TRNAinf-euk-iso) on the CUA tRNA, its pseudogene copy and the genome's tRNA-Gln (`trna_check.fa`, `iso.tbl`,
`iso.out`).

## `erf1/`
`erf1_aln.fa`: *Bilabrum*'s eRF1 (partial, exons 2-4) aligned with human, *Cafeteria*, two other stramenopiles and
*Tetrahymena*; domain 1 reads TA**V**NIKS where human has TASNIKS (S60). The count over UniProt reference proteomes
(2,714 of 2,776 with serine, none valine) was made by the screening agent; its script and reference table were not
kept. The pre-post check re-derived it by two other methods (2,708 of 2,752 and 2,792 of 2,875 serine; valine 0 in
all), see `checks/prepost_check.md`.

## `ssu/`
`ssu_queries.fa`: the NCBI 18S sequences used to confirm the organism's identity (chromosome 12 carries *Bilabrum
latius* 18S at 99.6-99.8%; the host dinoflagellate's 18S is on one unplaced contig).

## `screen/`
The eukaryote-wide screen that found it (1,927 units scored, 1,785 non-duplicate). `controls.tsv`: the 67 units of
known code the call rules were set on; `known_codes.tsv`: the ledger of known nuclear stop reassignments;
`summary_all.txt`: coverage, calls and recall of the full run; `GCA_965178125.1_TAG.txt`: the stage-2 tests for
*Bilabrum* TAG against the control codon CAA.

## `checks/`
- `evidence.md`: the screening agent's evidence log, in order (11 Nov).
- `factcheck.md`: the first fact-check and literature search by a fresh agent (11 Nov).
- `prepost_check.md`: the pre-post check by a third fresh agent, which re-derived every number (12 Nov).
- `tweet_check.md`: the check of the roundup tweet (21 Nov).

These reports were written in the lab's working tree; paths in them (`results/bilabrum/`, `review/`, `rna_map/`,
`scratchpad/`) refer to it. `results/bilabrum/review/` there is `review/` here, and `results/bilabrum/rna_map/` is
`rna/`.
