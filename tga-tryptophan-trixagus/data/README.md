# Data

Small tables only. The read-through results for the 10,593 genomes (about 100 MB), the genomes themselves, the
GTDB R232 quality-control table (5.6 MB, `qc_failed_r232.tsv` from GTDB's `auxillary_files/`), the per-genome score
table and the raw-read DIAMOND hits of the RF2 scan (15 MB) are not included. Genomes are GenBank assemblies; the
seven throscid Shikimatogenerans genomes are from BioProject PRJNA1062330 (Wierz, Malm, Kaltenpoth and Engl, Max
Planck Institute for Chemical Ecology; paper listed as unpublished). Times are US Eastern, November 2026.

Short names: Tser *T. sericeus* GCA_040210545.2, Tder *T. dermestoides* GCA_040211365.2, Ttur *T. turgidus*
GCA_040409785.2 (the readers); Tcar *T. carinifrons* GCA_040211755.2, Tduv *T. duvalii* GCA_040210955.2, Tmey
*T. meybohmi* GCA_037094695.1, Aul/AUS03 *Aulonothroscus* sp. GCA_040212085.2; S_bost *S. bostrichidophilus*
GCA_029851105.1 (bostrichid outgroup); *Nardonella* of *Metamasius hemipterus* GCA_023651715.1.

## `screen/`
- `tier1.txt`, `tier2.txt`, `tier3.txt`: the genomes screened, by tier (`src/tiers.py`). `controls.txt`: the positive
  controls.
- `rank_final.txt`: the scored flags per stop codon, with tiers and controls (`src/rank.py`).
- `labels.tsv`: each TGA flag with its NCBI lineage and code and whether that lineage's code is known (`src/label.py`);
  `labels_taa_tag_leads.tsv`: the sub-threshold TAA/TAG leads.
- `placement.tsv`, `placement_controls.txt`: GTDB placement of the 38 unlabelled flags (`src/place.py`).
- `tga_meaning.tsv`: the Swiss-Prot residue at read-through TGA, tryptophan or glycine (`src/tga_meaning.py`).
- `codetta_controls.log`: the Codetta power check on two known readers.
- `qc_failed_throscid.tsv`: the GTDB R232 QC rows of the seven throscid symbionts (marker percentages, CheckM).
- `qc_bias.tsv`, `trna_unread.tsv`: two side tables of the screen, not used in the README: how often genomes of
  genera with a known code fail GTDB QC (`src/qc_bias.py`), and codons with no tRNA reader in the 69 genomes NCBI
  excluded for "missing tRNA genes" (`src/trna_missing.py`).

## `symbionts/`
- `symbiont_profile.tsv`: per genome (the throscid and bostrichid symbionts, *Nardonella*, controls): genes under
  table 4, coding fraction under tables 11 and 4, TGA and TGG per 1,000 codons, gene ends, RF1/RF2 profile hits
  (`src/symbiont_profile.py`). Its RF2 column is wrong for the four non-readers (cutoff too strict; see `notes.md`).
- `ortholog_tga/A_vs_B.tsv`: at each in-frame TGA (and TGG) codon of A, the aligned codon in B's reciprocal-best
  ortholog (`src/ortholog_tga.py`).
- `read_support/`: raw-read support at in-frame TGA and TGG sites, per genome and run (`*.txt` summary, `*.sites.tsv`
  per site; `src/read_support_tga.py`).
- `codetta/`: the blinded Codetta re-check, per genome (`*.json`: the call; `*.codons.tsv`: per codon), plus the two
  power controls (`src/codetta_check.py`).

## `throscidae/`
- `proteomes/*.faa`: table-4 proteins of the eight symbionts (pyrodigal), as used for the tree and gene content.
- `concat.faa` and `concat.faa.*`: the 99-protein concatenated alignment and the IQ-TREE LG+F+G4 run
  (`src/throscidae_tree.py`).
- `au/`: the AU test of reader monophyly under LG (`CRITERION.txt` written before it ran). `c20/`: the stopped LG+C20
  search (log and starting tree). `c20au/`: the AU test under LG+C20.
- `gene_content.tsv`: genes in all readers and no non-reader, or the reverse (`src/gene_content.py`): only RF2.
- `tga_sites.tsv`: every in-frame TGA of the eight, classed by where it sits relative to the Swiss-Prot alignment;
  `trna_trp.txt` and `trna/`: the tryptophan tRNAs (tRNAscan-SE) (`src/throscid_trp.py`).

## `rf2_readscan/`
The second fact-check's scan of each run's first 4 million read pairs for RF2 (`src/rf2_readscan/`).
`rfdb.faa` is the protein set searched (all eight proteomes plus *E. coli* RF2); `summary.txt` is
`src/rf2_readscan/rf2_scan_summary.py` run on the per-run hit tables (not included).

## `checks/`
- `factcheck1.md`, `factcheck2.md`, `factcheck3.md`: three rounds of sentence-by-sentence fact-checks by fresh agents
  (11-12 Nov), in order.
- `recheck_mintaka.md`, `recheck_mintaka_tga.out`: an independent recount of codon usage (20 Nov; script
  `src/recheck_mintaka_tga.py`, on prodigal `-g 4` gene calls not included).
- `recheck_sabik.md`: an independent recount of the GenBank pseudogenes and release-factor annotations (20 Nov).
- `priorart_altair.md`: the prior-art check (13 Nov).
- `tweet_check.md`: line-by-line check of the roundup tweet (21 Nov).

These reports were written in the lab's working tree; paths in them (`P`, `results/...`, `work/...`) refer to it.
