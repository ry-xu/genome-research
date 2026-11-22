# Data

Small tables only. The census files (about 22,700 mitochondrial genomes), the profile models, the genomes, the read
runs, the mass spectra and the protein search database are not included. Accessions are GenBank/RefSeq. Times are US
Eastern, November 2026. "Strict" = at least 90% of aligned mitochondrial genomes, each weighted by 1/(genomes in its
family), have one amino acid at the column.

## `census/`
- `trna_vs_usage_fungi.tsv`: all 4,638 fungal mitochondrial genomes: AGA at strict Ala and strict Arg columns, class
  (`ala`: at least 1 Ala and 0 Arg; `arg`: at least 3 Arg and 0 Ala), number of tRNA(UCU), how many have G3:U70 and
  G3:U70 with A73 (`src/aga_ala/trna_vs_usage.py`). Source of the figure, the 32 genomes, 28 of 32 and 0 of 4,322.
- `aga_ala_Fungi.tsv`: the codon contrast per genome (`src/codes_census/codon_contrast.py`): AGA at strict Ala
  (`codon@Y`) and Arg (`codon@STD`) columns, its synonyms CGN/AGG at the same, one-sided Fisher p; the last row pools
  all genomes. `aga_ala_Tremellomycetes.tsv`, `aga_ala_Exobasidiomycetes.tsv`, `aga_ala_Cystobasidiomycetes.tsv`: the
  same per class (these pooled p-values count near-identical genomes separately; see `lineage_nosel.txt`).
- `lineage_nosel.txt`: the selection-free contrast per order, each strict column counted once
  (`src/aga_ala/lineage_nosel.py`). The README table.
- `all_columns.txt`: column consensus under AGA, CGN, AGG and GCN at all aligned columns, in the 32 (`ala_class`) and in
  5 same-class genomes that read AGA as Arg (`arg_sameclass`) (`src/aga_ala/all_columns.py`; genome lists
  `ala_class.txt`, `arg_sameclass.txt`).
- `aga_usage.tsv`: AGA, AGG and CGN use in some carriers and relatives. `agg_ala.txt`: the AGG contrast per class.

## `trna/`
- `fungi_ucu.aragorn.tsv`: every ARAGORN tRNA call with anticodon UCU (TCT) in the 4,638 genomes, with acceptor stem,
  the 3:70 pair and base 73 (`src/codes_census/org_trnas.py`; the full table of all tRNAs, 15 MB, is not included;
  `trna_vs_usage.py` reads only the UCU rows).
- `basidio.aragorn.tsv`: all tRNA calls in the basidiomycete genomes. `trna_origin_basidio.tsv`: each tRNA(UCU)'s
  identity to the same genome's tRNA-Ala and other tRNA-Arg (`src/aga_ala/trna_origin.py`).

## `reads/`
Codon tallies from each assembly's source run at every aligned codon (`src/codes_census/reads_codons.py`; columns:
gene, codon index, codon, profile column, column consensus, consensus share, strict, reads, reads matching the
assembly, other triplets). *Tremella mesenterica* BK074411.1 (SRR11469946), *Tilletia horrida* BK074391.1
(SRR28015936), *Trichosporon asteroides* BK074449.1 (SRR1521503), *Salmacisia buchloeana* BK074121.1 (SRR24299923).

## `proteome/`
PXD025163 (iProX), three TMT 4-plex sets of *T. asahii*. Per set: `tally.txt` (PSMs at spectrum q ≤ 0.01 by which
versions of the mitochondrial ORFs they need; `src/prot/tally_aga.py`), `mito_psms_q01.tsv` (those PSMs),
`target_peptides_anyq.tsv` (the same at any q), `results.json` (Sage's run record), `reporters.txt` (TMT reporter
ions). Set A also has the peptide's PSMs and fragment matches (`ala_peptide_*.tsv`), the cox2 scans and
`peptide_in_6_assemblies.tsv` (`src/prot/pep_in_genome.py`). `sites.tsv`: the AGA and CGN positions in each ORF of the
database. `sage_tmt.json`, `sage_tmt_iso.json`: the Sage configurations (set A with isotope errors −1 to 3, sets B
and C with 0).

## `checks/`
- `audit/report.md`: the code audit (fresh agent), with an independent reproduction without the census
  (`star.*` are its tables; `ala_class.tsv` its list of the class).
- `blind_recheck/report.md`, `key.tsv`: the blinded re-check of 8 genomes and the key.
- `aga-ala-prepost-check.md`: the pre-post fact-check (fresh agent). `postcheck-additions-2026-11-12.md`: the check
  of the lines added after it (it also covers additions to other posts). `tweet_check.md`: line-by-line check of the
  tweet and figure (21 Nov).

These reports were written in the lab's working tree; paths in them (`results/aga_ala/`, `work/...`, `scratch/`) refer
to it or to the agents' temporary folders. `results/aga_ala/` there holds what `census/`, `trna/`, `reads/`,
`proteome/` and `checks/` hold here.
