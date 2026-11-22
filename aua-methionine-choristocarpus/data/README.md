# Data

Small tables only. The genome downloads (Genoscope assemblies, GenBank records), the RNA-seq reads (ENA ERR12954301)
and the BLAST work files are not included. Times are US Eastern, November 2026. *Choristocarpus tenellus*: Genoscope
assembly GCA_964200185.2 (mitochondrial contig CAXNAD020000071.1) and GenBank ON552980 / RefSeq NC_085326.
*Discosporangium mesarthrocarpum*: GCA_964200235.2 (contig CAXNAE020000396.1).

## `codons/`
- `sites_all.parquet`: one row per target codon and reference protein (best tblastn HSP), for the two *Choristocarpus*
  sequences, *Discosporangium* and five brown-algal controls (`src/kill/k1_audit.py`). `src/check_numbers.py` recounts
  the README's numbers from it.
- `b_read.out`: AUN codons at conserved Met and Ile positions per genome, and the call (`src/read/b_read.py`).
  `b_read_context.out`: the same for *Sphaerosorus* (positive control), *Schizocladia* and *Nematochrysopsis*.
- `a_contigs.out`: the mitochondrial contigs found in the two Genoscope assemblies (`src/read/a_contigs.py`).

## `rna/`
- `k12_rnaseq.out`: the RNA test reported in the README: 13,678,105 reads of ERR12954301 mapped to NC_085326, read
  bases at codon position 3 per site class (`src/kill/k12_rnaseq.py`). `k13_rna_vs_dna.out`: depth by region.
- `c_rna.out`, `c_rna_sites.tsv`: the first, smaller RNA run (4.5 million reads) under the pre-registered rule
  (`src/read/c_rna.py`; no call by the rule, 270 A and 0 G at Met-position AUA).

## `composition/` (Altair's alignment-free check)
- `comp.out`: per stramenopile mitochondrial genome (212), AUG, AUA and AUU+AUC per 1,000 codons of annotated genes
  (`src/composition/comp.py`).
- `gene.out`: *Choristocarpus* per gene against the brown-algal median (`gene.py`). `groups.out`: the same pooled into
  respiratory genes and ribosomal-protein genes plus *tatC* (printed by an inline run; no script saved).
- `cons.out`: the residues the relatives have where *Choristocarpus* has AUA, AUG or AUU, per gene group (`cons.py`).
  `fish.out`: AUA facing leucine against AUG (`fish.py`).

## `screen/`
- `RULES_ata.md`, `NOTE.md`, `d_test.out`, `flagged.tsv`: the screen of 20,336 organelle records declared to read AUA
  as isoleucine, its 19 flags and the recall on records declared AUA = Met (`src/screen/`).
- `priorart_k1.md`: the first prior-art check (on *Choristocarpus* and *Chloroparvula*).

## `checks/`
- `RULES.md`: the rules, written before the downloads were read, with timed amendments. `RESULT.md`: the results log.
- `kill/`: the kill agent's report (`REPORT.md`) and its outputs (`k1`-`k13`).
- `priorart.md`: the prior-art search. `dbcheck/`: the reference-database check (`DBCHECK.md`, `d*.out`) and the raw
  tRNAscan-SE and ARAGORN runs (`trna/`).
- `factcheck.md`, `factcheck2.md`, `factcheck3.md`: fact-checks of the write-up by fresh agents, in order.
- `tweet_check.md`: line-by-line check of the thread (21 Nov).

These reports were written in the lab's working tree; paths in them (`checks/chor_wgs/...`, `work/...`, `scratchpad/`)
refer to it. `checks/chor_wgs/` there holds what `codons/`, `rna/` and `checks/` hold here.
