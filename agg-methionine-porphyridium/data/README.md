# Data

Small tables only. The census files (about 22,700 mitochondrial genomes, about 215 MB), the profile models and the read
runs are not included. Accessions are GenBank. Times are US Eastern, November 2026. The two *P. aerugineum* isolates:
SAG 1380-2 (NC_062302.1) and CCMP1948 (MT483995.1).

## `census/`
- `rhodo_agg.tsv`: per red algal mitochondrial genome in the census (390), the codons AGG, ATG, AGA, CGN and ATA at
  strict methionine, arginine and isoleucine columns (`codon@residue`), and each codon's aligned sites (`_n`)
  (`src/newleads/codon_by_clade.py`). AGG at strict Met: 7, all in NC_062302.1; at strict Arg: 1,418.
- `porph_agg_columns.tsv`: at each of SAG 1380-2's 52 aligned AGG codons, the profile column, its consensus residue,
  whether it is strict, and the codon of each of five related red algae there (`src/newleads/same_column.py`).
- `census_check.txt`, `red_per_genome.tsv`: the audit's own recount of the census numbers (`src/audit/census_check.py`,
  `src/audit/red_by_file.py`): SAG 1380-2's codons at strict Met and Arg columns, each relative's AGG use, and the 84 of
  233 at SAG's AGG columns.

## `genomes/`
- `NC_062302.1.gb`, `MT483995.1.gb`: the two GenBank records, as fetched on 11 November 2026. `src/check_numbers.py`
  and `src/newleads/codon_pair.py` run on them.

## `isolates/`
- `isolate_pair.txt`: for each codon of SAG 1380-2 (AGG, ATG, AGA, CGN, ATA), the codon CCMP1948 has at the aligned
  position, and the reverse (`src/newleads/codon_pair.py`).
- `isolate_pair_AGG_sites.tsv`: the 58 AGG sites one by one (gene, codon, CCMP1948's codon, the two protein contexts).
- `audit_nt/fwd/`, `audit_nt/rev/`: the audit's nucleotide-level pairing in both directions (`src/audit/pair_nt.py`).
- `nt_identity.txt`: nucleotide and amino-acid identity of the 14 shared genes (87% and 88%; `src/audit/nt_identity.py`).

## `reads/`
Codon tallies from reads at every aligned AGG, ATG, AGA, CGA, CGT and ATA of SAG 1380-2 (`src/codes_census/reads_codons.py`, mapping
quality 0 allowed, reads with more than 5% mismatches dropped; columns: gene, codon index, codon, profile column,
consensus, share, strict, reads, reads matching the assembly, the two commonest other triplets).
- `paer_SRR15114507.*`: DNA (whole-genome Illumina, PRJNA744153, van Beveren et al.).
- `paer_SRR1296983.*`: RNA (MMETSP0313, poly-A transcriptome).
`*.summary.txt`: the AGG rows, then totals per codon and consensus residue (reads, reads matching, share).

## `trna/`
- `porph_trnas.tsv`: ARAGORN tRNA calls for SAG 1380-2 and the five relatives (anticodon, acceptor stem, 3:70 pair,
  N73; `src/codes_census/org_trnas.py`). `pa2_trnas.tsv`: the same for CCMP1948.

## `blind/`
- `blind_prereg.md`: the pass rule, fixed before launch. `blind_key.tsv`: the code key (not shown to the agent).
- `calls.md`: the agent's calls and method; `report.md`: the unblinding (PASS). Its scripts are in `src/blind_recheck/`.
  The coded genomes it was given are the GenBank records in `blind_key.tsv` and are not included.

## `checks/`
- `audit.md`: the code audit and sentence-by-sentence fact-check of the first draft (fresh agent, 11 Nov). The four
  errors it found were fixed before posting; the numbers in the README are the fixed ones.
- `literature.md`: the literature search of 11 Nov. The sponge precedent (*Leucosolenia*) was added on 20 Nov after
  a later prior-art check; see `../notes.md`.

These reports were written in the lab's working tree; paths in them (`results/agg_met/`, `scratchpad/`, `work/...`)
refer to it or to the checking agents' scratch folders. `results/agg_met/` there holds what `isolates/`, `reads/`,
`trna/`, `blind/` and `checks/` hold here.
