# Data

Small tables only. The census files they come from (`mito.census.npz`, 18,015 RefSeq mitochondrial genomes;
`gbmito.census.npz`, 1,275 more GenBank species; `gbmitoorf.census.npz`, 3,408 GenBank genomes without CDS
annotation, translated from ORFs: 22,698 genomes, 14–168 MB each) are not included; `../src/README.md` says how they
were built. The column tables in `sites/` take all three as input; with fewer files the family weights change and a column or two
crosses the 90% line. Times are US Eastern, November 2026.

## `sites/`
- `pneumo_cols.tsv`: per *Pneumocystis* genome (`src/tga_arg/pneumo_cols.py Pneumocystis out.tsv <3 census files>`):
  UGA (TGA) at strict (≥ 90%) Arg and strict Trp columns, UGA at all placed sites and at Arg-majority ones, UGG at
  strict Trp and Arg, and regular Arg codons at strict Arg and Trp. The figure is drawn from it
  (`src/tga_arg/fig_pneumo.py`).
- `pneumo_same_cols.tsv`: at the 23 strict-Arg columns where *P. jirovecii* has UGA, the codon each of the other 8
  *Pneumocystis* genomes has (`src/tga_arg/pneumo_same_cols.py`).
- `isolates.tsv`: the codon at each UGA site of the reference in the reference (JX499143.1 = NC_020331.1) and the 18
  more *P. jirovecii* isolates, by the flank method (`src/tga_arg/isolates.py`). The flank method misses one site
  (JX855938 atp8:41, UGA by whole-genome alignment); the counts in the README include it.
- `pneumo_cun.tsv`: CUN, UUR and ACN codons at well-conserved columns, by consensus residue, per genome
  (`src/pneumo_cun.py`).

## `rna/`
Per-run codon tallies at every placed site of the *P. jirovecii* reference (10 runs of ENA PRJEB3063,
ERR135872–ERR135881) and of *P. carinii* NC_013660.2 (SRR2156857): gene, codon index, the genome's codon, profile
column, column consensus and share, strict flag, reads covering, reads with the genome's codon, other codons seen.
`summary.txt` sums the 10 *P. jirovecii* runs at the UGA sites (`src/tga_arg/tally_rna.py`). Reads were mapped by
`src/codes_census/reads_codons.py` with `MIN_MAPQ=0 MAX_DIV=0.05`.

## `genomes/`
- `fig_genomes.fa`: the two mitochondrial genomes in the README figure, from GenBank: *P. jirovecii* NC_020331.1 and
  *Pneumocystis* sp. 'macacae' MT726214.1 (host *Macaca mulatta*). `src/fig_uga_arg_rows.py` checks that each 27-nt
  stretch it draws occurs once in them, at the base it gives.

## `checks/`
- `blinded_panel.md`: the blinded re-check: 10 coded genomes, the pass rule fixed before launch, the agent's calls
  and the score.
- `prepost_check.md`: the pre-post check by a fresh agent with its own code (11 Nov): every census and RNA number
  re-derived, the circularity test, the isolate recount, the cox1 heme-binding pair (99.5% of 6,188 genomes).
- `lit_search.md`: the first literature search by a fresh agent (11 Nov).
- `citation_check.md`: the 20 Nov fact-check: all 10 citations opened and quoted, the searches for new work and for
  *Pneumocystis* proteomes, and the text fixes (applied).
- `tweet_check.md`: line-by-line check of the thread and its figure (21 Nov). The figure script it describes was
  not kept; `src/fig_uga_arg_rows.py` is a re-write that redraws the posted figure from `genomes/`.

These reports were written in the lab's working tree; paths in them (`W/`, `work/new-stops-euk/`, `results/tga_arg/`)
refer to it. `results/tga_arg/` there holds what `sites/` and `rna/` hold here.
