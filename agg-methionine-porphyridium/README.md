# A red alga whose mitochondria read AGG as methionine

A gene is read three bases (one codon) at a time, and each codon names one amino acid or a stop. Many mitochondria
read a few codons differently. Red algal mitochondria use a code in which AGG is arginine.

In the mitochondria of the red alga ***Porphyridium aerugineum***, a single-celled freshwater alga, AGG most likely
means **methionine**. Two isolates have had their mitochondrial genomes sequenced, by two different labs: SAG 1380-2
(NC_062302.1) and CCMP1948 (MT483995.1). In both, AGG sits where related proteins have methionine, and where the two
isolates' genes differ, AGG swaps with ATG (the methionine codon) and **never with an arginine codon**. No red alga
was known to read AGG that way.

![Bars for three codon classes in isolate SAG 1380-2, split by the codon the other isolate, CCMP1948, has at the same aligned position. At SAG's AGG: 37 same, 17 ATG, 4 other, 0 arginine. At its ATG: 42 same, 16 AGG, 13 other. At its other arginine codons: 41 same, 34 another arginine codon, 3 other.](figures/agg_isolate_pairs.png)

*Codon at each aligned position of 14 shared mitochondrial genes, SAG 1380-2 against CCMP1948. AGG trades places
with ATG; the other arginine codons trade places with each other. Drawn by `src/agg_met/fig_agg.py`.*

---

- **AI-generated.** The census, the analysis and the write-up were done by Claude (Anthropic) agents, run and curated
  by Ryan Xu.
- **Checked by AI, not yet by a human expert.** A separate agent audited the code and every number (and reproduced
  the isolate pairing with three alignment methods); another, given 9 coded genomes and told only to infer each one's
  code, called AGG = methionine in both isolates and in none of 5 related red algae, under a pass rule fixed before it
  started. Further agents checked the citations.
- **Computational only.** No protein data for *P. aerugineum* were found, and which tRNA reads AGG is unknown.
- Corrections and criticism are very welcome: ryan@ryanxu.net.

### What was found

| | count |
|---|---|
| SAG 1380-2, conserved methionine positions | ATG 13, AGG 7 |
| SAG 1380-2, conserved arginine positions | AGG 0, AGA 13 |
| 390 red algal mitochondrial genomes, AGG at conserved methionine positions | 7 (all 7 in this genome) |
| the same, AGG at conserved arginine positions | 1,418 |
| five related red algae, at the positions where SAG 1380-2 has AGG | ATG 84 of 233, an arginine codon once |
| CCMP1948, where SAG 1380-2 has AGG (58 sites) | AGG 37, ATG 17, something else 4, arginine codon 0 |
| CCMP1948, where SAG 1380-2 has AGA or CGN (78 sites) | a different arginine codon 34, ATG or AGG never |

"Conserved" = at least 90% of about 22,700 mitochondrial genomes, weighted by family, have one amino acid there. The
two isolates have the same codon at 66% of aligned positions. Changed AGG sites go to ATG 17 times and to an arginine
codon never, while changed arginine sites go to another arginine codon 34 times and to ATG never (one-sided Fisher
test p 7e-14). Recount: [`src/check_numbers.py`](src/check_numbers.py).

### Evidence

- **Positions.** AGG at 7 conserved methionine positions and 0 conserved arginine positions in SAG 1380-2, while
  AGA, the other AG arginine codon, stays at arginine positions (13). Across 390 red algal mitochondrial genomes, all
  7 AGG codons at conserved methionine positions are in this genome.
- **Two isolates, two labs.** Where the isolates' codons differ, AGG swaps with ATG and never with an arginine codon,
  while the arginine codons swap among themselves (figure above). Two genomes from two labs agreeing, together with
  the raw reads for SAG 1380-2, make an assembly error unlikely.
- **Raw reads** (the sequencing data behind SAG 1380-2). DNA: each of the 7 AGG codons at conserved methionine
  positions reads AGG in 99.7% or more of 1,253 to 1,837 reads (the lowest is 1,249 of 1,253); across all AGG codons
  at methionine-consensus positions, 30,855 of 30,909 reads (99.8%) say AGG, and of the other 54 at least 4 say ATG.
  RNA (a poly-A transcriptome with few reads on these genes): 13 of 13 reads say AGG at the 5 covered conserved
  sites. So AGG is in both the DNA and the messenger RNA; it is not made by messenger-RNA editing.
- **Blinded re-check.** A separate agent with its own code (Codetta, a code-prediction tool, plus its own likelihood
  against Swiss-Prot mitochondrial proteins), told only to infer each of 9 coded genomes' code, called AGG =
  methionine in both *P. aerugineum* isolates, in none of 5 related red algae, and in the known AGG = methionine green
  alga (*Chromochloris*) used as a positive control.
- **The tRNA** is unresolved. SAG 1380-2 has a tRNA with AGG's anticodon (CCU) and none with AGA's (UCU), which its
  relatives have. CCMP1948 has neither. The decoder may be imported from the cytoplasm.

### What was known

- **AGG = methionine elsewhere.** No NCBI code (tables 1–33) has it. It is reported in the mitochondria of a
  calcareous sponge, *Leucosolenia complicata*, where AGA and AGG both mean methionine (Lavrov et al. 2016,
  *Curr. Biol.*, as restated in Lavrov & Pett 2016 and Lukeš et al. 2025), and of one green alga, *Chromochloris
  zofingiensis* (Noutahi et al. 2019, *Mol. Biol. Evol.*; Žihala & Eliáš 2019, *Genome Biol. Evol.*), and in
  uncultivated Bacilli and two phage assemblies from the same samples (Shulgina & Eddy 2021, *eLife*).
- **Red algae.** We found no report of a changed AGG in any red alga. Neither *P. aerugineum* genome paper (van
  Beveren et al. 2022, *Genome Biol. Evol.*; Kim et al. 2022, *BMC Biol.*) mentions a code change.

### Open points

- **Whether the proteins carry methionine.** We found no protein (mass-spectrometry) data for *P. aerugineum*; the
  two sets in PRIDE are for *P. purpureum*, which rarely uses AGG.
- **Which tRNA reads AGG** (see above), and so whether AGG is read only as methionine or sometimes still as arginine.
  SAG 1380-2's CCU tRNA could be an arginine tRNA.
- **How far it reaches.** *P. purpureum* and the other relatives sampled rarely use AGG; where they do, it sits mostly
  at arginine positions and never at conserved methionine ones, and the blinded re-check did not call it methionine
  in any of them. Other *P. aerugineum* strains are unsequenced.
- **Few RNA reads, and no reads for CCMP1948.** The RNA check rests on 13 reads at conserved sites, and CCMP1948 has
  no public reads, so its codons are assembly-only. The two isolates are 87% identical at the nucleotide level, so
  they may be two species.

### What would settle it

Peptides. A mitochondrial proteome of *P. aerugineum* SAG 1380-2 (the strain is in the SAG culture collection) would
show methionine or arginine at these positions.

### Files

- `notes.md`: the agents' full write-up (method, evidence, checks, references)
- `figures/`: the figure above
- `data/`: census tables, the isolate pairing, read tallies, tRNA calls, the blinded re-check and the check reports;
  see [`data/README.md`](data/README.md)
- `src/`: the code; see [`src/README.md`](src/README.md). `uv run python -I src/check_numbers.py` recounts the numbers
  above from `data/`, and `uv run python -I src/agg_met/fig_agg.py out.png` draws the figure. The census itself needs
  files not included here (about 215 MB).
