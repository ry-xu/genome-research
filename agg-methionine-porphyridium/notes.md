# Notes: a red alga whose mitochondria read AGG as methionine

The agents' full write-up (method, evidence, checks, open points). The short version is in [`README.md`](README.md).
Times are US Eastern, November 2026. File paths are in this folder.

**Summary.** In the mitochondria of the red alga *Porphyridium aerugineum*, the codon AGG most likely means
methionine, not arginine. No red alga was known to read AGG that way. **No protein data exist for this alga, and the
tRNA that reads AGG is unknown, so the case rests on where AGG sits, in two genomes, with reads.**

## Background

Genes are read three letters (one codon) at a time. In the standard genetic code, 61 codons each name one of 20 amino
acids and 3 are stops. Some lineages, many of them mitochondria, read a few codons differently; NCBI lists 26
alternatives to the standard code. Red algal mitochondria use a code in which AGG is arginine (NCBI table 4; in the
census, 372 of 390 red algal records declare table 4, 17 table 1 and 1 table 3, all with AGG = arginine).

*P. aerugineum* is a single-celled freshwater red alga. Two isolates have had their mitochondrial genomes sequenced,
by two different labs: SAG 1380-2 (NC_062302.1; van Beveren et al. 2022) and CCMP1948 (MT483995.1; Kim et al. 2022).

## Finding

In both isolates, AGG sits where related proteins have methionine, and where the two isolates' genes differ, AGG
swaps with ATG (the methionine codon) and never with an arginine codon.

- **Where AGG sits.** At conserved methionine positions (at least 90% of mitochondrial genomes, weighted by family,
  have methionine there), SAG 1380-2 has ATG 13 times and AGG 7 times (plus AGT 1, ATA 1). At conserved arginine
  positions it has AGG 0 times; AGA, the other AG arginine codon, is there 13 times. Of its 52 aligned AGG codons, 19
  sit at columns whose most common amino acid is methionine (then leucine 12, isoleucine 4). Across 390 red algal
  mitochondrial genomes, AGG sits at conserved methionine positions 7 times (all 7 in this genome) and at conserved
  arginine positions 1,418 times. (`data/census/rhodo_agg.tsv`, `data/census/census_check.txt`)
- **In relatives.** At the positions where SAG 1380-2 has AGG, five related red algae (*P. purpureum* MT483996.1,
  *Erythrolobus coxiae* NC_062401.1, *Dixoniella grisea* NC_062400.1, *Rhodella violacea* NC_062303.1, *Neorhodella
  cyanea* NC_062397.1) have ATG 84 of 233 times and an arginine codon once (a CGA in *Rhodella* rps3, at a column that
  is not strict). (`data/census/porph_agg_columns.tsv`)
- **The two isolates** have the same codon at 66% of aligned positions (2,850 of 4,320, over 14 shared genes). Where
  SAG 1380-2 has AGG, CCMP1948 has AGG 37 times, ATG 17, something else 4, and an arginine codon 0 times (58 sites).
  Where SAG 1380-2 has ATG (71 sites), CCMP1948 has ATG 42 times and AGG 16. Where SAG 1380-2 has one of its other
  arginine codons (AGA, CGN; 78 sites), CCMP1948 has the same codon 41 times, a different arginine codon 34 times,
  something else 3 times, and ATG or AGG never. So changed AGG sites go to ATG 17 times and to an arginine codon never,
  while changed arginine sites go to another arginine codon 34 times and to ATG never (one-sided Fisher test p 7e-14,
  on the table [[17, 0], [0, 34]]). The reverse direction (CCMP1948's 57 AGG codons) gives AGG 37, ATG 16, arginine 0.
  (`data/isolates/isolate_pair.txt`, `isolate_pair_AGG_sites.tsv`)
- **Raw reads** (both runs are strain SAG 1380-2). DNA (SRR15114507, the van Beveren et al. project): each of the 7
  AGG codons at conserved methionine positions reads AGG in 99.7% or more of 1,253 to 1,837 reads (the lowest is
  1,249 of 1,253); across all AGG codons at methionine-consensus positions, 30,855 of 30,909 reads (99.8%) say AGG,
  and of the other 54 at least 4 say ATG (the table lists only the two commonest other triplets per site). RNA
  (SRR1296983, MMETSP0313, a poly-A transcriptome with few reads on these genes): 13 of 13 reads say AGG at the 5
  covered conserved sites, and 33 of 33 at methionine-consensus sites. So AGG is in both the DNA and the messenger
  RNA; it is not made by messenger-RNA editing. (`data/reads/`)
- **A blinded re-check** by a separate agent with its own code (Codetta, a code-prediction tool, plus its own
  likelihood against Swiss-Prot mitochondrial proteins), told only to infer each of 9 coded genomes' code: AGG =
  methionine in both *P. aerugineum* isolates (Codetta M at 25–26 columns; its own method LLR against arginine 72–74),
  in none of 5 related red algae, and in the known AGG = methionine green alga (*Chromochloris*) used as a positive
  control. It also recovered the *Scenedesmus* control's known code (AGG = alanine, TAG = leucine, TCA stop-like;
  NCBI table 22). The pass rule was written at 15:00 on 11 November, before launch. (`data/blind/`)
- **The tRNA** (transfer RNA, the adaptor that reads a codon) is unresolved. SAG 1380-2 has a tRNA with AGG's
  anticodon (CCU; GenBank calls it tRNA-Arg) and none with AGA's (UCU), which 4 of its 5 relatives have. CCMP1948 has
  neither. The SAG 1380-2 genome encodes only 9 tRNAs, so most must be imported from the cytoplasm anyway, and the
  AGG decoder may be among them. (`data/trna/`)

## Method and checks

- **Census.** About 22,700 mitochondrial genomes (all of RefSeq plus one per additional species from GenBank;
  18,015 + 1,275 + 3,408 in the three census files) were translated with their declared codes, every protein aligned
  to profile models, and each codon placed at its profile column (`src/codes_census/census_org.py`). A column is
  strict ("conserved") if at least 90% of genomes, weighted by 1/(genomes in their family), have one amino acid there
  (`src/codes_census/sitescan_org.py`). A scan over all codons flagged AGG in NC_062302.1 as fitting methionine's
  positions. The census was built to look for new stop codons; this came out of it. The census files (about 215 MB)
  and the profile models are not included. CCMP1948 is not in the census (one genome per species).
- **Isolate pairing** (`src/newleads/codon_pair.py`). Each shared CDS of the two GenBank records is translated with
  table 4 (in which AGG = arginine, which if anything biases the alignment against the claim), aligned globally at the
  protein level (BLOSUM62), and the codons paired through the aligned blocks.
- **Raw reads** (`src/codes_census/reads_codons.py`). Reads mapped with minimap2 to the genes; reads with more than 5%
  mismatches dropped; mapping quality 0 allowed; each read counted only if it covers all three bases of the codon.
- **Audit** (fresh agent; `data/checks/audit.md`, scripts in `src/audit/`). No bug in the pairing code. It reproduced
  the pairing with three other routes (MAFFT nucleotide alignment, a Biopython nucleotide alignment, and a protein
  alignment with AGG masked in both isolates): arginine at SAG's AGG sites is 0 in every route and both directions;
  ATG is 16–17 depending on the alignment (one rps3 site sits in an indel). It re-counted the census numbers with its
  own code and found four errors in the first draft (a wrong red algal total, a wrong first author, a Fisher p printed
  next to the wrong fractions, and an overclaim about *P. purpureum*); all were fixed before posting. It also found
  that the CCMP1948 cox1 annotation drops a codon at an exon junction where the genome most likely has AGG too, which
  would make the shared AGG count 38 of 59.
- **Literature** (`data/checks/literature.md`, plus the audit's and later prior-art checks).

## What was known

- **AGG = methionine elsewhere.** No NCBI code (tables 1–33) has it. It is reported in the mitochondria of a
  calcareous sponge, *Leucosolenia complicata*, where AGA and AGG both mean methionine (Lavrov et al. 2016,
  *Curr. Biol.*, as restated in Lavrov & Pett 2016 and Lukeš et al. 2025), and of one green alga, *Chromochloris
  zofingiensis* (Noutahi et al. 2019, *Mol. Biol. Evol.*; Žihala & Eliáš 2019, *Genome Biol. Evol.*), and in
  uncultivated Bacilli and two phage assemblies from the same samples (Shulgina & Eddy 2021, *eLife*). The 2016
  *Curr. Biol.* paper itself was not opened (its abstract does not mention the code); the sponge sentence rests on the
  same group's restatement, Lavrov & Pett 2016, *Genome Biol. Evol.* (PMC5633667): "in the calcaronean sponge L.
  complicata, where AGR codons have been reassigned to Methionine", and on Lukeš et al. 2025, *Curr. Biol.*, Table S2.
- **Red algae.** We found no report of a changed AGG in any red alga. Neither *P. aerugineum* genome paper (van
  Beveren et al. 2022, *Genome Biol. Evol.* 14:evac037; Kim et al. 2022, *BMC Biol.* 20:2) mentions a code change; Kim
  et al. annotate with table 4.

## Open points

- **Whether the proteins carry methionine.** We found no protein (mass-spectrometry) data for *P. aerugineum*; the
  two sets in PRIDE (PXD052401, PXD011275) are for *P. purpureum*, which rarely uses AGG.
- **Which tRNA reads AGG,** and so whether AGG is read only as methionine or sometimes still as arginine. SAG
  1380-2's CCU tRNA could be an arginine tRNA, so part of AGG may still be read as arginine.
- **How far it reaches.** *P. purpureum* and the other relatives sampled rarely use AGG (2, 8, 6, 11 and 9 aligned
  sites); where they do, it sits mostly at arginine positions and never at conserved methionine ones, and the blinded
  re-check did not call it methionine in any of them. Other *P. aerugineum* strains are unsequenced.
- **RNA reads are few** (13 at conserved sites), and **CCMP1948 has no public reads**, so its codons are
  assembly-only.
- **The isolates are 87% identical at the nucleotide level** (88% in amino acids; `data/isolates/nt_identity.txt`),
  so they may be two species. That strengthens the pairing but leaves the naming loose.

## What would settle it

A mass-spectrometry run on SAG 1380-2 (the strain is in the SAG culture collection): peptides with methionine or
arginine at these positions.

## References

- Kim D, Lee J, Cho CH, Kim EJ, Bhattacharya D, Yoon HS (2022). *BMC Biol.* 20:2.
- Lavrov et al. (2016). *Curr. Biol.* (not opened; see above).
- Lavrov & Pett (2016). *Genome Biol. Evol.* (PMC5633667).
- Lukeš J, Eliáš M, Kachale A, van der Gulik PTS, Speijer D (2025). *Curr. Biol.* 35:R1104–R1126, Table S2.
- Noutahi et al. (2019). *Mol. Biol. Evol.* 36:766.
- Shulgina & Eddy (2021). *eLife* 10:e71402.
- van Beveren et al. (2022). *Genome Biol. Evol.* 14:evac037.
- Žihala & Eliáš (2019). *Genome Biol. Evol.* 11:2992.
