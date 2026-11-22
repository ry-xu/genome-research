# Kill check: "the Choristocarpus tenellus mito reads AUA as Met"
Written by the kill subagent, finished at about 10:15 on 14 Nov. Lyra saved it here because the subagent could not write files, and condensed it lightly. Every number comes from a printed .out file in this directory.

## Verdict
**The claim survives every check that could be run. Nothing is killed.** Two pieces of supporting wording need correcting:
- **Independence.** NC_085326 and CAXNAD020000071 are two separate sequencing and assembly efforts of one culture, strain KU-1152 / RCC7273. They are not two isolates.
- **tRNA mechanism.** A missing trnI(CAU) does not separate Choristocarpus from Discosporangium. Discosporangium also lacks the syntenic trnI(CAU), yet it reads ATA as Ile.

## 1. Code bugs: none found
- **Independent reimplementation** (k1_audit.py, using Biopython reverse complement and table 4): 0 translation mismatches over 67,301-74,019 codons per target, and its counts match b_read.out exactly.
- **Pipeline without tblastn** (k2_annot.py: annotated CDS, ATA→X, global alignment):
  - NC_085326: ATA at Met/Ile 83/5, ATG 25/3.
  - Controls: ATA 1/110, 0/117 and 4/108.
- **Hand-sliced codons**, including the minus-strand sites, all match (k9). For example, cten - 3060 is ATA at Desmarestia rpl16 ALAAR[M]GKGKG, and Saccharina has ATG in that column.

## 2. Column bias: none
- **Share of Met columns.** cten 138/5839 (0.0236), against 0.0213-0.0237 in the controls.
- **Position in the protein.** ATA@M sites by protein-position quintile: 20 22 14 20 12. Only 3 of 88 fall in the first 30 residues.
- **Spread across genes.** 20 of 21 genes have ATA@M sites. The largest single-gene share is 0.148, and the top two genes hold 0.295.
- **They are conserved Met columns.** At the 88 cten ATA@M columns, other genomes have ATG at: Saccharina 88/88, Petalonia 86, Sargassum 85, Discosporangium 82, Phaeostrophion 69/86, Protohalopteris 61/76.
- **Reverse view.**
  - At Discosporangium's 113 ATG@M columns, cten has ATA 82 and ATG 27.
  - At Discosporangium's 50 matched ATA@I columns, cten has ATT 34, TTA 7, ATC 7 and ATA 0.

## 3. Reverse view against non-brown stramenopiles (k4)
blastp of Choristocarpus proteins (ATA→X) against 2432 proteins from 62 genera that played no part in the calls (Oomycota, diatoms, other ochrophytes, Bigyra). Share of Met and Ile residues aligned at each codon:

| codon | sites | Met share | Ile share |
|---|---|---|---|
| ATA | 217 | 0.390 | 0.134 |
| ATG | 50 | 0.437 | 0.093 |
| ATT | 559 | 0.020 | 0.505 |
| ATC | 35 | 0.042 | 0.439 |

Each subgroup gives the same picture: ATA Met share 0.363-0.427.

## 4. Met→Ile substitution instead of reassignment: rejected (k9)
- If ATA were still Ile, the ATA share at Met columns should follow Choristocarpus's own Ile usage.
  - At Ile columns, ATA/ATT/ATC = 6/245/13, so ATA's share is 0.023.
  - Met columns carry 88/3/2. About 2.1 ATA would be expected; binomial P = 1.1e-137.
- Ile-type codons fill 0.674 of Choristocarpus Met columns, against 0.066 in Discosporangium, 0.052 in Phaeostrophion and 0.021 in Saccharina.
- ATG+ATA fill 115 of 138 Choristocarpus Met columns. That is close to the controls (Discosporangium 116/137, Saccharina 138/142).

## 5. Shared artefact or reference-guided assembly
- **Alignment.** All 32,351 positions align (the WGS contig is rotated: NC position 7912 is WGS position 1). There are 27 differences:
  - 24 deletions in the WGS, 1 insertion and 2 adjacent substitutions;
  - 24 of the 27 sit in homopolymer runs of 4 or more, mostly poly-T runs of 7-17, which fits different sequencing chemistries;
  - 1 falls in a consensus Met/Ile codon.
- **Strain.** The Genoscope BioSample SAMEA115803497 names KU-1152 / RCC7273 (Villefranche, 1980). It notes the culture was "originally thought to be Discosporangium mesarthrocarpum strain KU-2346 … probably identical to KU-1152". NC_085326 is also KU-1152 from Villefranche (01 Sep 1980). So these are two sequencings of one strain.
- **The Discosporangium comparison genome is genuine** (k10):
  - its mito contig matches GenBank Discosporangium cox1 at 99.812% (AB899176) and 100% (LM995414);
  - its plastid matches rbcL AB252654 at 99.857%;
  - Choristocarpus matches these at 84.4%, 85.3% and 87.6%.

## 6. tRNAs (k6-k8)
- **Standard brown algae** have three CAU tRNAs: Met-e, Ile-CAU and Met-f. Phaeostrophion has an extra Met-f.
- **Choristocarpus** has two, the same in both sequences:
  - Met-f at 2625 (identity 0.831);
  - a tRNA at 2859 closest to Met-e (0.712; Ile-CAU 0.562; best blastn hit Cutleria trnW, 78.8%/66).
  - No Ile-CAU-class tRNA was found, and the syntenic slot between trnW and trnQ is empty.
  - Recall check: the same search finds Ile-CAU in that slot in Phaeostrophion (94.2%) and Protohalopteris (96.8%).
- **Discosporangium** (reads ATA as Ile) also has an empty slot. Its CAU tRNAs:
  - Met-f;
  - Met-e;
  - a third CAT-anticodon locus, X (33026-33095), closest to trnW 0.669, Ile-CAU 0.612, Met-e 0.629. Choristocarpus 2859 matches X at 0.646.
- **Reading:** the canonical trnI(CAU) was lost before the two genera split. Discosporangium keeps one more divergent CAU tRNA than Choristocarpus does. That fits the claim, but identity scores cannot tell Ile charging from Met charging. Choristocarpus has no extra trnM copies.
- **Starts.** All 35 annotated starts are ATG, so they are uninformative (k11).

## 7. RNA editing (k12, k13)
- **Data:** Choristocarpus KU-1152 RNA-seq ERR12954301, the first 600 MB of each read file (13,678,105 reads), mapped to NC_085326 with minimap2.
- **Coverage:** 1414 reads with MAPQ ≥ 20; 21,728 of 32,351 positions at depth ≥ 3.
- **Third codon base in reads:**

  | Site class | Sites with reads | Third-base reads |
  |---|---|---|
  | ATA@M | 86 of 88 | A 796, T 1, G 0 (0 of 74 sites at depth ≥ 3 show any G) |
  | ATG@M | 26 | G 287, T 3 |
  | ATT@I | 239 | T 2031, C 1 |
  | ATA@I | 4 | A 42 |

  The five deepest ATA@M sites (depth 22-27) read 100% A.
- **RNA or DNA?** Mean depth is CDS 7.43, intergenic 5.18, rRNA 3.34, and per-gene depth ranges 13-fold (rps13 1.40 to nad7 18.16). That fits RNA, but a DNA background cannot be ruled out. Treating the deep sites as mostly RNA is the subagent's inference.

## Could not check
- Proteomics.
- tRNA charging, lysidine, and structure (no tRNAscan or aragorn).
- Whether Genoscope polished its assembly against NC_085326.
- A second isolate.
- Other Discosporangiales or close outgroups.
- Whether the library is stranded or poly-A.
- The full RNA run.

## Files
- k1_audit … k13_rna_vs_dna (.py/.out).
- sites_all.parquet, annot_sites.parquet, consensus.parquet.
- Downloads (about 1.1 GB): scratchpad/chor_kill/.
