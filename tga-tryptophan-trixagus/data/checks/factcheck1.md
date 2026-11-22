# Fact-check: post.md (qc-rejected-codes, Trixagus TGA = Trp)

Checked 2026-11-11, around 18:05. P = project root lab/work/qc-rejected-codes. F = the factcheck scratch dir; my recount scripts are in F (orth.py, codetta.py, fam.py, pseudo.py, pstops.py, rf2aln.py, rfscan.py, pdfgrep.py) and downloads are in F/dl, F/dl2 and F/dl3. Placeholders ([C20 SENTENCE], [image], [TIME]) are ignored.

Summary: 25 items are not OK, made up of 7 WRONG, 13 OVERSTATED, 4 UNSUPPORTED and 1 UNVERIFIABLE. Three are new findings that are not in the README:
- Tmey's and Tduv's RF2 reading frames are broken.
- The Tmey RF2 TGA sits at an *E. coli* Trp position.
- Two 2026 GTDB-team papers are missing from "Known before".

## Lead and intro

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 1 | "symbionts of three *Trixagus* beetles read the stop codon TGA as tryptophan" | OK | P/results/ortholog_tga (recount: 145/159, 151/162, 88/104 at Trp positions), read_support (97.5/92.9/100 %), codetta JSON (W for all 3) | – |
| 2 | "each has deleted RF2" | OK | README 17:33 (flanks 113/514/124 bp vs 997 bp; tblastn finds nothing in Tder's 514 bp at e ≤ 10). The flank distances exist only in the README; I found no results file with them | – |
| 3 | "RF2, the release factor that ends genes at TGA" | OK (minor) | RF2 reads UGA and UAA | optional: "that ends genes at TGA (and TAA)" |
| 4 | "GenBank annotates these genomes with the standard code" | OK | CP157893.2, CP157895.2, CP158689.2: `/transl_table=11` | – |
| 5 | "19-25 genes in each are filed as pseudogenes" | OK | GenBank PGAP notes (F/pstops.py): pseudogenes with internal stops are Tser 21, Tder 25, Ttur 19; every one contains TGA, and 63 of 65 contain only TGA. Total pseudogenes are 27, 44 and 24 | – |
| 6 | "Bacteria end their genes at TAA, TAG or TGA." | OK | – | – |
| 7 | "A few bacteria with tiny genomes, mostly insect symbionts and *Mycoplasma*, read TGA as … tryptophan … and have lost RF2" | OK (minor nuance) | Choi 2026 Table 1: prfB is absent in all except *Providencia*, which has a pseudogene. Other cases are Eggerthellaceae gut MAGs (Parks et al. 2026) and Minisyncoccia (gTranslate 2026) | optional: "…and have usually lost RF2" |
| 8 | "Code changes like this are usually found in the genome sets that GTDB curates." | UNSUPPORTED | *Hodgkinia*, *Nasuia*, *Zinderia* and *Stammera* came from single-genome papers. Shulgina & Eddy 2021 (eLife 10.7554/elife.71402) screened ">250,000 … genome sequences in GenBank". GTDB-based screens are recent: gTranslate (PeerJ 2026, 10.7717/peerj.21825) and Parks et al. 2026 (Microb Genom 10.1099/mgen.0.001767) | "Recent computational screens for code changes run on the genomes GTDB curates." |
| 9 | "GTDB keeps a genome only if its standard-code genes look complete enough" | OVERSTATED | README "GTDB QC numbers": all 7 throscid symbionts fail on marker %, and so do the standard-code non-readers (35.8-36.7 %). README conclusion: "these tiny genomes would fail on marker percentage anyway". qc_bias.tsv: TGA=Trp Mollicutes fail QC *less* often (0.7 %) than bacteria overall (2.7 %) | "GTDB rejects genomes whose marker genes look incomplete. A genome that reads a stop as an amino acid scores lower still, because standard-code gene finders cut its genes short." Then add: "The *Trixagus* symbionts would fail on genome size anyway; their standard-code relatives fail too." |
| 10 | "a genome that reads a stop … has its genes cut short by standard-code gene finders" | OK | symbiont_profile.tsv coding fraction under table 11 vs table 4: Tser 0.744/0.907, Tder 0.665/0.825, Ttur 0.778/0.917. Markers 22.5-25.8 % vs 35.8-36.7 % | – |
| 11 | "So I screened the genomes GTDB had rejected." | OK | – | – |

## Screen bullets

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 12 | "10,593 … failed … R232 quality control: 8,746 good assemblies with low completeness, 1,347 that NCBI excluded …, and 500 at random" | OK | rank_final.txt: T1 8746, T2 1347, T3 500 (sum 10,593); README freeze 11:35 | – |
| 13 | "The test asks whether genes keep matching known proteins past a stop codon." | OK | README method | – |
| 14 | "No genome reads TAA or TAG as an amino acid" | OVERSTATED | rank_final.txt shows TAA flags only for the 2 synthetic controls and TAG flags only for the 4 controls. labels_taa_tag_leads.tsv lists 32 stage-2 TAA/TAG leads (mostly *E. coli* and *Shigella*), and none flags. README lists a recall limit for frameshift-rich assemblies. This is absence of a flag, not proof | "No genome flags for TAA or TAG, although synthetic and known controls…" |
| 15 | "although synthetic and known controls passed through the same pipeline flag" | OK | rank_final.txt; README controls. TAA has only synthetic controls | – |
| 16 | "141 genomes read through TGA." | OK | rank_final.txt: TGA 141 (T1 96, T2 42, T3 3) | – |
| 17 | "103 belong to lineages already known …: *Hodgkinia*, *Nasuia*, *Mycoplasma* and relatives, and Gracilibacteria (glycine)" | OK | labels.tsv: Hodgkinia 49, Nasuia 24, Mollicutes 17, code 25 12, NCBI code 4 1 = 103. These were labelled by NCBI lineage, not checked one by one | – |
| 18 | "Placing the other 38 on the GTDB tree showed that 34 of them belong to known lineages too." | OVERSTATED | README 17:17 / placement.tsv: the method is "nearest representative by summed best bitscore" on 59 markers, which is not tree placement. The *Nasuia* (GCA_964565395.1) placed with a *Zinderia*-like rep at 40.3 % and was identified only by whole-proteome DIAMOND. The 34 are 27 Gly-type + 4 Mycoplasmatales + 1 Minisyncoccia + 1 *Nasuia* + the "Chlamydia", and the Chlamydia is not a known-reader lineage | "Matching the other 38 to GTDB species representatives on 59 marker genes put 33 in lineages with a known code. One '*Chlamydia*' assembly…" |
| 19 | "One '*Chlamydia*' assembly from a cervical sample reads through TGA only in contaminating *Mycoplasma*-like DNA." | OK | README: 247 read-through hits on 151 contigs; top hits are Mollicutes, plus 5 mitochondrial. GTDB lists it as contaminated (CM1 74 %). It places in *Chlamydia* at 99.2 % | – |
| 20 | "Left: three symbionts of *Trixagus* beetles and one of a weevil." | OK | 38 − 34 = 4. This only adds up if the 34 includes the Chlamydia, so make that explicit (see #18) | – |

## *Trixagus* section

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 21 | "*Ca.* Shikimatogenerans is a beetle symbiont with a genome of 130-180 kb and 12-14% GC." | WRONG (for the genus) | symbiont_profile.tsv: *S. bostrichidophilus* is 117-200 kb with GC 13.1-17.2 %; *S. silvanidophilus* JKI-2014 is 308 kb with GC 16.2 %; the throscid symbionts are 133-168 kb with GC 11.8-12.9 %. IMPRS 2025 abstract: "~150 kbp and ~12% GC" | "*Ca.* Shikimatogenerans is a beetle symbiont; in Throscidae its genome is 130-170 kb and about 12 % GC." |
| 22 | "*T. sericeus* from the USA, *T. dermestoides* from Germany, *T. turgidus* from Japan" | OK | GenBank/BioSample host and country: Tser USA (36.08N 121.59W), Tder Germany, Ttur Japan | – |
| 23 | "TGA appears 2.65-3.94 times per 1,000 codons … 65-84% of all tryptophan codons" | OK | symbiont_profile.tsv: 3.40/0.98 (77.6 %), 3.94/0.77 (83.7 %), 2.65/1.43 (65.0 %) | – |
| 24 | "In the symbionts of four other throscid beetles, TGA makes up 1.4-4.8%." | OK | Tmey 1.4, Tduv 2.3, Tcar 3.2, AUS03 4.8 % | – |
| 25 | "the matching gene of a sister symbiont (from *T. carinifrons*)" | WRONG (minor) | ML tree: Tcar is sister to Tder only (UFBoot 50), not to Tser or Ttur | "of a related symbiont (from *T. carinifrons*)" |
| 26 | "has tryptophan at 85-93% of sites, as often as at their own TGG codons (89-91%)" | OVERSTATED (minor) | Recount (F/orth.py): 145/159 (91.2 %), 151/162 (93.2 %), 88/104 (84.6 %) vs TGG 41/45 (91.1 %), 25/28 (89.3 %), 51/56 (91.1 %). Ttur is 6.5 points lower, and the TGG controls are small (28-56 sites) | "about as often as" or "nearly as often as" |
| 27 | "Raw reads carry TGA at 93-100% of covered sites (107-168 sites per genome)" | OK, with a caveat | read_support: Tser 118/121 (depth 4), Tder 156/168, Ttur 107/107. The reads come from the same Illumina metagenome (SRP482428) that the Illumina-only assemblies were built from, so agreement is partly built in | add "(from the same sequencing data the genomes were assembled from)", or note it in Risks |
| 28 | "A blinded re-check. A separate agent ran Codetta … on 15 genomes it knew only by number." | OVERSTATED | README 16:33 says "labels G01-G15", but no G-label mapping exists in P. src/codetta_check.py takes accessions, and the outputs are results/codetta/GCA_*.json (16:22-16:32). The FASTA headers it downloads name the organism ("MAG: Candidatus Shikimatogenerans sp. Tser chromosome"), and the GFFs carry the host. Codetta's call is automatic anyway | "A separate agent, without the README or my results, ran Codetta…" (drop "knew only by number"), or keep "blinded" and add that the call is automatic at the tool's fixed threshold |
| 29 | "Codetta … infers a genome's code from conserved protein domains with no gene calls" | OK | README: six-frame translation against Pfam | – |
| 30 | "called TGA = tryptophan in all three, in a *Nardonella* … and in the *Hodgkinia* and *Nasuia* controls" | OK | codetta JSON: Tser n 17 (margin 35.8), Tder 16 (36.7), Ttur 14 (46.6), Nardonella MH 28 (100.4), Hodgkinia 35 (67.7), Nasuia 33 (86.2) | – |
| 31 | "no TGA call for any of the other nine genomes, the relatives below among them" | OK | 0 columns for Tcar, Tduv, Tmey, AUS03, S.bost, S.silv, Nardonella Rf and Pinf, and Sulcia | – |
| 32 | "RF2 is the only gene present in all four of the other throscid symbionts and missing from all three readers." | OK (exploratory) | gene_content.tsv has 1 row (Tmey 0_26, RF2, best hit RF2_ALIB4, e 1.4e-18). README says it is exploratory | – |
| 33 | "In each reader it is deleted outright; the genes on either side … 113-514 bp apart." | OK (README only) | README 17:33. No separate results file holds the flank coordinates | – |
| 34 | "The four other symbionts keep RF2" | OVERSTATED (new) | Only Tcar (XBT18661.1, 345 aa) and AUS03 (XBT18921.1, 335 aa) have a full-length RF2 reading frame under code 11. Tmey: the table-4 gene 0_26 (343 aa) has an in-frame TGA at codon 129, and GenBank CP146597.1 annotates only a 208-aa C-terminal "peptide chain release factor-like protein" (WXB47382.1) that starts at the Met just after it. Tduv: XBT18340.1 is 207 aa and starts at the same Met ("MYKLWC…"), with TAA/TAG stops in frame upstream (F/rfscan.py) | "All four still carry an RF2 gene, but in two (Tmey, Tduv) the N-terminal part is broken by a stop, leaving a ~207-aa fragment; and in every copy the SPF motif … is gone." |
| 35 | "in every copy the SPF motif … is gone, in three of them by a short deletion" | OK (README only) | README: Tcar "F--", Tduv "F--", Tmey "S--", AUS03 "SSI". The MAFFT alignment is not saved in results/, and AT-rich low-complexity sequence makes gaps uncertain | optional: "in three of them the alignment has a gap there" |
| 36 | "A more distant relative, from a bostrichid beetle, keeps an RF2 with SPF." | OK | README; GCA_029851105.1, host *Dinoderus porcellus* (Bostrichidae) | – |
| 37 | "None of these five relatives uses TGA as tryptophan in conserved protein sequence: of their 34 in-frame TGA codons, 1 sits inside a stretch that aligns to known proteins." | OVERSTATED | My tga_sites.tsv recount gives 34 sites with 1 inside an alignment (Tmey RF2, codon 129). My alignment (F/rf2aln.py) puts that TGA at *E. coli* RF2 **Trp144** (Tmey "GGNDSYNW·IYLLYRMY" / *E. coli* "GGTEAQDWASMLERMY"); Tcar has Leu there. So the only testable site sits at a Trp position, and under code 11 it truncates Tmey's RF2 | "Only 1 of their 34 in-frame TGA codons sits in conserved sequence: inside Tmey's own RF2, where it breaks the standard-code reading frame." |
| 38 | "Most of the rest sit at gene ends, with a TAA or TAG close behind." | OK | 22 of 33 sit past the alignment end. A TAA/TAG follows 1-56 codons later, and within 10 codons for 19 of the 22 | – |
| 39 | "On a tree of 99 proteins …, the three readers do not group together." | OK | concat: 99 families (F/fam.py), 34,439 columns, 10,080 parsimony-informative. ML tree `(S_bost,(Ttur,(Aul,(((Tcar,Tder)50,Tmey)100,(Tser,Tduv)95)100)93))` | – |
| 40 | "(from *T. carinifrons* and a throscid the record does not name)" | WRONG | BioSample SAMN40188933 (Tmey, CP146597.1) is named "Shikimatogenerans throscidophilus Tmey of Trixagus meybohmi" and was derived from a *T. meybohmi* metagenome. src/throscid_figure.py uses the same wrong label ("host not named") | "(from *T. carinifrons* and *T. meybohmi*)". Fix the figure label too |
| 41 | "the *T. sericeus* symbiont is sister to the one from *T. duvalii*" | OK | UFBoot 95 | – |
| 42 | "Forcing the readers into one group makes the tree far worse (AU test, p = 5e-58 …)" | OVERSTATED (precision and a missing risk) | au/test.iqtree: logL −311,572.8 vs −315,531.2, ΔlogL 3,958.4, p-AU 5.49e-58. A RELL-based p this small is not meaningful. concat.faa.log: "6 sequences failed composition chi2 test", so compositional heterogeneity is not addressed. The LG+C20 rerun is still running (iqtree3 PID 56481, 58 min at 18:03), and C20 models site heterogeneity, not across-lineage composition | "(AU test, ΔlogL ≈ 3,960, p < 1e-10)". Add the composition failure to Risks |
| 43 | "the threshold was set before the test" | OVERSTATED | au/CRITERION.txt was written at 17:27, but the ML tree (concat.faa.treefile and .iqtree) finished at 17:03 and already showed the readers apart. Only the α = 0.05 rule came before the AU run | drop the parenthetical, or say "the AU test was set up after the ML tree showed the readers apart" |
| 44 | "So RF2 was deleted, and TGA taken up as tryptophan, at least twice." | OK (inference) | Rests on Tder in the (Tcar,Tder,Tmey) clade (UFBoot 100) and Tser with Tduv (UFBoot 95). Parks et al. 2026 make the same "at least twice" argument for Eggerthellaceae, which is worth citing | – |
| 45 | "All three records, and NCBI's taxonomy for them, use the standard bacterial code" | OK | `/transl_table=11`. NCBI Taxonomy gives gcode 11 for taxids 3158568, 3158566, 3158569 and 3158565 | – |
| 46 | "NCBI's pipeline files 19-25 genes per genome (11-14%) as pseudogenes with internal stops, against none in the *T. carinifrons* symbiont" | OK (state the denominator) | 21/25/19. Internal-stop ÷ (CDS + internal-stop) gives 10.8/13.7/11.0 %. Internal-stop ÷ (CDS + all pseudogenes) gives 10.4/12.4/10.7 %. Tcar has 0, and so do Tduv, Tmey and AUS03 | optional: "(10-12 % of protein-coding genes)" |

## Weevil paragraph

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 47 | "*Nardonella* from the palm weevil *Metamasius hemipterus*" | WRONG (minor) | EPPO and UF/IFAS call it the West Indian cane weevil or silky cane weevil; it attacks sugarcane, banana and palms. "Palm weevil" usually means *Rhynchophorus*, which is the source of the comparison *Nardonella* (GCA_004296495.1), so the label is confusing | "from the weevil *Metamasius hemipterus*" |
| 48 | "(GCA_023651715.1) has no RF2" | OK | README: after the relaxed recheck, tblastn of the sister's RF2 hits only the MH RF1 locus | – |
| 49 | "the blinded Codetta run called its TGA tryptophan" | OK (but see #28 on "blinded") | codetta JSON: W, n 28, margin 100.4 | drop "blinded", as in #28 |
| 50 | "At its TGA sites a related *Nardonella* has tryptophan 72% of the time." | OK | recount: 61/85 (71.8 %); TGG control 7/10 | – |
| 51 | "its genome is a 7-contig Nanopore draft" | OK | symbiont_profile: 7 contigs, 180 kb. Palmieri et al. 2022 MRA (10.1128/mra.00738-22): Nanopore, Canu/Flye, Illumina used for Polypolish | – |
| 52 | "reads cover only 25 of its 141 TGA sites (all 25 carry TGA)" | OK | read_support/GCA_023651715.1_SRR21424116.txt | – |
| 53 | "It stays a candidate." | OK | README: rule outcome "insufficient coverage" | – |

## Known before

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 54 | "TGA = tryptophan is known in *Mycoplasma* and its relatives, some Gracilibacteria-related lineages (glycine), and a handful of insect symbionts." | WRONG / incomplete | Glycine is not tryptophan, so the sentence contradicts itself. It omits UGA=Trp in **Actinomycetota** (Eggerthellaceae gut MAGs with RF2 loss, "at least twice"; Parks, Chaumeil, Chuvochina & Hugenholtz 2026, Microb Genom, 10.1099/mgen.0.001767), which is the closest precedent. It also omits UGA=Trp in **Patescibacteriota/Minisyncoccia** (Chaumeil, Hugenholtz & Parks 2026, PeerJ, 10.7717/peerj.21825), which the screen itself used to call GCA_035453665.1 known. The stop-codon-census notes (RESULT.md:556) already cite Parks 2026 | "TGA as a sense codon is known in *Mycoplasma* and relatives, in Gracilibacteria-related lineages (as glycine), in some Minisyncoccia (gTranslate 2026), in Eggerthellaceae gut bacteria that also lost RF2 at least twice (Parks et al. 2026), and in a handful of insect symbionts." |
| 55 | "Choi et al. 2026 (*ISME J*) list six animal endosymbionts: …" | OK (nuance) | 10.1093/ismejo/wrag259, PMC13712303 Table 1 is titled "Features of **select** animal symbionts using the genetic code with UGA: stop -> tryptophan" and lists those six | optional: "list six animal endosymbionts (a selection)" |
| 56 | "*Cerorastrobacter* is the only one in the phylum Bacteroidota, the phylum of Shikimatogenerans." | OK | Table 1: Cerorastrobacter is in Flavobacteriia. GenBank lineage of Tmey/Tcar: Bacteroidota; Flavobacteriia | – |
| 57 | "Wierz, Malm, Kaltenpoth and Engl (Max Planck Institute for Chemical Ecology) sequenced them." | OK | GenBank authors and affiliation (Dept. Insect Symbiosis, MPI-CE). The abstract adds JGU Mainz | – |
| 58 | "Their paper … is listed as unpublished in the GenBank records" | OK | Reference title "Diversity, functionality, and evolutionary history of bacterial symbionts in false click beetles (Coleoptera, Throscidae)", Unpublished. Europe PMC today finds no published version | – |
| 59 | "a 2025 conference abstract on the work does not mention the code" | OK | IMPRS Symposium 2025 program, p. 30, poster 5, "Small and strict or large and loose? Symbioses in false click beetles" (Wierz, Malm, Engl, Kaltenpoth): "~150 kbp and ~12% GC", nothing on the code. https://www.ice.mpg.de/492439/imprs_symposium_2025_program_and_abstracts1.pdf | – |
| 60 | "I found no report of a code change in Shikimatogenerans or *Nardonella*." | OK | Europe PMC rerun today: Shikimatogenerans AND code terms gives 4 hits, none on code. Nardonella AND code terms gives 1 hit (Cixiidae). Nardonella AND Metamasius gives 2 hits, none on code. Wierz author search shows nothing relevant | – |

## Done / Not known

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 61 | "Tiers, test and pass bar were written down before scoring, with controls that had to pass first." | OK (nuance) | README freeze 11:35; tier files 11:26; rank_final 16:06. The stage-2 amendment and control gate came at 11:31. The near-miss list was added at 11:52 after seeing 2 misses (the bar itself is unchanged) | – |
| 62 | "The checks on the *Trixagus* genomes (orthologs, reads, Codetta) each had their pass rule written before they ran." | UNSUPPORTED (Codetta) | README 16:19 fixes a numeric rule for reads and a qualitative prediction for orthologs ("predicts mostly TGG/TGA there"). No Codetta pass rule is written anywhere; the only bar is "the tool's own 0.9999 bar", named in the results | "The read check had a numeric pass rule and the ortholog check a written prediction before they ran. Codetta's call uses the tool's own fixed threshold." |
| 63 | "The RF2 gene-content comparison, the tree and its test came after the first three checks." | OK | Tree 17:03, gene_content.tsv 17:24, AU 17:47; the first three checks ran 16:20-16:33 | – |
| 64 | "The gene-content comparison was exploratory; the tree test's pass rule was written before it ran." | OVERSTATED | Same as #43: the rule came after the ML tree had been seen | "…the AU test's threshold was written before it ran, but after the ML tree had shown the readers apart." |
| 65 | "All eight symbionts have one tryptophan tRNA, with the usual CCA anticodon and the same anticodon loop." | OK | trna_trp.txt: CTCCAAA in all. The AUS03 loop parse failed, but the sequence holds CTCCAAA at the same columns | – |
| 66 | "No protein data exist." | OK (as far as public data go) | PRIDE search gives 0 projects for "Shikimatogenerans" and 0 for "Throscidae". Unpublished data cannot be ruled out | optional: "No public protein data exist." |
| 67 | "these genomes cannot say" / "*T. turgidus* … long branch near the root …, less certain" | OK | Ttur terminal branch 0.391 vs 0.10-0.29 for the others | – |
| 68 | "two known TGA readers among them were missed by the bar and only show up as near-misses" | WRONG | rank_final.txt: "TGA: 9 near-miss (not called)", all known readers (*Mycoplasmopsis* ×5, *Hodgkinia* ×3, *Nasuia* ×1). The "two" are the first ones seen at 11:52. This is a limit of the read-through test, not of "the QC screen", and some near-misses are 106-118-contig assemblies, not frameshift-rich single contigs | "**The test's limit.** Fragmented or frameshift-rich assemblies can hide a code: nine known TGA readers missed the bar and appear only as near-misses." |

## Note for Ryan

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 69 | "Tiers, the read-through test and its bar were frozen before scoring, with … controls that had to flag first. TAA/TAG came back 0." | OK | as #12-15 and #61 | – |
| 70 | "Every TGA flag was placed: known lineage, or one of the four cases above." | OVERSTATED | The 103 were labelled by NCBI lineage, not placed. Only the 38 were matched to GTDB reps (method in #18). There are five cases above (Chlamydia, 3 *Trixagus*, *Nardonella*), not four | "Every TGA flag is accounted for: 103 by NCBI lineage, 33 of the other 38 by nearest GTDB representative; the rest are the Chlamydia contaminant and the four cases above." |
| 71 | "Three checks had pass rules written before they ran (16:19): orthologs, raw reads, and the blinded Codetta re-check" | UNSUPPORTED (Codetta) | as #62 and #28 | as #62 |
| 72 | "The tree test's rule was written at 17:27, before it ran." | OVERSTATED | as #43: ML tree done at 17:03 | add "but after the ML tree was seen (17:03)" |
| 73 | "After the fact … the SPF alignment, the TGA-site and tRNA comparisons, and the placement of the 38 were done afterwards and are described as such." | WRONG | The post body labels only the gene-content comparison (and the tree) as after-the-fact. The SPF, TGA-site, tRNA and placement bullets carry no such label. Placement did have controls fixed before its run (README 17:17) | label those bullets "exploratory" in the post, or change the note |
| 74 | "Prior work. Choi … Table 1, Europe PMC and Crossref searches, the GenBank records and a 2025 abstract …" | OVERSTATED (incomplete) | It misses Parks et al. 2026 (Eggerthellaceae UGA=Trp, RF2 loss, "at least twice") and gTranslate 2026 (Minisyncoccia) | add both and cite them in "Known before" |
| 75 | Risks: "The genomes belong to the Jena group … They may well have noticed the code." | OK, add one item | The BioSamples already carry proposed species names ("Shikimatogenerans throscidiphilus Tser …", "…throscidophilus Tmey of Trixagus meybohmi"), so a descriptive paper is clearly in preparation | add: "The BioSamples already carry a proposed species name, so their paper is far along." |
| 76 | Risks: "No protein data exist, so TGA = Trp rests on conservation, reads and Codetta" | OK | as #66 | – |
| 77 | Risks: "My earlier RF2 profile wrongly called RF2 absent … Both are corrected in the README" | OK | README 17:32 correction; the SPF paragraph is withdrawn | add: Tmey/Tduv RF2 truncation (#34) and the Tmey TGA at E. coli Trp144 (#37) |
| 78 | Risks: "The *T. turgidus* branch is long … 'at least twice' claim does not depend on it." | OK, add one item | – | add: "6 of 8 sequences fail IQ-TREE's composition test; the C20 rerun (still running) does not address this." |
| 79 | Recommend: "it meets the TAG standard" | UNSUPPORTED | "The TAG standard" is not defined in the post, the README or commons | name the standard, e.g. "the bar used for the TAG post (…)" |
| 80 | Would care: list, incl. "Parks & Hugenholtz (GTDB, for the QC angle)" | OK | They also wrote gTranslate and the Eggerthellaceae UGA=Trp paper, so the case for including them is stronger | optional: "(GTDB; authors of the 2026 Eggerthellaceae and gTranslate UGA=Trp papers)" |
| 81 | "Nothing has been sent." | UNVERIFIABLE (fine) | – | – |
