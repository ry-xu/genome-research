# Pre-post check: AGA read as alanine in basidiomycete mitochondria (new-stops-euk) (independent check, 2026-11-11)

Claim checked, from the outbox draft `outbox/2026-11-11-mito-aga-alanine/post.md` (written 20:44):
- In the mitochondria of about 30 basidiomycete fungi (*Trichosporon* and relatives, *Tremella*, *Tilletia*, *Cystobasidium*), the codon AGA most likely means alanine, not arginine.
- "No known fungal genetic code reads AGA as alanine."

What I read:
- the post and `image.png`;
- the README section "Lead (2026-11-11): mitochondrial AGA read as alanine";
- the code in `src/codes_census/` (`sitescan_org.py`, `codon_contrast.py`), `src/aga_ala/` (`trna_vs_usage.py`, `all_columns.py`) and `src/prot/` (`build_db.py`, `run_sage.sh`, `tally_aga.py`, `pep_in_genome.py`);
- everything in `results/aga_ala/`: tables, `audit/`, `blind/`, `reads/` and `prot/setA/`;
- the full set-A Sage output in the scratchpad, which I read without changing it.

How I checked:
- I wrote my own loader for the census site table and re-derived every number from it. The scripts are in `prepost/aga-ala_work/`.
- A separate agent ran the literature search. I re-read the key papers myself: Ahrendt 2026 (PMC full text), Barcytė 2025, Noutahi 2019 (abstract) and Turmel 2020 (PMC full text).
- I did not re-run the proteomics search or re-map any reads.
- I edited nothing else in `~/lab` and messaged no agent.

## Verdict: FIX, THEN POST

- **The signal holds.** Every headline number I re-derived matches (section 2).
  - The per-lineage test survives without the selection rule. I took every genome in each order, whatever its AGA use, and got p 5e-16, 9e-13, 2e-7 and 1e-5.
  - The tRNA is specific. 28 of the 32 genomes carry a tRNA(UCU) with G3:U70 and A73; none of the 4,322 arginine-reading genomes that have a UCU tRNA does.
  - Raw reads support AGA at every site tested.
- **Novelty holds for fungi** (section 1). No paper reads AGA or AGG as anything but arginine in any fungal mitochondrion. The route itself is not new: pedinophyte green algae read AGA and AGG as alanine through a tRNA(UCU) with G3:U70 and A73. The post already says so.
- **Two "Known before" citations are misdescribed and must be fixed** (fixes 1 and 2):
  - Ahrendt 2026 did not "find" codes 1, 3, 4 and 16. That list is in its introduction, cited from Knight et al. 2001.
  - The pedinophytes read both AGA and AGG as alanine, not "AGG, sometimes with AGA".
  - Noutahi 2019 is about AGG only, read by a tRNA(CCU).
- **Three factual slips in the body** (fixes 3, 4 and 7):
  - The raw-read genomes are not "one from each group".
  - The single conserved AGG position is not "in these genomes".
  - "None needed AGA = arginine" in the proteome had almost no chance of being otherwise.
- **The per-genome Fisher test is partly circular** (fix 6). The 32 were picked for having no AGA at arginine positions, so the AGA-at-arginine cell is 0 by construction. Swap in the version without the selection: the conclusion is the same, and the claim is then honest.
- **Proteomics:** one peptide is fine to mention, in two or three sentences, as preliminary (section 3).

## 1. Novelty

A separate agent searched by organism and by mechanism (web, Europe PMC, PMC full texts, NCBI). I re-checked the sources marked "read".

- **Ahrendt et al. 2026**, *Nucleic Acids Res.* 54(2):gkaf1419, doi:10.1093/nar/gkaf1419, PMC12802936. I read the full text.
  - "Only codes 1, 3, 4 and 16" is in the **Introduction**, citing Knight, Landweber & Yarus 2001 (*J. Mol. Evol.*). It is background, not a finding of the survey.
  - The tRNA sentence in the paper is: "Except for Chytridiomycota, which lack tRNAArg genes, all Fungi have tRNAArg(UCU)". It also reports "a replacement of tRNAArg(UCU) with tRNAArg(ACG)" at the base of the Ascomycota.
  - tRNAs were predicted with "tRNAscan-SE with the organellar option", so "Arg" is a label taken from the anticodon. That is exactly why a reassigned tRNA(UCU) would be missed.
  - The analysed set is 328 high-quality mitogenomes (352 with outgroups); "~400" comes from the introduction.
- **Barcytė, Žihala & Eliáš 2025**, *PLoS Genet.* 21(10):e1011901, PMC12574857 (read): "the Arg-to-Ala reassignment of AGA/AGG codons was shown to be apomorphic for the whole order Marsupiomonadales".
  - The "apomorphic" wording means the change is shared by the whole order.
  - The marsupiomonadalean tRNA(UCU) carries "the G3·U70 pair and adenosine at the position 73". This is the same route as here, so the post's "the same route" is right.
  - The paper also has dinoflagellate plastids with AGA/AGG "(most likely alanine)".
- **Turmel, Otis & Lemieux 2020**, *Mitochondrial DNA B* 5(1):548–550, doi:10.1080/23802359.2019.1710283, PMC7748506 (read): "Our protein sequence alignments predict that AGR codons are decoded as alanine instead of arginine" (*Marsupiomonas*). This is the first report of AGA = Ala in a mitochondrion, and the post doesn't cite it.
- **Žihala & Eliáš 2019**, *Genome Biol. Evol.* 11(10):2992–3007: "AGG, sometimes together with AGA, encodes alanine … in diverse sphaeroplealeans". This is correct for Sphaeropleales.
- **Noutahi et al. 2019**, *Mol. Biol. Evol.* 36(4):766–783. I read the abstract: "decoding of AGG codons (normally arginine) as alanine … by tRNA(CCU)". It covers **AGG only**, through a CCU tRNA, and is not a source for AGA or for UCU.
- **Ling et al. 2014**, *Nucleic Acids Res.* 42(1):499–508: *Ashbya gossypii* mitochondrial CUU/CUA = Ala, read by tRNA-Ala(UAG) with G3:U70. The post's description is correct. (Table 26, CUG = Ala in *Pachysolen*, is nuclear, so "the only alanine reassignment known in fungal mitochondria" holds.)
- **Li et al. 2023**, *IMA Fungus* 14:7 (*Apiotrichum*): they annotate with code 4 and report "2 tRNAs with different anticodons coding for arginine". This is cited correctly.
  - The other mitogenome papers on these genera do the same: *T. inkin* (Liu & Wang 2021, *Mitochondrial DNA B*) and *Hannaella oryzae* (Li et al. 2021, *Front. Microbiol.*). Each annotates the UCU tRNA as trnR.
- **Codetta-line surveys** (Shulgina & Eddy 2021 and 2023; Melnykov 2026 KACI; Fakih 2026): none covers fungal mitochondria.
- **NCBI tables** (gc.prt v4.6, tables 1–33): no table has AGA or AGG = Ala. AGR is reassigned only to stop, Ser, Gly or Lys. "No NCBI code … reads AGA or AGG as alanine" is correct.
- **Not searched:** the Ahrendt supplement, bioRxiv full texts, and the *Tilletia* mitogenome papers. A fungal case could still sit unnoticed in a supplement, and the Note for Ryan already says so.

**Bottom line:** "No known fungal genetic code reads AGA as alanine" holds. It is a new fungal instance of a known mechanism (tRNA(UCU) gaining the alanine identity elements), on two or three independent fungal branches. Keep the claim to fungi; it is not the first AGA = Ala in any mitochondrion.

## 2. Correctness: the headline numbers

Every number below was re-derived with my own loader (`census_check.py`, `analyze.py`, `trna_extra.py`, `reads_check.py`, `prot_check.py`, `extra2.py`, `extra3.py`), not taken from the README.

| post claim | my result | check |
|---|---|---|
| about 22,700 mitochondrial genomes; 4,638 fungal | 22,698 (RefSeq 18,015, GenBank 1,275, ORF-only 3,408); 4,638 fungal | match |
| 32 genomes, about 29 species | 32 (Trichosporonales 18, Tilletiales 7, Cystobasidiales 4, Tremellales 2 = *Tremella mesenterica* + *Hannaella oryzae*, Sporidiobolales 1 = *Sporobolomyces salmonicolor* QB252587.1); *T. asahii* appears 4 times, giving 29 species; all on table 4 | match |
| 78% alanine; rest mostly Gly and Ser | 124 of 159 AGA codons at strict columns are Ala; the others are G 15, S 7, P 5, L 4, T 2, V 2 | match |
| 1 to 9 Ala positions per genome; 124 sites at 24 distinct positions | 1–9; 124; 24 distinct across all 32 (39 when each lineage counts its own) | match |
| synonyms at AGA's Ala positions 0; at 35 to 53 Arg positions per lineage | 0; 35–53 | match |
| Fisher p 5e-16, 3e-12, 3e-8, 3e-7, 0.028 | per lineage, distinct columns [AGA@Ala, AGA@Arg; syn@Ala, syn@Arg]: Trichosporonales [16,0,0,53] 5.2e-16; Tremellales [11,0,0,49] 2.9e-12; Tilletiales [6,0,0,51] 2.8e-8; Cystobasidiales [5,0,0,50] 2.9e-7; Sporidiobolales [1,0,0,35] 0.028 | match, but see "Is the test circular?" below |
| All positions: 30% of 1,321 AGA at A, 1% at R; CGN at R 74% | AGA: A 0.297, R 0.010 (13 codons) of 1,321; CGN: R 0.740 | match |
| relatives that read AGA as Arg: AGA at R 76% | 5 same-class arginine-reading genomes (MT556837.1, BK072956.1, BK073411.1, BK074202.1, CP096881.1): AGA R 0.763 (59 sites); their CGN 0.713 | match |
| tRNA: 28 of 32 with UCU + G3:U70, all with A73 | 28; I re-read G3, U70 and A73 from the tRNA sequences with 0 disagreements against the stem parse. The 4 without: NC_039747.1 *Hannaella*, NC_068697.1 *Pascua*, QB252587.1 *Sporobolomyces*, BK074121.1 *Salmacisia* | match |
| 0 of 4,322 arginine-reading genomes with a UCU tRNA have G3:U70 | arginine class (≥3 Arg, 0 Ala) is 4,387 genomes; 4,322 of them have a UCU tRNA; 0 have G3:U70 (22 have no ARAGORN output) | match |
| (figure) 4,387 genomes, 0 carry the tRNA | as above | match |
| blinded re-check: 4 of 4 test genomes, 0 of 4 controls | `blind/key.tsv`: flagged NC_054273.1, BK074411.1, NC_009880.1, MZ439919.1; controls BK072956.1, MT556837.1, BK074202.1, BK072594.1; report calls AGA = Ala in exactly the flagged 4 | match; wording caveat (fix 8) |
| audit: 71% vs 93%; AGA at Arg once in 6 test genomes; 17 of 18 in 3 relatives | `audit/report.md`: 92/129 = 71% vs GCN 93%; 1; 17/18 | match (not re-run; numbers read from the audit's own tables) |
| raw reads "one from each group" | BK074411.1 *Tremella*, BK074449.1 *Trichosporon*, BK074391.1 *Tilletia horrida* and BK074121.1 *Salmacisia*: two are Tilletiales and none is Cystobasidiales | **mismatch** (fix 3) |
| *T. mesenterica* 9 of 9 sites, 54 of 54 reads | strict Ala sites: 9, 54/54 (minimum depth 2) | match |
| *T. horrida* 25, *T. asteroides* 25, *Salmacisia* 26 | these counts are AGA at **any** Ala-consensus column, not strict: 25 (2,163/2,170), 25 (1,320/1,322), 26 (531/531). Strict counts are 6 (416/416), 11 (585/586) and 5 (83/83). For *T. mesenterica* the any-column count is 30 (164/164) | numbers match, **units mixed** within one sentence (fix 3) |
| "over 99.7% of reads agreeing" | *T. horrida* is 2,163/2,170 = 99.68% | **mismatch, minor** (fix 3) |
| AGG "at just one conserved position in these genomes (a methionine one)" | AGG sits at 0 strict columns in the 32, and only 12 AGG codons are aligned across all 32. The one strict-Met AGG is in BK074393.1, a *Tilletia* that is **not** among the 32 | **mismatch** (fix 4) |
| shared positions: 9, about 4 expected | 9 columns used in ≥2 subphyla | count matches; the expectation and P ≈ 7e-4 were not re-derived |
| proteome: 1 peptide in 2 spectra; 20 control spectra all cox2; 0 CGN = Ala; 0 AGA = Arg | my tally of the full set-A Sage file at spectrum q ≤ 0.01 (71,436 target and 713 decoy PSMs): class A = 2 (TPAPNPDHVPEDHTR, fractions A13/A14, RT 15.1/15.2 min); needs CGN = Arg = 20, all cox2; needs CGN = Ala = 0; needs AGA = Arg = 0 | match; but "none needed AGA = arginine" has almost no power (fix 7) |
| peptide flanks are CGN; Arg version never matched | in ORF CM001777.1_+1_24946 (638 aa, 18 AGA): CGT before the peptide, AGA at the Ala, CGA after. TPRPNPDHVPEDHTR was in the database and never matched; the missed-cleavage TPRPNPDHVPEDHTRAPK appears only at q 0.228 | match |
| peptide found only on mito contigs of 5 *T. asahii* assemblies, always AGA | 6 hits in 5 assemblies (CBS 8904 JH925097.1; M6011 32.6 kb contig; Z1 62.8 kb scaffold ×2; CBS 2479 CP116791.1; CM197668.1), all CCT **AGA** CCA … | match |

**Is the test circular?** Partly, as written. The 32 were picked for AGA at ≥1 strict Ala column and 0 strict Arg columns, so in the post's Fisher table the AGA-at-Arg cell is 0 by construction. The synonym counts are not selected, but the cell that matters most is. I re-ran the test with no selection: every genome in each order, whatever its AGA use, distinct columns.

| order (genomes) | AGA@Ala | AGA@Arg | syn@Ala | syn@Arg | p (one-sided Fisher) |
|---|---|---|---|---|---|
| Trichosporonales (22) | 16 | 0 | 0 | 53 | 5.2e-16 |
| Tremellales (32) | 11 | 0 | 0 | 55 | 9.3e-13 |
| Tilletiales (8) | 6 | 1 | 0 | 51 | 1.7e-7 |
| Cystobasidiales (5) | 5 | 3 | 0 | 50 | 1.2e-5 |
| Sporidiobolales (21) | 1 | 0 | 0 | 51 | 0.019 |

- The 1 AGA@Arg in Tilletiales is BK074393.1 at cob:115, where the other 7 Tilletiales genomes have CGA or CGG.
- The 3 AGA@Arg in Cystobasidiales all come from MT366950.1, a second assembly of the same species as BK072865.1, which has none. That is worth a look: a duplicate assembly that disagrees at 3 sites suggests an assembly or annotation difference in one of them.
- A second check that doesn't depend on the selection: AGA at columns that are *not* strict (none of which was used to pick the genomes) sits at A 23.1% and R 1.1% (n = 1,162). CGN at the same kind of columns sits at R 48% and K 28%.
- **Verdict:** the signal is not an artefact of the selection, but the post should show the selection-free version.

## 3. Weak parts

What a specialist would push on, most serious first.

1. **AGA is a much noisier "alanine codon" than GCN.**
   - At strict columns, AGA sits at Ala 78% of the time (124/159). GCN in the same 32 genomes sits at Ala 97.5% of the time (1,841/1,889).
   - Glycine alone: 15 of 139 AGA codons at Ala or Gly columns are at Gly, against 27 of 1,868 for GCN (p 3.8e-8).
   - Counting distinct columns per lineage: AGA 39/53 = 74% Ala, GCN 301/317 = 95% (`extra3.py`).
   - In Tilletiales, AGA's distinct columns are 6 Ala and 5 non-Ala, against GCN's 63 and 0. The non-Ala ones are concentrated: nad5:468 (Ser, 7 genomes), nad5:440 (Pro, 5), nad4:368 (Gly, 4), cob:46 (Leu, 4).
   - The post's "the rest mostly glycine and serine, which other alanine codons also stand in for" hides this: the other alanine codons almost never do. Arginine is still clearly ruled out (AGA at Arg ≈ 1%). But an expert will ask whether AGA is read ambiguously, or is still changing, in some lineages. The post should state the comparison (fix 5), not explain it away.
2. **Six genomes rest on a single site.**
   - Four *Cutaneotrichosporon* (which do carry the tRNA), *Pascua* NC_068697.1 and *Sporobolomyces* QB252587.1.
   - The last two have no G3:U70 tRNA either.
   - *Sporobolomyces* has only 3 aligned AGA codons in the whole genome, and its order has 5 in total. Its p 0.028 means nothing after scanning 4,638 genomes. It should not be listed as a lineage with its own p-value (fix 6).
3. **There is no direct evidence that the tRNA is made or charged with alanine.** G3:U70 and A73 are predicted from sequence.
   - Alternatives a reviewer will raise: tRNA editing, or a nuclear-encoded tRNA-Arg imported into the mitochondrion.
   - Barcytė 2025 (Fig. S18) shows a G3·U70 tRNA(UCU) in *Tetradesmus*, which is not thought to reassign AGA. So the tRNA alone is not decisive.
   - What carries this point is the specificity of the tRNA together with the codon placement: 0 of 4,322 arginine-reading genomes, and 28 of 32 genomes picked on codons alone. The post already says "most likely", which is the right strength.
4. **The proteomics is one peptide.** See the subsection below.
5. **How many origins.**
   - The three subphyla share 9 AGA-at-Ala columns. I confirmed the 9 but did not re-derive "about 4 expected" or P ≈ 7e-4.
   - No tree has been built. The post says this under "Not known", which is enough.
6. **The blinded re-check was not fully blind.**
   - Its working folder was named `blind_aga`, and it said so.
   - It had species names and tested all 64 codons.
   - Its test set had no *Cystobasidium*.
   - It also called ATA = Met in the control *Ceraceosorus* (Codetta logP(Ile) −54), which the post does not mention. That is a separate lead for the README, not for this post.
   - The result still stands, but "told only to infer each genome's code" overstates the blinding (fix 8).
7. **Raw reads for *T. mesenterica* are thin.** 9 strict sites with 54 reads in total and minimum depth 2. The other three genomes have minimum depth 11 to 34 at those sites. The *Tremella* reads are consistent but carry little weight on their own.

**Proteomics, set A: is one peptide enough to mention?** Yes, as a short preliminary line. Not as evidence for the conserved-site call.

- **The identification is solid as far as set A can make it.**
  - Two PSMs (fractions A13 and A14), with discriminant scores 1.218 and 0.954.
  - Decoys at or above those scores: 30 and 225 out of 68,589 decoys.
  - Class-A peptides are only 0.07% of decoy matches (45/68,589). That gives about 0.02 and 0.15 expected false alanine-only matches at the two spectra's own scores, and about 0.5 across all of q ≤ 0.01.
  - The best decoy in the alanine-only class (WFRPLGAMTELK) sits at q 0.017, just outside the cut.
  - The arginine version was searchable and never matched.
- **What limits it.**
  - It is one peptide at one site, in an ORF with no known relatives, so it says nothing about the conserved sites behind the call. The post already says this.
  - The two spectra come from adjacent fractions at the same retention time, so they are probably one elution split across two fractions, not two independent observations. The post's "in 2 spectra" is fine, but don't lean on the 2.
  - TMT co-isolation can mislead.
  - The dataset's isolates are not sequenced. The 6-frame check uses other *T. asahii* assemblies, which all agree.
  - Dataset: iProX-hosted ProteomeXchange PXD025163. Per its OmicsDI record, it holds two *T. asahii* strains isolated 15 years apart from one patient. "Clinical isolates" is fine.
- **The negative control does not test anything.** "None needed AGA = arginine" has close to zero power. The only mitochondrial protein detected with any depth is cox2. Its ORF (269 aa) has just two AGA codons, at positions 6 and 254, and the detected cox2 peptides lie at about 176–220. So no detected peptide could have shown AGA = Arg. The claim is literally true but reads as a control that passed (fix 7).

## 4. Image

`outbox/2026-11-11-mito-aga-alanine/image.png`: a scatter plot. x is AGA at conserved Arg positions, y is AGA at conserved Ala positions. Blue points have the G3:U70 tRNA(UCU).

- **Correct.** The 32, 28, 4,387 and 0 all match my counts. The largest x among fungi is 49 and the axis runs to 50, so nothing is clipped. The footnote's definition of "conserved" matches the code.
- **Reads at a glance.** Yes. A blue column on the y-axis against a long grey band on the x-axis carries the claim without the text. It works as the lead image.
- **Two blue points are off the y-axis**, at about (1, 5) and (3, 3).
  - (1, 5) is BK074393.1, a *Tilletia* whose one AGA@Arg is cob:115.
  - (3, 3) is MT366950.1, the second *Cystobasidium* assembly discussed above.
  - They are honest, but next to "32 genomes put AGA only at alanine positions" they look like exceptions nobody explained. Add a few words to the caption: "two more tRNA carriers (blue, off the axis) each put AGA at 1–3 arginine positions".
- **The 4 of the 32 without the tRNA can't be seen.** They are grey points drawn under the blue ones at x = 0. That's fine, because the annotation states "28 of them".
- **Several blue points sit at (0, 0).** These are carriers with no AGA at any conserved column: *Vanrija* ×3, *Auricularia*, *Filobasidium* and *Tremella fuciformis*. They fit the post's "the codon may have disappeared first", but a reader will wonder. Optional: one clause in the footnote.
- **The six grey points at y ≈ 1, x 35–42** are genomes with AGA at both kinds of position. They are honest and need no comment.
- **The figure does not show the synonym control,** which is the actual test. That's acceptable for a lead image, because the body carries it.

## 5. Wording fixes

Required:
1. **Ahrendt sentence (line 18).**
   - Replace: "A 2026 survey of about 400 fungal mitochondrial genomes (Ahrendt et al., *Nucleic Acids Res.*) found only codes 1, 3, 4 and 16 and states that all fungi have tRNA-Arg(UCU)."
   - With: "A 2026 study of 328 fungal mitochondrial genomes (Ahrendt et al., *Nucleic Acids Res.*), which named tRNAs by their anticodon, reports a tRNA-Arg(UCU) in all fungi except chytrids."
   - Why: the code list is cited background (Knight et al. 2001), not a survey result, and the tRNA statement has exceptions. Optionally add "(named from the anticodon, so a reassigned one would still be called Arg)".
2. **Green algae sentence (line 19).**
   - Replace: "AGG, sometimes with AGA, reads as alanine in the mitochondria of some green algae: Sphaeropleales (Žihala & Eliáš 2019, *Genome Biol. Evol.*; Noutahi et al. 2019, *Mol. Biol. Evol.*) and pedinophytes, where Barcytė et al. (2025, *PLoS Genet.*) traced it to a tRNA(UCU) with G3:U70, the same route as here."
   - With: "AGA and AGG both read as alanine in the mitochondria of the pedinophyte green algae Marsupiomonadales (Turmel et al. 2020, *Mitochondrial DNA B*; Barcytė et al. 2025, *PLoS Genet.*), read by a tRNA(UCU) with the alanine marks G3:U70 and A73, the same route as here. In Sphaeropleales it is AGG, sometimes with AGA (Žihala & Eliáš 2019, *Genome Biol. Evol.*), and in some of them AGG is read by a tRNA(CCU) with the alanine mark (Noutahi et al. 2019, *Mol. Biol. Evol.*)."
3. **Raw reads (line 12), and the Note for Ryan (line 35).**
   - Replace: "**Raw reads** (the sequencing data behind 4 of the assemblies, one from each group) show AGA at every alanine site: *Tremella mesenterica* 9 of 9 sites (54 of 54 reads), *Tilletia horrida* 25, *Trichosporon asteroides* 25 and *Salmacisia buchloeana* 26 (over 99.7% of reads agreeing)."
   - With: "**Raw reads** (the sequencing data behind 4 of the assemblies: *Tremella mesenterica*, *Trichosporon asteroides*, and *Tilletia horrida* and *Salmacisia buchloeana* from the *Tilletia* group; none yet from *Cystobasidium*) show AGA at every site where relatives most often have alanine: 30, 25, 25 and 26 sites, with 99.7% to 100% of reads agreeing."
   - In the Note, replace "raw reads for 4 assemblies (one per group)" with "raw reads for 4 assemblies (2 *Tilletia* group, 1 *Tremella*, 1 *Trichosporon*; no *Cystobasidium*)".
4. **AGG (line 28).**
   - Replace: "**AGG**, the other codon read by a UCU tRNA, sits at just one conserved position in these genomes (a methionine one), so its meaning can't be called."
   - With: "**AGG**, the other codon a UCU tRNA can read, is almost absent from these genomes (12 aligned AGG codons across all 32, none at a conserved position), so its meaning can't be called."
5. **Intro (line 5).**
   - Replace: "in which AGA sits where related proteins have alanine (78% of the time; the rest mostly glycine and serine, which other alanine codons also stand in for) and never where they have arginine."
   - With: "in which AGA, at well-conserved positions, sits where related proteins have alanine 78% of the time (the rest mostly glycine and serine; the usual alanine codons, GCN, manage 97%) and never where they have arginine."
6. **Fisher bullet without selection, and drop *Sporobolomyces* (lines 6–7).**
   - In line 6, replace "(that is how the 32 were picked; the next line is the test)" with "(that is how the 32 were picked; the next line tests it without that rule). Six of the 32 rest on a single site."
   - Replace line 7, from "**Against its synonyms:**" through "0.028 (*Sporobolomyces*, 1 genome).", with: "**Against its synonyms, in every genome of each group:** taking all mitochondrial genomes in each order, whatever their AGA use, and counting each position once, AGA sits at 16, 11, 6 and 5 alanine positions and at 0, 0, 1 and 3 arginine positions; the five other arginine codons (CGN and AGG) sit at alanine positions 0 times and at 50 to 55 arginine positions. A one-sided Fisher test gives p 5e-16 (*Trichosporon* and relatives, 22 genomes), 9e-13 (*Tremella* and relatives, 32), 2e-7 (*Tilletia*, 8) and 1e-5 (*Cystobasidium*, 5)."
7. **Proteins (line 13).**
   - Replace "(PRIDE/iProX PXD025163;" with "(iProX, PXD025163;".
   - Replace "As controls, 20 spectra matched peptides that need CGN = arginine (all from the cox2 protein), none matched peptides that need CGN = alanine, and none needed AGA = arginine." with "As controls, 20 spectra matched peptides that need CGN = arginine (all from the cox2 protein), and none matched peptides that need CGN = alanine."
   - Why: no detected peptide covers an AGA, so the AGA = arginine "control" could not have failed.
8. **Blinded re-check (line 10).**
   - Replace: "told only to infer each genome's code"
   - With: "told to infer each genome's code (it tested all 64 codons; its folder name did mention AGA)".

Optional:
- **Drop *Pascua* and *Sporobolomyces* from the count.** Both rest on one AGA codon with no tRNA. That leaves 30 genomes and about 27 species. If you do, update the intro, the figure annotation and the "28 of the 32" lines.
- **Add the selection-free all-columns figure to line 8:** "at positions not used to pick the genomes, AGA still sits at alanine 23% of the time and at arginine 1%".
- **Line 5, "three branches of the basidiomycete tree" followed by four names.** Make it "three branches (subphyla) of the basidiomycete tree: *Trichosporon* and relatives and *Tremella* on one, *Tilletia* on another, *Cystobasidium* on the third". This also matches "three subphyla" in "Not known".
- **Line 3, "NCBI lists 26 alternatives to the standard code".** Table 11 differs from the standard code only in its start codons, so "25 alternatives that change what a codon means" is more exact. Either is defensible.
- **Line 21, "all of RefSeq".** It is 18,015 of 18,031 RefSeq records.
- **Define terms on first use:** tRNA ("the adaptor that brings an amino acid to a codon"), CGN ("N = any base"), and a one-line gloss for "one-sided Fisher test". ARAGORN can stay a name in "Done".
- **Cut "Done" to two lines.** It repeats the evidence bullets almost item for item.
- **Shorten the audit and blinded bullets** to one sentence each. The intro paragraph and the image carry the post.
- **Proteins, "in 2 spectra"** → "in 2 spectra (adjacent fractions)".
- **Note for Ryan, outside the post.**
  - The README's peptide statistics (hyperscore 58.0/41.3, gap 22.9/10.8, near-hit q 0.45) don't match the set-A results file (delta_next 25.8/16.2; the R-version missed-cleavage peptide at q 0.228). They probably come from an earlier run. The post doesn't quote them, but the README should say which run they come from.
  - The blinded agent's ATA = Met call in *Ceraceosorus* is an unrelated lead worth logging.

Not errors:
- "No NCBI code (tables 1–33) reads AGA or AGG as alanine".
- "AGA is arginine in every code used by fungi" (tables 1, 3, 4, 12, 16 and 26 all have AGA = Arg).
- The Ling 2014 sentence.
- The Li 2023 citation.
- "Tilletia (wheat bunt and smut fungi)": the 7 include *T. indica*, Karnal bunt of wheat.
- The Note's risk list, which is accurate.
