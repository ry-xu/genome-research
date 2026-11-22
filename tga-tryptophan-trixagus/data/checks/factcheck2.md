# Fact-check 2: post.md (qc-rejected-codes, *Trixagus* TGA = Trp)

Checked 2026-11-11, 18:20-19:50.

**Abbreviations**
- P = lab/work/qc-rejected-codes
- T = P/results/throscidae_tree
- F2 = this directory, which holds my scripts (ident.py, findgap.py, bs.py, xml2txt.py, junction_reads.sh, junction_check*.py, holes.py, depth_window.py, clip_tails.py, edge_reads.py, rf2_readscan*.sh, rf2_scan_summary.py, rf2_hit_gc.py, lowgc_rf2.py) and downloads (dl_ncbi/, dl_lit/, dl_lit2/, dl_lit3/, dl_tax/, dl_reads/, dl_reads2/, dl_reads_SRR28034348/, dl_reads_SRR28034353/, dl_reads_SRR28034362/)
- S = the session scratchpad. Earlier GenBank flat files are in S/gbk_dl, and the first fact-check is S/factcheck/FACTCHECK.md (FC1).

**Sources read in full text (Europe PMC XML)**
- Parks et al. 2026 Microb Genom (PMC13462651)
- gTranslate, Chaumeil, Hugenholtz & Parks 2026 PeerJ (PMC13695967)
- Choi et al. 2026 ISME J (PMC13712303)
- Shulgina & Eddy 2021 eLife (PMC8629427)
- Choi et al. 2025 MBE (PMC12188211)

## Summary

20 items are not OK. Most serious first:

1. **WRONG (record).** The C20 "iteration-4 tree" is the search's starting tree. C20's own search returned the LG topology (#41, #45, #63).
2. **OVERSTATED.** "Every copy is damaged" (#16). Shared SPF loss also means "TGA taken up at least twice" (#22) need not mean decoding changed twice.
3. **OVERSTATED.** "A second case, not the first" (#33, #68). gTranslate 2026 already shows code-11 genomes nested among code-4 relatives (*Stammera*, a Minisyncoccia genus).
4. **OVERSTATED (repeat of FC1 #28/#49).** "Blinded" and "under coded labels" (#13, #29).
5. **OVERSTATED (repeat of FC1 #62).** The ortholog check had a written prediction, not a pass rule (#36, #44).
6. **OVERSTATED (scope).** "It meets the brief's bar" holds for TGA = Trp, not for the "at least twice" story the Recommend line proposes to post (#69).
7. **WRONG (minor).** *Nardonella* "7-contig": it is 26 contigs in 7 scaffolds. "Has no RF2" in a gapped draft is missing its caveats (#28, #29).
8. **WRONG (minor).** The species name is spelled "throscidophilus"; 6 of the 7 BioSamples spell it "throscidiphilus" (#62).
9. **OVERSTATED.** The intro's mechanism. GTDB's CheckM/CheckM2/GTDB-Tk switch to table 4 by coding density, and these readers would qualify (#4).
10. **UNCLEAR (minor).** Tder's "514 bp" (#15). Reads do not support 298 bp of that junction, so its true length is unresolved; the README's "514 vs 113 bp fits separate events" (README l.606) leans on it. The deletion itself holds: my raw-read protein scan finds 0 RF2 reads in each reader's run, against 27-206 expected. The same scan finds 66 of 77 expected in the Tcar control.
11. Smaller items: #7 (the 103 list), #24 (GTDB wording), #25 (image filename), #31 (Known-before omits protist symbionts), #38 (tRNA context), #41 (Ttur "less certain"), #48 ("all fixed"), #2 and #20/#21 (optional nuances), #3 (Shulgina & Eddy).

## Lead, intro and screen

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 1 | Lead: three *Trixagus* symbionts read TGA as Trp; each has deleted RF2; GenBank uses the standard code; 19-25 genes per genome are pseudogenes | OK (see #15 for Tder's RF2) | Tser/Tder/Ttur PGAP "Pseudo Genes (internal stop)" 21/25/19 (S/gbk_dl/*.gb, header block); `/transl_table=11` | – |
| 2 | "Bacteria end their genes at TAA, TAG or TGA. A few bacteria … read TGA as … tryptophan instead, and have lost RF2." | OK (minor nuance) | Choi 2026 Table 1: *Providencia* keeps a prfB pseudogene. gTranslate 2026: some code-4 *Stammera* assemblies contain prfB | optional "usually lost RF2" |
| 3 | "Recent screens for such code changes run on the genomes GTDB curates." | OK (nuance) | True of gTranslate (GTDB R09/R10, plus GlobDB filtered at CheckM2 ≥ 50 %) and Parks 2026. But the largest Codetta screen, Shulgina & Eddy 2021, used "251,571 assembled genomes from GenBank, including partial assemblies" (F2/se.txt l.45). So QC-failed genomes up to 2020 *were* screened. The *Trixagus* genomes (2024) postdate it | optional: add "(Shulgina & Eddy's 2021 screen took all of GenBank but predates these genomes)" |
| 4 | "a genome that reads a stop as an amino acid looks less complete than it is, because standard-code gene finders cut its genes short" | OVERSTATED | gTranslate 2026 intro (F2/gtrans.txt l.12): CheckM uses a coding-density rule (CD4 − CD11 > 5 % and CD4 > 70 % → table 4), "adopted by … GTDB-Tk, CheckM2 and SingleM". The readers' own coding fractions pass it easily (P/results/symbiont_profile.tsv: Tder 0.665 → 0.825, Tser 0.744 → 0.907, Ttur 0.778 → 0.917). So GTDB's tools probably scored them under table 4. The post's own GTDB bullet (#24) concedes the code is not why they failed | "…can look less complete than it is when gene finders use the standard code (GTDB's tools switch codes only when coding density clearly improves)" |
| 5 | Screen: 10,593 = 8,746 + 1,347 + 500; test = homolog matching past a stop | OK | P/results/rank_final.txt; tiers.tsv | – |
| 6 | No TAA/TAG flags; synthetic and known controls flag | OK | FC1; README freeze/controls | – |
| 7 | "141 … read through TGA. By NCBI lineage, 103 belong to groups already known to read it: *Hodgkinia*, *Nasuia*, *Mycoplasma* and relatives, and Gracilibacteria (which read TGA as glycine)." | OK, list incomplete | P/results/labels.tsv: Hodgkinia 49, Nasuia 24, Mollicutes orders 17, Gracilibacteria/"Ca. Altimarinota" 11 (code 25), **Absconditabacterales 2** (code 25), **Ca. Pinguicoccus supinus 1** (code 4, ciliate symbiont, GCA_015697925.1). Sum 103 | "…and Gracilibacteria and Absconditabacteria (which read TGA as glycine), plus one ciliate symbiont (*Pinguicoccus*)" |
| 8 | "Matching the other 38 … put 33 in known groups too. One '*Chlamydia*' assembly … contaminating *Mycoplasma*-like DNA. Left: three symbionts of *Trixagus* beetles and one of a weevil." | OK | 33 + 1 + 4 = 38; README placement section; FC1 | – |

## The *Trixagus* symbionts

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 9 | Throscid Shikimatogenerans are 133-168 kb and 12-13 % GC; *T. sericeus* USA, *T. dermestoides* Germany, *T. turgidus* Japan | OK | symbiont_profile.tsv rows 34-40: 133,256-168,095 bp, GC 0.118-0.129 (Tcar 11.8 % rounds to 12). BioSample geo_loc (F2/dl_ncbi): USA, Germany, Japan | – |
| 10 | TGA 2.65-3.94 per 1,000 codons; 65-84 % of Trp codons; four other throscids 1.4-4.8 % | OK | profile: Ttur 2.65, Tser 3.40, Tder 3.94. Share 65.0/77.6/83.7 %. Others: Tmey 1.4, Tduv 2.3, Tcar 3.2, AUS03 4.8 % | – |
| 11 | Orthologs: Tcar has Trp at 85-93 % of TGA sites vs 89-91 % at TGG | OK | P/results/ortholog_tga (FC1 recount 88/104, 145/159, 151/162) | – |
| 12 | Raw reads: TGA at 93-100 % of covered sites (107-168 sites) | OK | read_support/*.txt: Tder 156/168 = 92.9 %, Tser 97.5 %, Ttur 100 % | – |
| 13 | "A blinded re-check. A separate agent, given the genomes under coded labels and told not to look them up, ran Codetta on 15 genomes." | OVERSTATED (FC1 #28 not fully fixed) | README l.338-341 (corrected 18:11): "The accessions were in the brief, and the downloaded genome headers can name the organism". "Coded labels" implies the agent did not have the identities; it did | "A separate agent, without my README or results, ran Codetta on 15 genomes (by accession; Codetta's call is automatic)." Drop "blinded" in the heading too |
| 14 | Codetta: TGA = W in the 3 readers, *Nardonella*, *Hodgkinia*, *Nasuia*; no call for the other 9 | OK | P/results/codetta/*.json; README table at l.343 ff. | – |
| 15 | "In each reader it is deleted outright; the genes on either side of it now lie 113-514 bp apart." | OK for "deleted outright" (now checked against raw reads) / UNCLEAR for "514 bp" | Deletion: a raw-read protein scan (DIAMOND blastx code 4 of the first 4 M spots of each run against all 8 Shikimatogenerans proteomes plus E. coli RF2) finds **0 reads** whose best hit is any throscid RF2, against an expected 206 (Tder), 177 (Ttur) and 27 (Tser) at each genome's median protein read density. The positive control (the Tcar run, RF2 present) gives 66 reads, against 77 expected. Junction length: Tder's 514 bp (Illumina SPAdes/MEGAHIT assembly) contains the 10-bp estimated gap at 39,816-39,825, which reads do span (an 11-bp deletion, MAPQ 60). But 39,457-39,754 (298 bp, GC 0.30 against 0.12 genome-wide, IUPAC ambiguity codes) has 0-2 read depth against 10-23 on both sides, and no read pair spans it. Reads stop abruptly at 39,456 and are soft-clipped into poly-T/poly-A sequence. Of the six such zero-depth stretches in the genome, this is the only one outside the gene-free 3′ end (156.8-161.6 kb). Tser and Ttur have no internal gap and short junctions (113/124 bp) | Keep "deleted outright", and consider adding the read evidence: "no RF2 reads in any reader's raw data, while the same search finds them in a relative's". Change "113-514 bp" to "113-124 bp in two; the third junction is not fully resolved", or drop the numbers. README l.606 should not use 514 vs 113 bp as evidence for separate events |
| 16 | "The four other symbionts still carry an RF2 gene, but every copy is damaged. The SPF motif … is gone in all four, and in two the front part of the gene is broken as well." | OVERSTATED ("damaged") | README l.466-480: GGQ intact. Tcar (345 aa) and AUS03 (335 aa) are full length. Only SPF is gone ("F--", SSI), from one MAFFT alignment. Tmey's "front break" is an in-frame TGA (codon 129) that breaks it only under code 11; under code 4 it runs full length (README l.476-479). Whether any copy still terminates is unknown (README l.482-483) | "…every copy has lost SPF … and in two the 5′ part is broken (in *T. meybohmi* by an in-frame TGA)" |
| 17 | "A more distant relative, from a bostrichid beetle, keeps an RF2 with SPF." | OK | README l.468 | – |
| 18 | Relatives: 1 of 34 in-frame TGA in conserved sequence (Tmey RF2); most of the rest at gene ends with TAA/TAG close behind | OK | T/tga_sites.tsv recount: 34 TGA sites in AUS03/S_bost/Tcar/Tduv/Tmey. 1 "inside alignment" (Tmey); 22 "past alignment end" (codons to stop: median 3.5, 18/22 ≤ 7) | – |
| 19 | "On a tree of 99 proteins … the three readers do not group together." | OK | README l.542-552; T/concat.faa.iqtree | – |
| 20 | Tder 84 % identical to Tcar and Tmey vs 74.5 % and 65 % to the other readers; Tser sister to Tduv | OK (nuance) | F2/ident.py: 83.8 / 84.1 / 74.5 / 64.7. Tser–Tduv is UFBoot 95 on the LG tree, but by raw identity Tser is *least* similar to Tduv (72.5 %) among Tder/Tcar/Tmey/Tduv (README l.567) | optional: "on the tree, Tser is sister to Tduv" |
| 21 | "Forcing the readers into one group makes the tree far worse (AU test, p < 1e-10), also under … LG+C20." | OK, with caveat missing | T/au/test.iqtree p-AU 5.5e-58. T/c20au/test.iqtree tree 3 Δ3,503.5, p 2e-53. Under C20 the constrained tree is the LG-optimised one, not re-optimised under C20 (README l.589-592; c20au/CRITERION.txt) | add "(constrained tree optimised under LG)", or leave it to the Note |
| 22 | "So RF2 was deleted, and TGA taken up as tryptophan, at least twice." | OVERSTATED (half) | RF2 deletion ≥ 2 follows from the tree and the no-regain argument (README l.603-611). But all four non-readers share SPF loss, while the bostrichid keeps SPF. That points to one ancestral change in RF2 before the split, so TGA decoding could have begun once (ambiguously) and been completed separately. This is Parks 2026's "hybrid explanation", which they test with tRNA synteny (F2/parks2026.txt l.593). Here all eight tRNA-Trp are identical (CCA), so that test cannot be run | "So RF2 was deleted at least twice. Whether TGA reading began once, in an ancestor whose RF2 had already lost SPF, or separately, is open." |
| 23 | GenBank: records and taxonomy use the standard code; 19-25 pseudogenes (11-14 %), none in Tcar | OK | NCBI Taxonomy efetch (F2/dl_tax/ef.xml): gcode 11 for all 13 Shikimatogenerans taxa. 21/(174+21) = 10.8 %, 25/(158+25) = 13.7 %, 19/(153+19) = 11.0 % (README l.391-397). Tcar 0 | – |
| 24 | GTDB: all seven fail on markers; readers 22-26 % vs 36-37 %; "So GTDB's QC is not what hid this code: these genomes would have been rejected anyway." | OK numbers / UNCLEAR wording | P/data/qc_failed_r232.tsv.gz: markers 22.50/23.33/25.83 vs 35.83-36.67, all `Failed marker percentage`=1. QC did hide them from GTDB-based screens; what is meant is that the code is not *why* they were rejected. Parks 2026 makes the same point for *Hodgkinia*, *Nasuia*, *Stammera* and *Zinderia* (F2/parks2026.txt l.627) | "So the code is not why GTDB rejected them: these tiny genomes fail its marker check anyway, as Parks et al. 2026 note for other code-4 symbionts." |
| 25 | `![…](image.png)` caption: tree, TGA share, RF2 state; readers in orange do not group | OK content / WRONG filename | The file is T/figure.png. Figure labels Tmey "SPF lost, TGA inside" and Tduv "SPF lost, front gone", which matches #16 | fix the path when publishing |

**Read check on Tder's RF2 junction (#15).** Method: SRR28034346, first 4,000,000 spots, streamed. minimap2 -ax sr -k15 -w5 against Tder CP157895.2 plus the Tcar and Tmey RF2 regions (F2/junction_reads.sh, F2/dl_reads/). Result:
- **Gap.** 9 reads span the 10-bp N run directly (for example 60M11D91M at 39,756, MAPQ 60; 13M11D138M at 39,803, MAPQ 60). The gap hides nothing; it is a spacer between two pieces of sequence that reads show to be adjacent.
- **Stretch with no reads.** 39,457-39,754 is a separate problem; see #15 (F2/holes.py, F2/depth_window.py, F2/clip_tails.py).
- **Reads near the RF2 regions.** Reads near the junction, mates included (106 records), have no blastx hit to RF2 above noise (2 hits of 21-28 aa at e = 0.005). The 67 reads that minimap2 placed on the Tcar/Tmey RF2 regions are 25-30-bp AT-rich seed matches; 1 of the 67 has a weak blastx hit to RF2.
- **Whole-run protein scan.** This decides the question (F2/rf2_readscan.sh, rf2_readscan_others.sh, rf2_readscan_control.sh, rf2_scan_summary.py):

| run | own-genome reads | expected RF2 reads (median / 10th pct) | reads best-hitting a throscid RF2 |
|---|---|---|---|
| Tder SRR28034346 | 28,117 | 206 / 131 | **0** |
| Ttur SRR28034353 | 23,660 | 177 / 125 | **0** |
| Tser SRR28034362 | 4,242 | 27 / 15 | **0** |
| Tcar SRR28034348 (control) | 11,895 | 77 / 43 | **66** (64 to Tcar RF2, median 100 % identity) |

The readers' other RF2 hits are not Shikimatogenerans:
- E. coli RF2: Tder 892 reads (median 81 % identity, GC 0.44); Tser 38 (94 %).
- Ttur: 275 low-GC reads (GC < 0.20) whose best hit is the E. coli or *S. bostrichidophilus* RF2. Re-searched against the six RF2s alone, all 275 still best-hit E. coli (median 52 %), not any throscid RF2. They come from an AT-rich proteobacterium in the library (F2/lowgc_rf2.py).
- *S. bostrichidophilus* RF2: a few reads at 37-53 % identity.

## *Nardonella*

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 26 | *Nardonella* of *Metamasius hemipterus* (GCA_023651715.1) | OK | BioSample: host *M. hemipterus*, Brazil | – |
| 27 | "At its TGA sites a related *Nardonella* has tryptophan 72 % of the time." | OK | README l.299-300: 61/85 (72 %), and 7/10 at TGG (70 %) | optional: add "(70 % at its own TGG)" |
| 28 | "has no RF2" | OK as found; caveat missing | README l.489-490: RF1 only after the relaxed recheck. But the assembly is 26 contigs in 7 scaffolds, with 19 gaps (~2 kb of N) (F2/dl_ncbi/GCA_023651715.1.json; S/gbk_dl/GCA_023651715.1.gb). Its coding fraction is only 0.445 under table 11 and 0.546 under table 4, against 0.91-0.93 for other *Nardonella* (profile row 41; README l.293), which signals many indels | "no RF2 is found in its gapped draft" |
| 29 | "the blinded Codetta run called its TGA tryptophan … its genome is a 7-contig Nanopore draft, and reads cover only 25 of its 141 TGA sites (all 25 carry TGA)" | WRONG (contig count) / OVERSTATED ("blinded") / caveat missing | NCBI: 26 contigs, 7 scaffolds (RagTag-scaffolded, README l.320-321). 25/141 is correct (P/results/read_support_nardonella.log). The Illumina run was used for Polypolish, so read agreement is partly built in (README l.309-310, 318-319) | "…a 7-scaffold (26-contig) Nanopore draft with low coding density; the reads (also used to polish it) cover only 25 of 141 TGA sites" |
| 30 | "It stays a candidate." | OK | – | – |

## Known before / Done / Not known

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 31 | "TGA as a sense codon is known in *Mycoplasma* …, Gracilibacteria-related lineages (as glycine), Minisyncoccia (gTranslate 2026), Eggerthellaceae gut bacteria, and in a handful of insect symbionts." | OK, incomplete | gTranslate (F2/gtrans.txt l.132): UGA = Trp in Minisyncoccia family JAKLIH01 and genus GCA-2747955. Confirmed. gTranslate l.10 also lists ciliate symbionts (*Pinguicoccus*, *Organicella*) and a dinoflagellate symbiont (Fastidiosibacteraceae). The screen itself flagged *Pinguicoccus* (#7) | "…in a handful of insect symbionts, and in symbionts of ciliates and dinoflagellates" |
| 32 | Choi et al. 2026 list six animal endosymbionts (*Cerorastrobacter*, *Hodgkinia*, *Nasuia*, *Zinderia*, *Stammera*, *Providencia*); only *Cerorastrobacter* is Bacteroidota, the phylum of Shikimatogenerans | OK | Choi 2026 Table 1 (F2/choi.txt l.55-110): exactly those six, Cerorastrobacter class Flavobacteriia. Text: "only six endosymbionts, including this first case in Bacteroidota" (l.121) | – |
| 33 | Parks 2026: "TGA read as tryptophan, RF2 lost, and the change made at least twice, with a lineage that still stops at TGA between the readers. The *Trixagus* symbionts are a second case, not the first." | OK for Parks / OVERSTATED for "second" | Parks abstract (F2/parks2026.txt l.6): UGA = Trp in Eggerthellaceae, "loss of release factor 2 (prfB)", "occurred at least twice", Tapirivita (code 11) between Equivita and Gorillivita. All fair. But gTranslate 2026 (l.7, 92, 134) already reports (a) a basal code-11 lineage of *Stammera capleta* that keeps prfB, inside a code-4 species, and (b) three code-11 genomes "in a lineage otherwise comprised of" code-4 genomes in Minisyncoccia genus GCA-2747955. Neither paper claims "at least twice", but an ordinal "second" is not established | "Parks et al. 2026 made the same at-least-twice argument for Eggerthellaceae, and gTranslate 2026 shows code-11 genomes nested among code-4 relatives in *Stammera* and a Minisyncoccia genus. The *Trixagus* symbionts add an insect-symbiont case." |
| 34 | Wierz, Malm, Kaltenpoth, Engl (MPI Chemical Ecology) sequenced them; paper listed as unpublished; 2025 abstract does not mention the code; no report of a code change in Shikimatogenerans or *Nardonella* | OK | CP157895.2 REFERENCE 1 "Diversity, functionality, and evolutionary history of bacterial symbionts in false click beetles", JOURNAL Unpublished, MPI CE Jena. Europe PMC "Shikimatogenerans" gives 15 hits, and with (UGA OR "genetic code" OR prfB) 3 hits: none report a code (Choi 2025 MBE checked). Crossref and Europe PMC title searches find no published throscid-symbiont paper | – |
| 35 | Done: tiers, test and bar written before scoring, controls first | OK | README freeze 11:35, control gate 11:31 | – |
| 36 | "For the *Trixagus* genomes, the ortholog and read checks had their pass rules written before they ran; the Codetta re-check used the tool's own call rule." | OVERSTATED (ortholog) | README l.254-262 (16:19): the reads had a numeric rule. The ortholog check had a qualitative prediction ("TGA = Trp predicts mostly TGG/TGA there … A frameshift/indel artefact predicts no consistent residue"), with no threshold | "the read check had a numeric pass rule and the ortholog check a written prediction, both before they ran" |
| 37 | Done afterwards, without pass rules: matching of the 38, RF2, TGA-site, tRNA, identity, tree; the tree test's rule written before it ran but after the tree | OK | README timestamps: RF2 17:33, tree 17:48, AU criterion 17:27:57 after ML tree 17:03, identity 18:11 | – |
| 38 | Not known: one tRNA-Trp, CCA, same anticodon loop in all eight; "The Eggerthellaceae readers carry a tRNA with a UCA anticodon; these do not." No public protein data | OK, context missing | T/trna_trp.txt: CTCCAAA in all 8. Parks Table 2: tRNA-Trp(UCA) in 32 of 34 MGR-GC4 MAGs (most also keep CCA, 27/34). Missing context: *Zinderia*, *Stammera* and *Providencia* read UGA as Trp with only tRNA-Trp(CCA) (Choi 2026 Table 1; gTranslate l.92). So "no UCA tRNA" is a known pattern | add "as in *Zinderia*, *Stammera* and *Providencia*" |
| 39 | "Whether the relatives' RF2 still works. Every copy has lost SPF, and two are broken at the front." | OK | README l.466-483 (Tmey: broken under code 11 only) | – |
| 40 | "The order of events … these genomes cannot say." | OK | – | – |
| 41 | "The *T. turgidus* symbiont sits on a long branch, and neither model settles where it joins the tree, so whether it is a third, separate loss is less certain." | OK as stated; see #63 / UNCLEAR | C20's own search returned the LG placement, not the alternative (#63). The alternative (tree 2) was never a C20 result, and AU p 0.109 only fails to reject it. README l.609-611 says both open placements imply a third loss | optional: "neither model rejects the other placement; both make Ttur a third loss, so only a Ttur grouped with another reader, which nothing supports, would not" |
| 42 | "nine known TGA readers missed the bar and show up only as near-misses" | OK | rank_final.txt: 9 TGA near-misses | – |

## Note for Ryan

| # | claim | verdict | evidence | fix |
|---|---|---|---|---|
| 43 | Screen frozen with controls; TAA/TAG 0; every TGA flag accounted for (103 + 33 + *Chlamydia* + 4; one *Nasuia* by protein search) | OK | 103 + 33 + 1 + 4 = 141 | – |
| 44 | "Two checks had pass rules written before they ran (16:19): orthologs and raw reads." | OVERSTATED | as #36 | "the read rule and the ortholog prediction were written at 16:19" |
| 45 | AU rule 17:27 after ML tree 17:03; "The C20 version's rule was written at 18:15, before it ran." | OK, but see #63 | au/CRITERION.txt 17:27:57; c20au/CRITERION.txt 18:15:27. The c20au rule mislabels its tree 2 as "the C20 search's best tree at iteration 4" | – |
| 46 | After the fact: RF2 comparison, SPF alignment, TGA-site, tRNA, identity comparisons and matching of the 38; "the post says so" | OK | post l.30 | – |
| 47 | Prior work: Choi 2026 Table 1, Parks 2026, gTranslate 2026, Europe PMC/Crossref, GenBank, 2025 abstract; "NCBI Taxonomy gives code 11 for all of these taxa" | OK | F2/dl_tax/ef.xml: all 13 Shikimatogenerans taxa gc 11 | – |
| 48 | "A fresh agent checked … found 25 problems, all fixed here." Examples: Parks 2026 missed, *T. meybohmi* host, the two broken RF2 genes | OVERSTATED | FC1 summary: 25 non-OK, confirmed. Not all are fixed: FC1 #28/#49 ("blinded": still "under coded labels" and "the blinded Codetta run", see #13/#29) and FC1 #62 (ortholog "pass rule", see #36/#44). The examples are correct (FC1 #40, Parks, RF2 5′ breaks) | "…found 25 problems; this draft fixes them except the blinding wording" (or fix those) |
| 62 | "The BioSamples already carry a proposed species name ('Shikimatogenerans throscidophilus')" | WRONG (spelling) | F2/dl_ncbi BioSample titles: 6 of 7 throscid samples (Tser, Tduv, Tder, Tcar, AUS03, Ttur) read "Shikimatogenerans **throscidiphilus** …". Only Tmey's sample name has "throscidophilus" (FC1 #40/#75) | "('Shikimatogenerans throscidiphilus')" |
| 63 | "The *T. turgidus* branch is long, and its position may be an artefact; C20 does not settle it either (AU p 0.11 for the alternative). The 'at least twice' claim does not depend on it." | OK in words / rests on a WRONG record | T/c20/concat.faa.treefile was written 17:17:07, the same time as the ckp (curIteration 0, finishedModelInit), and has Ttur sister to the other *Trixagus*. Its logL under C20 in c20au is −310,811.40, the same as the log's initial "Optimal log-likelihood: -310811.498" (log l.48). The search then found "BETTER TREE … iteration 1: -310765.610" and "iteration 4: -310765.122" (log l.71, 74). That matches c20au **tree 1, the LG topology** (−310,765.39). So README l.581-584 ("Its best tree at iteration 4 was … Ttur sister to the other *Trixagus*") and c20au/CRITERION.txt ("the C20 search's best tree at iteration 4") are wrong: tree 2 is the start tree. C20's partial search supports the LG placement. The AU numbers themselves are right | Note: "C20's partial search returned the LG tree; the other Ttur placement (the search's starting tree) is not rejected (AU p 0.11)." Correct the README C20 paragraph, the table row label "C20 iteration-4 tree" and c20au/CRITERION.txt line 2 |
| 64 | Risk: Jena group, unpublished, may have noticed the code, scooping | OK | #34 | – |
| 65 | Risk: no protein data; reads from the same sequencing, so read support is partly built in; ortholog and Codetta do not share that | OK | README l.309-310 | – |
| 66 | "6 of the 8 sequences fail IQ-TREE's composition test … The 'at least twice' claim also rests on raw protein identity (84 % vs at most 74.5 %), which needs no model." | OK | T/concat.faa.log (only Tser, Tmey pass); F2/ident.py 83.8/84.1 vs 74.5/64.7 | – |
| 67 | Risk: earlier RF2 profile wrongly called RF2 absent; earlier SPF suggestion; the RF2 story was built after looking | OK | README l.248, l.470-471 | – |
| 68 | Recommend: write to the Jena group first; post as a second case of repeated RF2 loss after Parks 2026, not a first | OVERSTATED ("second case") | as #33 | "post as another case of repeated RF2 loss, after Parks et al. 2026" |
| 69 | "It meets the brief's bar: raw reads, controls fixed before looking, and an independent re-check by a different method." | OVERSTATED (scope) | The brief (lab/agents/organizer/brief-codon-free.md item 5) asks for "a clear signature, checked against raw reads, with controls that are fixed before you look, and an independent re-check by a subagent using a different method". TGA = Trp meets it. The RF2-deletion and "at least twice" findings, which the Recommend line proposes to headline, are post hoc (post l.30; Note l.49) and were not part of the Codetta re-check | "The TGA = Trp call meets the brief's bar …; the repeated-loss story is exploratory." |
| 70 | Would care: Engl, Kaltenpoth, Wierz, Malm; Parks, Chaumeil, Chuvochina, Hugenholtz; Choi et al.; Shulgina & Eddy; NCBI annotation staff | OK | Parks 2026 author list (F2/dl_lit/parks_meta.json); gTranslate authors Chaumeil, Hugenholtz, Parks | – |
| 71 | "Nothing has been sent." | OK (not checkable; no outbox or mail action seen) | – | – |

Row numbers 49-61 are unused, to keep the Note rows' numbers fixed from my working draft.

**Which FC1 fixes held.**
- **Held (FC1 row numbers):** #9, #14, #18, #21, #25, #26, #37, #40, #42, #43, #47, #64, #68, #70, #72, #73, #74.
- **Held, with a smaller issue left:**
  - #8: the GTDB-screens sentence, now a nuance (here #3);
  - #54: Known before still omits the protist symbionts (here #31);
  - #79: the bar is now named, but its scope is the problem (here #69).
- **Not held:**
  - #28/#49: "blinded" and "coded labels" (here #13, #29);
  - #62/#71: the Codetta part is fixed, but the ortholog check is still called a "pass rule" (here #36, #44);
  - #34: "keep RF2" was replaced by "every copy is damaged", an overstatement in the other direction (here #16).
- **#81:** fine.
