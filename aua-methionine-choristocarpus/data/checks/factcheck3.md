# Fact-check 3: database-check edits to outbox/2026-11-14-choristocarpus-aua-met/post.md

Fresh checker, 2026-11-14. Scope: only the sentences changed between
`scratchpad/chor_post_pre_dbcheck.md` and `post.md` (diff: post lines 3, 5, 27-30, 46-48, 74, 75).
Source of truth: `checks/chor_wgs/dbcheck/` (DBCHECK.md, d*.out, trna/, data/). Paths below are relative to `dbcheck/` unless noted.

## Independent recounts I ran (read-only; my scripts are in the session scratchpad `fc3/`)
- **tRNA finder raw files (trna/\*.txt):** CAT-anticodon tRNAs per run. *Choristocarpus* NC_085326 and CAXNAD020000071.1 have 2 in each of the 5 runs (tscanO, tscanG, tscanB, aragorn, aragorn_ps), so 10 of 10 runs. *Discosporangium* CAXNAE020000396.1 has 3 in each of its 5 runs. There are 0 TAT in any of these runs.
- **trnW-trnQ slot:** in all 5 *Choristocarpus* runs, trnW (7563/7564-7636/7637) is followed directly by trnQ (7775). In all 5 *Discosporangium* runs, trnW (21714/21715) is followed directly by trnQ (21920/21921). In the controls, a CAT tRNA sits between trnW and trnQ: Saccharina 8353, Phaeostrophion 22257 and Desmarestia 9300 (d7_cat_context.out:22,30,40).
- **CAT tRNA identity, my own global alignment with end gaps penalised:**
  - *Choristocarpus* CAT1 (2626-2700): Saccharina fMet 63/75 = 0.840, then Desmarestia fMet 0.747.
  - *Choristocarpus* CAT2 (2859-2931): Desmarestia Met-e 51/73 = 0.699, Saccharina Met-e 0.671 and Phaeostrophion Met-e 0.630. Against the three Ile-CAU tRNAs it scores 0.534-0.589. This matches d7.
  - Caveat: d6_trna_table.out:145-147 uses a different metric (matches over aligned columns, with end gaps free). It ranks Desmarestia trnI(cau) first for CAT2 (40/48 = 0.833), because that alignment is shifted. With the raw match counts, Met-e still wins (51 against 40). The post follows the full-length figure, and I think that is right. A reader of DBCHECK, though, will see the conflicting table.
- **NC_085326 / ON552980 GenBank CDS:**
  - 35 CDS each, with no transl_table.
  - Stops: TAA 26, TGA 6 (rps8, rps19, rps11, nad2, cob, nad3) and TAG 3.
  - 0 internal TGA, and 227 internal ATA.
  - Every stored /translation equals the table-1 translation.
  - The Genoscope contig has 6 TGA stops and 229 ATA. This matches d3_genbank.out:9,18,60 and d13_starts.out:1-2.
- **tRNAdb (data/tdb_r2.json):** total 61 entries, all anticodon CAT, all "Met", from 61 organisms, with exactly 1 entry per organism.
- **Same submitters:** the Direct Submission author lists of ON552980 and ON552993 are identical (Choi S.-W., Choi J.W., Graf, Kawai, Jung, Choi C.G., Yoon), both submitted 19-MAY-2022. ON552993 22257..22330 is gene="trnM-CAU" and sits between trnW-CCA and trnQ-UUG (data/ON552993.1.gb:245-253).
- **New finding (not in DBCHECK):** the lab's own GenBank pull (`work/organelle-reading/data_gb/records_mito.parquet`) has 3 brown-algal mitogenomes that declare **transl_table=4** on every CDS (110 CDS in all):
  - PX764281 (*Ishige foliacea*), PX764282 (*Ishige okamurae*) and PX764283 (*Petroderma maculiforme*).
  - Submitted 22-DEC-2025 by Hwang, Kim, Choi, Kawai and Yoon; I checked this in the raw gbff (`data_gb/raw/mito_protist_*.gbff.gz`).
  - Table 4 reads AUA as Ile but UGA as Trp, so it is a non-standard code, though not AUA = Met.
  - The other 358 brown-algal GenBank records and all 96 RefSeq ones declare no table.
- **NCBI table 23:** the table text in data/wprintgc.html says it "has been created for the mitochondrial genome of the labyrinthulid *Thraustochytrium aureum*". So NCBI does not give every labyrinthulid the standard code. Only the Labyrinthulomycetes node (35131) has MGC 1 (d2_codes.out:10).

## Sentence verdicts

### post.md:3
**S1. "NCBI gives *Choristocarpus* the standard code, but it also gives that code to the yellow-green algae and labyrinthulids, so the assignment is a default."**
**FIX (minor).**
- Supported:
  - *C. tenellus* MitoGenCode 1 (d1_taxonomy.out:3; d2_codes.out:2).
  - Xanthophyceae 2833 mito = 1 (d2_codes.out:6).
  - Labyrinthulomycetes 35131 mito = 1 (d2_codes.out:10).
- Two problems:
  - (a) NCBI has a dedicated mito table (23) for the labyrinthulid *Thraustochytrium aureum* (wprintgc.html table-23 text), so "gives that code to the labyrinthulids" is true only of the group's node.
  - (b) "is a default" is an inference. The inheritance flag was not checked (DBCHECK.md:28-29,153). What the evidence shows is that the assignment carries no information about this alga.
- Suggested: "NCBI gives *Choristocarpus* the standard code, but it also gives that code to the yellow-green algae and to the labyrinthulids as a group, where AUA for methionine is published, so the assignment is not evidence either way."

### post.md:5
**S2. "(NCBI has no such order and files the genus under Sphacelariales)"**
**OK.** Discosporangiales has 0 hits as Scientific Name and 0 as All Names (d14_tax_search.out:6-7). The lineage runs "...Phaeophyceae; Sphacelariales; Choristocarpaceae" (d1_taxonomy.out:5,9). Optional: "NCBI Taxonomy" rather than "NCBI".

### post.md:27
**S3. "Reference databases (a separate agent, 14 Nov): none records a non-standard code for any brown alga, and none contradicts this one."**
**FIX.**
- (a) This is wrong as written. GenBank holds 3 brown-algal mitogenomes declared with table 4, a non-standard code with UGA = Trp: PX764281-PX764283 (*Ishige* ×2, *Petroderma*; see the recount above).
- The database check did not survey brown-algal GenBank records beyond *Choristocarpus* and 3 controls (DBCHECK.md:34-36). Its own "or any non-standard mitochondrial code ... any brown alga" (DBCHECK.md:8-9) is also too broad.
- What is supported: no database records AUA = Met for any brown alga.
  - Phaeophyceae node MGC 1 (d1_taxonomy.out:57).
  - No NCBI code table names brown algae (d2_codes.out:44-55).
  - All 61 tRNAdb brown-algal CAT entries are labelled Met (d10_tdb.out:21-23).
  - All 457 brown-algal GenBank/RefSeq mito records declare table 1 or 4. Both tables read AUA as Ile.
- (b) "none contradicts this one" sits badly with the next sentence, which says every copy reads AUA as isoleucine. A lay reader will take that as contradiction. The intended point is that those readings are uncommented defaults.
- Suggested: "Reference databases (a separate agent, 14 Nov): none records AUA as methionine for any brown alga. Those that give a reading for *Choristocarpus* use the standard code with no note or check behind it."

### post.md:28
**S4. "NCBI Taxonomy gives the standard code."**
**OK.** d1_taxonomy.out:3 (MitoGenCode id=1).

**S5. "The 35 genes in GenBank and RefSeq, the 70 NCBI Protein records and the 35 UniParc sequences all read AUA as isoleucine, with no notes."**
**FIX (wording).**
- Numbers and the AUA claim are supported:
  - 35 CDS, no transl_table, 227 ATA, 0 translation mismatches (d3_genbank.out:9,18,51,60; my recount).
  - 70 proteins equal the table-1 translation and 0 equal AUA→Met (d12_proteins.out:1,11).
  - UniParc: 35 entries, 0 with a UniProtKB accession, 35 equal table 1 (d12_proteins.out:22).
- "with no notes" is literally false for NCBI Protein. Each set of 35 records carries 63 automatic CDD /note qualifiers (d12_proteins.out:7,17,29-30; data/nc_proteins.gp has 63 "/note"). The GenBank nucleotide records do have 0 /note (d3_genbank.out:10,52).
- "genes" could be read to include the 22 tRNA genes.
- Suggested: "The 35 protein-coding genes in GenBank and RefSeq, the 70 NCBI Protein records and the 35 UniParc sequences all read AUA as isoleucine, with no note about the code."

**S6. "UniProtKB has no entry."**
**OK.** uniprot_count_116065.hdr x-total-results 0; uniprotkb_q2.hdr 0; uniprot_376270.tsv header only (d12_proteins.out:24-28).

### post.md:29
**S7. "The tRNA databases cannot tell: tRNAdb labels every brown-algal CAU tRNA (61) methionine, and GtRNAdb has no brown algae."**
**FIX.**
- The 61 is right, but "every brown-algal CAU tRNA (61)" reads as if 61 were all the CAU tRNAs of brown algae. In fact tRNAdb holds exactly one CAT entry per species: 61 entries from 61 organisms (d10_tdb.out:21-23; my recount of data/tdb_r2.json gives 61/61, with 1 per organism).
- Standard brown algae have three CAU tRNAs each (d6_trna_table.out:65-132). That is why tRNAdb cannot tell, and the sentence hides the reason.
- The tRNAdb *Choristocarpus* set is also a filtered 13 of 22 (d10_tdb.out:1-17).
- GtRNAdb is supported: 610 eukaryote links, none a brown alga or stramenopile (d11_gtrnadb.out:2-25).
- Suggested: "The tRNA databases cannot tell: tRNAdb keeps only one CAU tRNA per brown alga (61 species) and labels each one methionine, and GtRNAdb has no brown algae."

### post.md:30
**S8. "No NCBI code table fits."**
**OK.** No table differs from table 1 only at ATA (d2_codes.out:77-83).

**S9. "Every table with AUA for methionine also reads UGA as tryptophan, and *Choristocarpus* ends 6 genes with UGA."**
**OK.**
- AUA = Met tables are 2, 3, 5, 13 and 21, and all have TGA:*->W (d2_codes.out:15-29,78-82).
- 6 TGA stops in NC_085326 and in the Genoscope contig (d13_starts.out:1-2; my recount: rps8, rps19, rps11, nad2, cob, nad3).
- My recount also finds 0 internal TGA in all 35 CDS. This supports UGA being a real stop, and the post could say so.

**S10. "Annotating it correctly would need a new table."**
**FIX (minor).** This is true only if the claim holds, and as written it presumes the result.
- Suggested: "If AUA does mean methionine here, annotating it correctly would need a new table."

### post.md:46
**S11. "A separate agent ran two tRNA finders (tRNAscan-SE, ARAGORN) on both genome sequences."**
**OK.** tRNAscan-SE 2.0.13 (-O, -G, -B) and ARAGORN 1.2.41 (default, -ps) were run on NC_085326 and CAXNAD020000071.1 (d4_run_trna.sh, d5_run_trnaB.sh; trna/).

**S12. "Every run finds two CAU tRNAs, one like the start-codon methionine tRNA and one like the other methionine tRNA."**
**OK, with one recommended addition.**
- Supported:
  - 2 CAT in each of 10 of 10 runs (my raw recount; d6_trna_table.out:10-43).
  - CAT1 is 63/75 = 0.840 to the Saccharina initiator (d7_cat_context.out:4).
  - CAT2 is 51/73 = 0.699 to the Desmarestia Met-elongator, against 40/73 = 0.548 to trnI(cau) (d7_cat_context.out:6). My recount with end gaps penalised agrees.
- Note that "like" means modest similarity (0.70). See also the conflicting d6 metric under the recounts.
- Omission: tRNAscan-SE -O also reports a weak "Ile" tRNA with anticodon AAU in both assemblies (score 20.2, malformed anticodon stem, best match 0.640; d6_trna_table.out:13,28; d8_aat_hit.out:1-9; DBCHECK.md:69-73). An AAU-anticodon tRNA bears directly on "which tRNA reads AUA", and anyone who reruns tRNAscan-SE -O will see it.
- Suggested addition after S12: "tRNAscan-SE's organelle mode also reports a weak isoleucine tRNA with anticodon AAU (score 20, malformed). The other runs do not find it, and it is probably not a real tRNA."

**S13. "The place where brown algae keep the isoleucine CAU tRNA (between the tRNAs for tryptophan and glutamine) is empty."**
**OK (minor wording).**
- Supported: trnW is immediately followed by trnQ in all 10 *Choristocarpus* runs (my recount; d7_cat_context.out:2,8).
- The pattern comes from 3 genomes, so optionally say "where the brown algae checked keep...".
- "Empty" means no tRNA is predicted there. The gap is 137 letters, enough room for a degenerate copy that the finders miss.

**S14. "Three standard brown algae have it there."**
**OK.** Saccharina 8353, Phaeostrophion 22257 and Desmarestia 9300 each sit between Trp-CCA and Gln-TTG, with 0.877-1.000 identity to trnI(cau) (d7_cat_context.out:25-26,33-34,43-44).

### post.md:47
**S15. "GenBank labels both *Choristocarpus* CAU tRNAs "trnM"."**
**OK.** gene=trnM-CAU, product=tRNA-Met at 2624:2701 and 2858:2931 (d3_genbank.out:20,22,62,64).

**S16. "That is not evidence: the same submitters also label *Phaeostrophion*'s isoleucine CAU tRNA "trnM"."**
**OK.** ON552993 22257..22330 is gene="trnM-CAU" in the trnW-trnQ slot (data/ON552993.1.gb:247-251). Its Direct Submission authors and date are identical to ON552980's (data/ON552980.gb and ON552993.1.gb REFERENCE blocks). It is identified as Ile-CAU by synteny and 0.877 identity to Desmarestia trnI(cau) (d7_cat_context.out:33-34).

### post.md:48
**S17. "The isoleucine CAU tRNA is also missing from that place in *Discosporangium*, which reads AUA as isoleucine."**
**OK.** Trp 21715 → Gln 21921, with nothing between them in all 5 runs (my recount; d7_cat_context.out:14).

**S18. "*Discosporangium* keeps one extra CAU tRNA that *Choristocarpus* lacks."**
**OK.** 3 CAT in each of 5 *Discosporangium* runs against 2 in *Choristocarpus* (my recount; d6_trna_table.out:44-64).

**S19. "One tRNAscan-SE mode calls it the AUA-reading isoleucine type, but that mode gave the same call for none of the three standard brown-algal isoleucine CAU tRNAs, so the label is weak."**
**OK.**
- The extra tRNA is called "Ile2" by tscanB (trna/CAXNAE020000396.1.tscanB.txt; d6_trna_table.out:55).
- The three control Ile-CAU tRNAs are called "Met", with best isotype Arg (d6_trna_table.out:75,98,122; raw tscanB.txt lines re-grepped).
- Optional strengthener: this mode's highest-scoring isotype model for the extra tRNA is lysine (Lys 67.2 against Ile2 64.5; trna/CAXNAE020000396.1.tscanB.iso.txt row trna20). My identity recount also puts it about equally close to the Met-elongators and the Ile-CAU tRNAs (0.699 against 0.685).

### post.md:74 (Note for Ryan)
**S20. "NCBI still files the genus under Sphacelariales; the post now says so and follows the current order (Discosporangiales)."**
**OK.** NCBI placement is from d1_taxonomy.out:5. "Now says so" is true (post.md:5). "Current order" was not changed in substance, so I did not re-check it.

### post.md:75 (Note for Ryan)
**S21. "Database check (15:22-15:33, separate agent): `checks/chor_wgs/dbcheck/DBCHECK.md`."**
**OK.** DBCHECK.md:3. The files are timestamped 15:22-15:34.

**S22. "NCBI Taxonomy, GenBank/RefSeq, NCBI Protein, UniProt/UniParc, tRNAdb (mitotRNAdb), GtRNAdb and the NCBI code tables; tRNAscan-SE and ARAGORN runs."**
**OK.** DBCHECK.md sections 1-4.

**S23. "Nothing records or contradicts the claim; the post's NCBI, order and tRNA sentences are changed to match it."**
**FIX (minor).**
- "contradicts" has the same problem as S3.
- Ryan should also know what the check did not cover (DBCHECK.md:152-157) and the new table-4 finding.
- Suggested: "Nothing records AUA as methionine for any brown alga. The standard-code readings of *Choristocarpus* are uncommented defaults. Not checked: the NCBI inheritance flag, the Phaeoexplorer portal, ENA and BOLD. Separately, 3 brown-algal GenBank mitogenomes (PX764281-3, *Ishige* and *Petroderma*, Yoon lab 2025) are declared with table 4 (UGA = Trp), so 'no non-standard code in brown algae' is not true of GenBank. The post's NCBI, order and tRNA sentences are changed to match."

## Overall
Overall: the numbers check out, but four sentences need fixing before posting. The most important is S3, which is false as written: GenBank holds 3 brown algae declared with code table 4. S5 ("no notes"), S7 ("every CAU tRNA (61)") and S1 (the labyrinthulid wording) are also wrong or misleading, and S10 and S12 need smaller fixes.
