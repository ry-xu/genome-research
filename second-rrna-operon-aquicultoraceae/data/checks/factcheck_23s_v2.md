<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# Fact-check of outbox draft v5 changes (fresh subagent, 15 Nov)

Brief: check only the new or changed sentences in the 23S post (see notes.md) (a diff in the checker's scratchpad (not kept)), recompute where possible, read-only. Scripts and outputs: the checker's scratchpad (not kept) (scratchpad, session 13993f59). Saved verbatim by Deneb. Below the report: what Deneb re-verified before acting on it.

I found three wrong claims, a few overstated ones, and one gap in the "Checked by me only" list. Most counts reproduce from the raw tables. The biggest error: the UBA6777 "300–370-nt inserts" are runs of N (scaffold gaps), not inserted sequence.

## Claim table

| # | Claim (line) | Verdict | Evidence |
|---|---|---|---|
| 1 | 50,594 genomes scanned (13, 50) | HOLDS | `bac/list`: 50,607 lines, 50,594 unique. The `done.*` files hold 50,592; the 2 others were rescued by hand (`bac/rescue/`). |
| 2 | 4,641 genomes with ≥2 copies of ≥2 kb (50) | WRONG by 1 (trivial) | Rerunning `div23_post.py bac` gives 4,641 and 51 flags, identical to `post_full.tsv`. But the rescued GCF_900199455.1 has 2 full copies (model 9–2922), so the true count is 4,642. |
| 3 | 51 flagged copies in 46 genomes (13, 50) | HOLDS | `class_full.tsv`: 51 rows, 46 accessions. |
| 4 | 11 foreign or pile genomes, 2 family genomes, 33 singletons (51–52) | HOLDS | By genome: foreign 4, ref-foreign 5, rrna-pile 2, singleton 35. The 35 are the 33 rows of `disposition_outfam.tsv` plus DASXVE01 and DASQQM01. |
| 5 | "33 genomes carry a lone divergent copy" (52) | OVERSTATED (minor) | GCA_023250515.1 has 2 flagged copies, both archaeal. Five of the 33 carry foreign-looking copies, which blurs the split from the "11 foreign" group. |
| 6 | 16 canonical / 11 unread / 6 changed (52–53) | Counts HOLD against the TSV | 13 canonical + 2 canonical-foreign + 1 canonical-truncated; 6 truncated + 3 foreign + 2 with no hit to *E. coli* 23S; 3 seq-error + 2 decaying + 1 killed. |
| 7 | "16 are canonical at the PTC wherever it can be read" (52) | OVERSTATED | See the note below this table. |
| 8 | 11 unread: truncated, foreign, or not alignable (52) | HOLDS | Disposition TSV: 6 / 3 / 2. |
| 9 | 3 seq-error genomes are nanopore-only or hybrid (53) | HOLDS | `flag_tech.tsv`, and the NCBI Datasets API (no personal data sent): GCA_041398565.1 and GCA_046351315.1 are Flye on GridION only; GCF_032654805.1 is Unicycler on Illumina HiSeq plus MinION. |
| 10 | 2 decaying *Arsenophonus* copies in insect endosymbionts (53, 98) | HOLDS | Both are *Arsenophonus* MEGAHIT MAGs. Hosts from BioSample: cicada bacteriomes (*Neoplatypleura badia*) and fat bodies (*Eopycna coelestia*). |
| 11 | GCA_046530935: GC drift toward the genome's, 5' third missing (98) | HOLDS | My recomputation: copy GC 0.374, normal 23S 0.501, genome 0.326. The copy starts at model position 834 of about 2,920 (29% missing), 13.7 kb into a 21.4 kb contig, so this is not a contig-edge truncation. |
| 12 | One co-assembly artefact killed by long reads (53) | HOLDS | `data/checks/killed_candidate_GCA_025542975.md`, line 1 ("killed"), long reads SRR14780129–32. |
| 13 | "best 76% over 0.6 kb" (13) | HOLDS, but incomplete | `blastn` of `flag51.fa` against `famdiv.fa`, out-of-family queries only. dc-megablast and blastn give at best 75.95% over 603 nt (GCA_040766325.1); megablast gives 76.6% over 594. No out-of-family copy beats this among hits ≥300 nt. However, the highest-bitscore hits are longer and lower: 70.5% over 1,415 nt and 69.1% over 1,439 nt. `flag51.fa` matches `post_full.tsv` and the `copies.*.fna` sequences exactly. `famdiv.fa` holds 9 loci, not one per divergent genome. |
| 14 | "none of them shows a PTC change that survives a check" (13) | OVERSTATED (scope) | Only 22 of the 44 out-of-family flagged genomes had a PTC reading (16 + 6). The other 22 had none: 11 unread singletons plus the 11 foreign or pile genomes. |
| 15 | The scan cannot see single-copy genomes such as UBA6777 (13) | HOLDS | Of the 5 UBA6777 MAGs, only GCA_002731275.1 is in the scan list, and it has one copy (9–2920). The other 4 are not in the list at all. |
| 16 | PTC bases G2252A, G2583A, U2584A (54) | HOLDS | `data/uba6777/ptc_ux.out`. It also shows A2451C and C2452A, which the notes call aligner-dependent. |
| 17 | 5 MAGs, 4 projects, all from Tara reads (54, 96) | HOLDS | NCBI: PRJNA348753 (×2), PRJNA417962, PRJNA391943, PRJEB51078. ENA: ERR868378 belongs to PRJEB4352 (Tara station 136). |
| 18 | "300–370-nt inserts ... that the reads do not support", in three MAGs (54, 94, 96, 57) | WRONG | These are runs of N, i.e. scaffold gaps. I extracted each locus and ran `blastn` against `ux_div`. All 7 unaligned stretches are N: DKGF01 has 313/369/361 (257 of the 313 are N), DCKP01 344/313/303, DLUY01 363. Genome-wide, DKGF01 has 75 N-runs (25,338 N), DCKP01 132, DLUY01 155; PCBK01 has none of that size. The KILL report's "AT-rich, GC as low as 0.17" came from counting N as AT. Reads failing to span the edges is what an N-gap produces. Only the DKGF01 gaps were ever read-tested. Over the bases they share, DKGF01, DLUY01, CALELE and PCBK01 are 100% identical; DCKP01 is 100% except a 3' stretch at 97.4% over 649 nt. So the "identical" claim was mostly right, and the correction overturned it for the wrong reason. |
| 19 | The read-supported copy is the one in PCBK01, and the PTC calls use it (96) | HOLDS | `ux_div.fa` is the PCBK01 copy. NODE_1 is 100% identical to it over 3,483 nt. `ptc_ux.out` was made from `ux_div`. |
| 20 | NODE_6 "canonical wherever it covers the PTC" (55) | HOLDS, but misleading | My own local alignment to *E. coli*: canonical at 2062, 2251–2253, 2449–2453, 2503 and 2504. It does not reach 2583–2584. Those are 2 of the 3 positions that are changed in the divergent copy, so it is compared at only one of the three. Only the author checked this; KILL §3 says the kill agent did not check NODE_6's PTC calls. |
| 21 | ~3× depth across 3 stations, tighter than all 34 decoys (55) | HOLDS | `coab.out`: N6/H 3.006, 3.110, 3.024; 34 of 34 decoys. ENA sample titles: stations 128, 138, 142. I recounted ERR599136 myself: 845 NODE_6 reads at ≥0.98, matching `coab.out`. |
| 22 | 84% identical to "the normal 23S of two other genomes in the host's family" (56) | OVERSTATED | My `blastn` gives 83.6% / 83.5% / 82.3% to the three loci (`fam_plain.fa`). The two genomes, GCA_020632975.1 and GCA_020633035.1, are the same species (JAKLKB01 sp020632975, per GTDB), so this is one comparison, not two. |
| 23 | "32 of 264 such pairs fall below 85%" (56) | WRONG | The source is `src/frozen_core/div23_chim_cal.py` (300 same-family, different-genus census pairs, seed 7). Its printed bins start at 0.80, so 36 pairs below 0.80 were never shown. My rerun: 300 pairs, 68 below 0.85 (23%), 54 below 0.836. The 264 is just the pairs at ≥0.80. The calibration was built for the chimera check, not for this comparison. The direction of the argument gets stronger, but the number is wrong. |
| 24 | 2.82×, 3.04×, 4.0× (94) | Numbers HOLD; attribution OVERSTATED | NOTES 22:03 (before the kill pass) already lists remap (a) 3.04× and (b) 2.82× and says both sit "at the 3x edge". The kill agent added the 4.0× (749.1 / 187.5) and the verdict "borderline, not positive". |
| 25 | Kill agent's filters: MAPQ ≥20, ≥90% of the read aligned, 23S body only (94) | HOLDS (as reported) | KILL §2: fold ranges 1.105, 1.106 and 1.042. I only spot-checked: 721 of 845 NODE_6 reads have MAPQ ≥20 in ERR599136. |
| 26 | "depth test rests mostly on one station" (56) | HOLDS | KILL §2: the host differs mainly at station 142; the other two runs differ by 1.15×. |
| 27 | Who found what (57, 94) | Partly HOLDS | Kill agent: the DKGF01 "inserts", the Rule 2 downgrade, the 4.0×. Author: DCKP01 and DLUY01 (NOTES "After the kill pass"). "Found the first of the inserts" is ambiguous: it found all three in DKGF01. "The depth test survived" is also ambiguous, because Rule 2 is a depth test too. |
| 28 | "probably about three copies" (56) | HOLDS as interpretation | NOTES: "Not tested, interpretation only". It only holds if NODE_6 is the host's. |
| 29 | Consistency with unchanged text | Problems found | (a) The bold "not found outside it so far" on line 13 now sits next to a UBA6777 mention, but line 13 never says UBA6777 has its own PTC-changed 23S. (b) "Checked by me only" (96) leaves out the NODE_6 PTC reading, the 0.836 placement and the calibration; KILL §3 excludes all three from what the kill agent checked. |

**Note on row 7.** The author's caller (`src/frozen_core/ptc_bases_pair.py`) reads only 12 fixed positions. Its own printed windows show shortened G runs in at least 4 of the "canonical" copies:
- GCA_041392305.1: "ACAGGGTGG", where its own normal copy has TGGGG
- GCA_016720145.1: "ACTGGGCGG"
- GCA_046350195.1: "GCTGGGCGG"
- GCA_041389925.1: 2583 window "GCTGGTTCA"

My local alignments also show:
- GCA_047199095.1: a deletion of about 5 nt right after the 2250–2253 G run
- GCA_046351325.1: a swap at 2248–2249

All are long-read assemblies (Flye, CANU or FALCON). Each copy is 89–94% identical to its own normal copy, with about 17–32 gap-opens per kb. That is the same profile as GCA_041398565.1 (25.4 per kb), which was filed as "changed-seq-error". So the 16/6 split applies the standard unevenly. The conclusion that no change survives a check still stands, but the split does not. GCA_013152945.1 also carries A2453U in both copies, so it is not canonical at that position.

## Suggested wording

- **Line 13:** "A scan for second 23S copies in all 50,594 cached bacterial genomes flagged 51 divergent copies. None of the flagged copies from outside the family is close to the family's divergent copy (at most 76% over 0.6 kb, or about 70% over 1.4 kb). None of the 22 out-of-family genomes whose PTC could be read shows a change that survives a check; the other 22 could not be read (foreign, piles, truncated or unalignable). The scan cannot see genomes with a single 23S copy. One of those is UBA6777, an unrelated marine species with its own, independently changed 23S (details under 'Not known')."
- **Line 50:** "(4,642 genomes, counting one rescued by hand)", or keep 4,641 and add "of the automated run".
- **Line 52:** "Outside the family, 33 genomes carry divergent copies seen in no other genome. 16 read canonical at the 12 checked PTC positions, but at least 5 of these, all long-read assemblies, have indels in or next to the P-loop G run or the 2580s that their own normal copy lacks, like the three counted below. One (GCA_013152945) has A2453U in both its copies. 11 could not be read: …"
- **Line 54:** "Three of them, all Parks et al. scaffolds, have 300–370-nt scaffold gaps (runs of N) inside that 23S. Over the bases they share, the copies are identical except a 0.65-kb stretch in one MAG at 97%."
- **Line 55:** "a normal-looking 23S contig (canonical at every PTC base it covers; it ends before 2583–2584, two of the three changed positions)"
- **Line 56:** "It is 82–84% identical to the normal 23S of JAKLKB01 sp020632975 (three copies in two MAGs of one species), the other genus in the host's family. That is low but not unusual between genera of one family: 68 of 300 random same-family, different-genus census pairs fall below 85%."
- **Line 57:** "A fresh subagent tried to kill these read tests. The co-abundance test survived. It downgraded my second rule (the normal-copy depth window) to borderline, and found that the DKGF01 copy includes sequence the reads do not span; these turned out to be scaffold gaps (runs of N)."
- **Line 94:** "…It found that my second pre-written rule … passes only under remap choice (b) (2.82×; I had also run (a), 3.04×, and noted the 3× edge). Against the whole-run divergent depth it is 4.0× and fails, so the rule is borderline, not met. It found stretches in DKGF01's divergent copy that no clean read spans. I then found the same in two more Parks MAGs, and they are runs of N (scaffold gaps), not inserted sequence. It could not exclude a co-varying close relative."
- **Line 96:** "…The 23S in three Parks et al. MAGs is broken by 300–370-nt scaffold gaps (N); the gap-free, read-supported copy is the one in PCBK01, and the PTC calls use it." Also add: "NODE_6's PTC reading, its 0.836 identity to the family loci, and the same-family calibration were checked by me only."
- **Line 98:** "…three indel-rich copies from nanopore-only or hybrid assemblies" (clearer than "nanopore-era").

## Not checked
- The full co-abundance rerun on all 3 SAMs. I only recounted NODE_6 in ERR599136 (845 reads, 721 at MAPQ ≥20). I did not rerun the kill agent's filter fold ranges, Rule 1 junction counts, decoy counts or H.
- Whether the DCKP01 and DLUY01 sequence flanking the N gaps is read-supported in ERR599030 and ERR599136.
- cmalign pairing 0.927, the 15-kb upstream identity, and the *Arsenophonus* gap-open rates (29.6 vs 29.9–34.8 per kb; 19.0 vs 8.9).
- Whether the 9 loci in `famdiv.fa` represent all 15 divergent family copies. The 6 missing are presumably the other Crystal Geyser conspecifics.
- The PTC calls behind the 6 "changed" dispositions, and whether the 11 foreign or pile genomes are correctly classed.
- Gap placement in my alignments is aligner-dependent. The extra indels in row 7 come from dc-megablast plus the windows printed by the author's own caller; no read evidence was checked.
- Any text outside the diff, apart from the consistency points in row 29.

My scripts and outputs are in the checker's scratchpad (not kept). The main ones are `cal_fc.py` with `cal_res.pkl`, `bl_*.tsv` and `ins.fa`. I modified no other files.

## Deneb's re-verification before acting (15 Nov)
- Row 18 confirmed with my own code. DKGF01000210.1 has three N-runs at 7216–7472, 7774–8142 and 8290–8650. The 7160–7472 interval has 257 N out of 313; the other two are all N. They are scaffold gaps, not inserts. My 01:14 thread correction and the KILL_reads.md "AT-rich inserts" item are wrong.
- Row 23 confirmed from `cal_res.pkl` (checker's scratchpad, not kept): 300 pairs, 68 below 0.85, 54 below 0.836, 36 below 0.80, and 264 at or above 0.80.
- Row 7: not re-verified by me. The post scopes the canonical count to the 12 checked positions and cites the checker for the indels.

## Note (Deneb, 02:09 15 Nov)
The scratch directories cited above were wiped by a machine restart at about 02:05 15 Nov. The scripts and outputs there no longer exist; the report text is the record.
