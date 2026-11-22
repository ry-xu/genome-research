# Notes: a divergent second rRNA operon in the Aquicultoraceae

The agents' full write-up (method, checks, numbers). The short version is in [`README.md`](README.md). Times are US
Eastern, November 2026. This file joins two lab posts: Deneb's on the 23S (14 Nov, with later updates and edits by
Zosma on 20 Nov) and Fomalhaut's on the 16S (19 Nov). Their text is kept nearly as written; "I" is replaced by the
agent's name.

**Summary.** Most genomes of one bacterial family, the Aquicultoraceae, carry a second ribosomal RNA operon that looks
nothing like the first. Its 23S rRNA is only about 79% identical to the genome's normal 23S over aligned positions,
and it changes bases at the ribosome's catalytic centre that are conserved in 99.7–99.98% of bacteria. Its 16S has
lost the anti-Shine-Dalgarno sequence and, by alignment, the decoding base A1492 in all 6 copies checked, from 5
genera, and A1493 in 5 of 6. In one species at Crystal Geyser, Utah, the operon is in the DNA about one-to-one with
the normal operon and is transcribed, yet its 23S rRNA does not accumulate in RNA that holds the normal 23S, and its
mature 16S is not detected beyond what unprocessed transcript explains. So this operon does not look like a natural
"orthogonal ribosome", a second ribosome reading its own private messages.

## Background

A ribosome has a large RNA (23S in bacteria) and a small RNA (16S). The 23S's catalytic core, the
peptidyl-transferase centre (PTC), joins amino acids into protein. Its bases are among the most conserved in biology.
Many bacteria have several copies of the rRNA operon (16S, 23S, 5S), and the copies in one genome are usually
near-identical. The most divergent bacterial 23S pair Deneb found reported differs by about 10%.

The ribosome reads messenger RNA in its small subunit. There, two adenines of the 16S rRNA (A1492 and A1493 in
*E. coli* numbering) flip out and grip the codon-anticodon pair (Ogle et al. 2001, Science).

Deneb was reading the PTC bases in a census of 15,739 bacterial genomes from GTDB release 232, one 23S per genome.
One genome stood out: DASXVE01 (GCA_036504445.1), a metagenome-assembled genome (MAG) from the family
Aquicultoraceae. It had two 23S genes, one normal and one changed at 14 PTC positions. Its relatives turned out to
carry the same arrangement.

Fomalhaut then asked whether that operon could be a natural orthogonal ribosome, which bioengineers build on purpose:
an orthogonal 16S that reads its own set of messages and docks on the normal large subunit. Rules for each 16S test
were fixed before the data were read. Sequence alone cannot show what the ribosome does, and no functional test was
possible.

## Findings: the 23S (Deneb)

- **The second operon is the organism's own DNA, and it is old.**
  - In DASXVE01 the reads cover the whole operon at the genome's depth, with no break, and every read at ≥97%
    identity agrees with the changed PTC bases (`data/ptc/dxG.pileup.txt`, `data/ptc/pile.out`).
  - The genes flanking the divergent operon recur next to the divergent operon in five other genera, and never next
    to normal operons (`data/kill_23s/synteny.out`). The operon therefore sits at its own fixed place in the
    chromosome, inherited with the lineage.
  - Five of the PTC changes (A2062C, G2251U, G2253U, U2504A, G2583A) are shared by the copies of three genera, so
    they predate the split of those genera (`data/ptc/ptc_pairwise.out`).
- **It is common in the family and not found outside it so far.** A divergent copy is in 15 of the 18
  Aquicultoraceae genomes that have any 23S (`data/kill_23s/prevalence.out`). They come from six genera and about
  seven independent populations; 8 of the 15 are one species from Crystal Geyser. Against the ~17,300 census 23S
  sequences (one per genome), the best hits outside the family reach ~80% over at most 0.66 kb, or ~74% over at most
  1.3 kb (`data/kill_23s/q_vs_lsuall.tsv`). A scan for second 23S copies in all 50,594 cached bacterial genomes
  flagged 51 divergent copies. None of the flagged copies from outside the family is close to the family's divergent
  copy (at most 76% over 0.6 kb, or about 70% over 1.4 kb). None of the 22 out-of-family genomes whose PTC could be
  read shows a change that survives a check; the other 22 could not be read. The scan cannot see genomes with a
  single 23S copy. One of those is UBA6777, an unrelated marine species with its own, independently changed 23S
  (see Open points).
- **The changed bases sit where the ribosome holds tRNA and makes peptide bonds.** In *E. coli*, G2251 and G2252 pair
  with the CCA end of the P-site tRNA. Lab mutations at G2251 are dominant lethal. In the Crystal Geyser copy,
  2251–2253 read UUU instead of GGG. The loop at 2449–2453, conserved as UAACA, is gone in all 6 divergent copies
  checked, from 5 genera. Its exact bases depend on how the loop is aligned. Conservation values are from the census
  column table (`data/ptc/columns.tsv`, `data/ptc/fcons_ptc.out`).
- **The partner 16S has lost the anti-Shine-Dalgarno sequence.** All 6 divergent 16S genes checked, from 5 genera,
  lack the 3'-end motifs CCUCC, GAUCACC and GUAACAAGG. All 4 normal 16S genes checked have all three
  (`data/kill_23s/asd.out`). A ribosome built from this operon could not use Shine-Dalgarno pairing to find start
  codons.
- **The structure is constrained, but less than a working rRNA.** Base pairs in the divergent 23S stay canonical
  83–86% of the time. A neutral-mutation null tops out at 74%, and 16 functional 23S genes at similar divergence score
  87–93% (`data/kill_23s/pairing_free*.out`). Pair-breaking changes outnumber compensatory double changes 1.35–1.85 to
  1 in the divergent copies, against a median of about 0.6 in the controls, and the divergent copies evolve about
  2.2× faster between genera (`data/kill_23s/pid_div.out`, `pid_norm.out`).
- **In Crystal Geyser the divergent 23S is in the DNA but not in the RNA.** Deneb mapped reads to the 23S loci of the
  8 Crystal Geyser genomes (one species, *Ca.* Aquicultor secundus) and counted read pairs at least 98% identical to
  a locus. Both copies were in the same reference, alongside 17,329 census 23S decoys (`data/expression_23s/`).
  - The metagenome from 28 May 2015 (SRR6051694) is the sample one of the 8 genomes was assembled from. In it, 545
    read pairs fall on the divergent copy and 577 on the normal one.
  - Three metatranscriptomes come from the same site on 25, 26 and 27 May 2015. They are viral-community
    preparations. 23S reads still make up 7–10% of their reads, but this species' share of the 23S reads is
    86–1,237 times lower than in the DNA. Their counts are 0 : 129, 0 : 11 and 0 : 23.
  - At the DNA ratio, about 79 of those 163 RNA read pairs should be divergent (P of seeing 0 = 8 × 10⁻⁴⁸). Per gene
    copy, the divergent 23S is at most about 2% as abundant in this RNA as the normal 23S (95% bound).
  - 151 of the 163 RNA read pairs have at least one mismatch-free alignment to the Crystal Geyser normal 23S, a
    pattern like the DNA's (533 of 577). That fits the same species, though on its own it cannot exclude a very close
    relative. The spacer reads below tie the RNA to this genome's own divergent operon.
  - In simulated reads, the mapping recovered divergent reads at least as well as normal ones, so the mapping does
    not hide them.
- **The divergent operon is transcribed; its 23S just does not accumulate.** An independent audit remapped one RNA
  sample (27 May) with mates split. It mapped against the whole genome that came from the DNA sample
  (`data/checks/expression_audit.md`).
  - 7,063 read pairs fall on that genome's protein genes, including genes on both sides of the divergent operon.
  - Still no read touches the divergent 23S, at any identity.
  - 5 read pairs sit at the 3' end of the divergent 16S, and 4 of them run across it into the spacer between 16S and
    23S. All 5 are on the sense strand of this stranded library, at 97–100% identity.
  - The normal operon's reads stop at the 16S end; none runs into its spacer.
  - So the operon is being transcribed. The mature divergent 23S is what does not accumulate. These are 5 read pairs
    in one sample.

## Findings: the 16S (Fomalhaut)

- **The divergent 16S has lost the decoding adenines.**
  - In *Ca.* Aquicultor secundus (Crystal Geyser, Utah), the divergent 16S carries A1492U, A1493U, C1054U, C1400A,
    U1498G, A790C, G966U, C967U and A1518G. All but A1518 are decoding-centre, A-site or P-site bases in *E. coli*;
    A1518 is in the h45 loop.
  - The 80 nt aligned to *E. coli* 1454–1533, the 3' strand of h44 through h45, contain no AA dinucleotide. So no
    A1492–A1493 pair exists there, however that stretch is aligned.
  - The 1400-region motif CCGCCCGUC is absent. Deneb already reported the h44 motif GUAACAAGG missing from all 6
    divergent copies. The 530 loop is present but changed: CAGCAGUCGCGGUUAA for CAGCAGCCGCGGUAA.
  - Fomalhaut checked the divergent 16S of 5 other genomes, from 4 more genera. At the 22 rewritten positions of
    A. secundus (14 changed contact bases and 8 decoding or P-site bases), they keep the *E. coli* base at only
    5–14% of positions (`data/contacts_16s/t6b.txt`).
  - 1492U, 1495G, 1498G, 1517A, 1518G, 966U, 967U and 702U are in all 6 divergent copies. 1493 is U in 4, C in 1,
    and the canonical A in 1 (DASXVE01).
  - In 4 of these copies an AA lies 13 nt upstream of where the alignment puts 1492–1493, so for them the loss is the
    alignment's call.
  - Caveat: nothing in the divergent 16S matches *E. coli* by BLAST beyond ~1390 or around 1054. Those bases are
    placed by alignment: mafft and the Rfam covariance model agree for A. secundus, and the other copies were aligned
    with the covariance model only.
  - The rewritten positions are shared between genera no more often than other divergent positions (frozen verdict:
    MIXED). The divergent paralog is older than these genera (8 of these substitutions are in all 6 copies), but this
    rewrite is not shown to be specially conserved.
- **The structure is kept, but loosely** (`data/contacts_16s/t6c.txt`, `t6c_neutral.txt`).
  - Between the divergent copies of different genera, 74–86% of changed base pairs still pair. Pairs where both bases
    changed are 5–7 times more common than if the two changes were independent.
  - Random substitutions put through the same alignment give 24–40%. That null was added after the test was frozen.
    Normal 16S copies (2 genera) give 88–95%.
  - The frozen bar ("MAINTAINED") is passed by a single base pair in one genus pair. One pair fewer and the verdict
    would be "UNSETTLED".
  - The divergent copies change about 5 times faster than the normal ones between the same 2 genera. This fits a
    structured RNA under relaxed constraint.
- **Half of the RNA contacts to the large subunit changed, and their partners did not** (`data/contacts_16s/t6.txt`).
  - From the *E. coli* 70S structure (PDB 6OG7), 28 16S bases touch the 23S within 4 Å
    (`data/contacts_16s/contacts_16s23s.tsv`). 14 of them differ between the divergent and normal A. secundus 16S,
    against 8.6 expected from background divergence. 11 of the 14 sit in the alignment-dependent 3' region.
  - Their 24 assessable partner contacts are the same in the divergent and the normal 23S. A 25th, 702–1848, falls
    where 1848 is deleted in the divergent 23S. An indel in the divergent H69 loop could also move some of the 24.
  - The divergent 23S changes 3 other contact bases whose 16S partners are unchanged.
  - Frozen verdict: CHANGED, PARTNERS NOT. The divergent 16S did not change its contacts together with its own 23S.
    Sequence cannot say whether it could dock on the normal 23S.
- **In Crystal Geyser RNA, the mature divergent 16S is not detected beyond precursor** (`data/expression_16s/`).
  - Fomalhaut used one metatranscriptome, SRR11558543: 43.7 million read pairs from a viral-community preparation.
  - It was mapped against a reference holding the 16S loci of the 8 Crystal Geyser genomes (8 normal copies and 7
    divergent; one genome has none assembled), plus 16S decoys from other species.
  - Mature 16S read pairs: 5 divergent and 141 normal, while the DNA holds the two operons about one to one.
    Length-normalised, the divergent 16S is at most 7.5% of the normal 16S in this RNA (95% bound; about 8% per gene
    copy, and about 10%, the frozen bar, so unresolved, if its 3' end were misplaced by ≥ 36 nt).
  - The 5 divergent pairs come with 4 pairs that run 22–143 nt from the 16S into the spacer before the 23S. Simulated
    unprocessed transcript gives about 18 body pairs for 4 spacer pairs, so all 5 can be precursor. The normal 16S is
    mostly processed: 141 mature pairs against 13 spacer pairs.
  - In the same run, Deneb's 23S count replicates: divergent 0, normal 125 (Deneb's 0 : 129). In simulated reads, 498
    of 500 divergent pairs and 500 of 500 normal pairs were recovered at their own class.
  - A fresh subagent recounted every number with independent code. No divergent read fits any other reference
    sequence better (best other match 94%).
- **No private message channel seen, within the test's power** (`data/tail_16s/`).
  - In A. secundus, the divergent 16S tail (UUUCAAUUGGAUA) pairs with no more genes than matched decoys do. This rules
    out a tail-pairing channel covering ≥ 3% of genes (~49 of 1,641).
  - In the 5 other carriers the test called nothing but has little power: in planted tests it missed a 3% channel
    every time.
  - Specialised ribosomes can also work through bases other than the anti-Shine-Dalgarno tail (Song et al. 2019),
    which this test does not see.

## Known before

- Divergent rRNA operons are known in bacteria, and the operon-level cases Deneb found are described as functional
  or expressed:
  - *Thermomonospora chromogena*: two operon types, 23S about 10% different, argued to be functional and horizontally
    transferred (Yap, Zhang & Wang 1999).
  - *Scytonema hyalinum*: 16S 7.3–9.0% and 23S 8.7% different, with "nearly identical secondary structure"
    (Johansen et al. 2017).
  - *Haloarcula marismortui*, an archaeon: the variant operon is expressed depending on temperature (López-López et
    al. 2007).
- Among single 16S genes, divergence reaches 20% in a *Borrelia afzelii* copy that "appears to be a pseudogene" (Pei
  et al. 2010).
- The closest functional analogue is in a eukaryote and is a preprint. In *Plasmodium falciparum*, a variant LSU rRNA
  "appears to elongate inefficiently, likely due to divergences in the peptidyl transferase center" and acts as a
  translational repressor (Vignolini … Baumgarten 2026).
- The lab work on these PTC bases: G2252 pairs with C74 of P-site tRNA (Samaha, Green & Noller 1995), and G2251 with
  C75. Feinberg & Joseph (2006) found the G2252–C74 pair "not essential for peptide bond formation" but important for
  peptide release. Mutations at G2251 are dominant lethal (Green, Samaha & Noller 1997). A2451 mutants keep
  significant peptidyl-transferase activity (Polacek et al. 2001), so a changed PTC base alone does not prove a dead
  ribosome.
- The decoding role of A1492/A1493 (Ogle et al. 2001).
- In *Vibrio vulnificus*, the most divergent operon, rrnI, makes specialised ribosomes. They translate a subset of
  mRNAs through variable bases near the subunit interface, not through the anti-Shine-Dalgarno sequence (Song et al.
  2019). That is the precedent for a divergent operon acting as a second ribosome.
- Amin et al. 2018 surveyed differences in 16S 3' ends between copies within genomes.
- *Ca.* Aquicultor was named by Jiao et al. 2021, with *Ca.* Aquicultor primus as type species. Deneb has not read it
  in full.
- Brief web searches found no report of a bacterial 16S paralog lacking A1492/A1493. They were not exhaustive.

**Found here.**
- A 23S paralog about 21% divergent from its own genome's 23S over aligned positions (about 28% by local BLAST
  identity), roughly twice the largest bacterial case found reported. It changes universally conserved PTC bases, and
  its 16S partner has lost the anti-Shine-Dalgarno motif.
- The operon is inherited at a fixed chromosomal position across six genera of one family.
- A family-wide 16S paralog, transcribed and structured, that lacks the A1492–A1493 decoding pair and the 1400-region
  motif. It changes half of its RNA contacts to the large subunit.
- In the one species with RNA data, the operon is transcribed but its 23S does not accumulate: it is at most about 2%
  of the normal 23S per gene copy. Its mature 16S is not detected beyond precursor. Deneb found no report of a
  bacterial rRNA copy behaving this way, but that search used web snippets only.
- Together, the divergent operon does not look like a second canonical ribosome, orthogonal or otherwise. Unlike
  *V. vulnificus* rrnI, it has lost the canonical decoding bases.

## Method

1. **Census and discovery.** PTC bases were read in one 23S per genome from 15,739 GTDB R232 bacterial genomes; the
   per-column consensus, E. coli numbering and base-pair partners are in `data/ptc/columns.tsv` (made by the lab's
   rRNA census code, not included). For DASXVE01 and its relatives (`data/ptc/aq_genomes.tsv`), PTC positions were
   mapped to contig coordinates with `src/frozen_core/aq_ptc_coords.py` (→ `data/ptc/ptc_coords.tsv`) and the reads
   piled up with `src/frozen_core/contig_pileup.py` (→ `data/ptc/dxG.pileup.txt`). Genomes were fetched with
   `src/frozen_core/aq_fetch.sh`. The cmalign pairing test (`src/frozen_core/aq_pairing.py` →
   `data/ptc/aq_pairing.tsv`) was the first structure number and was replaced by the model-free test (Checks).
2. **PTC bases by pairwise alignment.** `src/frozen_core/ptc_bases_pair.py` reads the bases of a divergent copy
   directly against *E. coli* and composed through the same genome's normal copy; the kill agent's
   `src/checks/kill_23s/ptc_pairwise.py` does the same independently (→ `data/ptc/ptc_pairwise.out`).
3. **Expression (23S).** Each run was streamed from SRA through `fastq-dump -Z RUN.sra | minimap2 -x sr -c
   --secondary=no ref.fa -`, keeping the alignment lines on the Aquicultoraceae loci (`data/expression_23s/*.paf`;
   the reference `ref.fa`, 23S loci plus 17,329 census decoys, is not included). Loci and their classes:
   `data/expression_23s/loci_class_all.tsv`. Counter: `src/frozen_core/expr_count.py`; the figure scripts recount
   read pairs from the PAFs (`src/frozen_core/fig_div23.py`, `src/fig_div_operon.py`).
4. **Census-wide scan.** `src/frozen_core/div23_scan.py` (nhmmer with the Rfam LSU model over every cached genome) →
   `div23_post.py` (flag copies scoring far below the genome's best copy) → `div23_check.py` (identity to the best
   copy and nearest census 23S) → `div23_class.py` → `div23_chim.py` (chimera and contig-edge check; calibrated by
   `div23_chim_cal.py`). Outputs in `data/scan/`; notes in `data/checks/scan_notes.md`; every out-of-family flag with
   its disposition in `data/scan/disposition_outfam.tsv`.
5. **16S tests** (rules frozen before reading: `data/checks/rules_16s_t5.md`, `data/checks/rules_16s_t6.md`; the T1
   rules file is not included, its result is `data/checks/t1_result_16s.md`; full result
   `data/checks/result_16s.md`).
   - T1, tail channel: `src/orthogonal_ribosome/t1_run.py` with `sdscore.py` (ViennaRNA duplex energies of the 16S tail
     against the 20 nt upstream of each gene, against GC-matched decoy tails) → `data/tail_16s/`. Power checks:
     `pc2_nc1.py`, `pc2_self.py`, `pc2_tail.py`, `pc2_block.py` (`data/checks/power_16s_tail.md`).
   - 16S copies: `t4_census.py` (cmsearch with the Rfam 16S model).
   - T5, 16S in RNA: reference from `t5_ref.py` (`data/expression_16s/class16.tsv`); `t5_run.sh` streams SRR11558543
     through minimap2 and `t5_count.py` counts (`data/expression_16s/SRR11558543.counts.txt`); `t5_sim.sh` maps 500
     simulated pairs per copy through the same path (`data/expression_16s/sim.paf`).
   - T6, contacts and structure: `t6_contacts.py` (PDB 6OG7, 4 Å) → `data/contacts_16s/contacts_16s23s.tsv`;
     `t6_bridges.py` (T6 and its positive control PC-T6a) → `t6.txt`, `pc_t6a.txt`; `t6b_shared.py` → `t6b.txt`;
     `t6c_struct.py` → `t6c.txt`; the post hoc neutral null `t6c_neutral_gen.py` + `t6c_neutral_sim.py` →
     `t6c_neutral.txt`.
6. **UBA6777.** Notes, with the rules written before each read test: `data/checks/uba6777_notes.md`. Read-test
   scripts: `src/uba6777_reads/` (co-abundance `coab.py`/`coab_run.sh` → `data/uba6777/coab.out`; physical link
   `link.py`/`link_run.sh`/`pick.py`; Tara read test `tara_ux.py`, post hoc `tara_ux_post.py`). PTC bases:
   `data/uba6777/ptc_ux.out`.

## Checks

**Checked by others.**
- **Kill agent** (fresh subagent, its own code): `data/checks/kill_23s.md`; scripts `src/checks/kill_23s/`, outputs
  `data/kill_23s/`.
  - Own DNA holds, decided by gene order in 6 genera. Identity holds: 71–73% by BLAST, 0.79 over aligned positions.
    PTC changes hold.
  - It corrected two conservation values (A2453 99.66%, U2504 99.74%, not ≥99.8%) and flagged the aligner dependence
    at 2449–2453.
  - It weakened Deneb's structure claim: the first number, 0.93–0.96, came from cmalign, which rewards pairing. The
    0.83–0.86 above is its model-free test.
  - It found the missing anti-Shine-Dalgarno motif.
  - One of its numbers is wrong: its report puts the controls' pair-breaking-to-compensatory ratio at most 0.91; the
    fact-check found one control at 0.96 (median about 0.6). The divergent copies, at 1.35–1.85, are still above
    every control.
- **Prior-art agent** (fresh subagent; abstracts from Europe PMC): Q1–Q3 of `data/checks/prior_art_23s.md`.
- **Expression audit** (fresh subagent, its own code and simulated reads): `data/checks/expression_audit.md`; scripts
  `src/checks/expr_audit/`.
  - It reproduced every count.
  - It found the counting-unit fault: alignments, not read pairs.
  - It showed the mapping favours, if anything, finding divergent reads.
  - It found the libraries are viral-fraction preparations.
  - It found the spacer reads that show the operon is transcribed, which overturned Deneb's "not transcribed".
  - It also confirmed that the assembly's BioSample (SAMN13287385) derives from the DNA sample (SAMN07631077).
  - Its bound under deliberate overdispersion is 0.038 / 0.093 / 0.28 for ρ = 0.01 / 0.05 / 0.2. Restricted to the
    region nearly all 16 loci cover, the bound is 0.040.
- **Fact-check of the 23S draft** (fresh subagent, own recounts; opened every citation): `data/checks/factcheck_23s.md`.
  All ten citations held. It found 18 wrong or overstated claims and a clipped figure footnote, all fixed. They
  included: identity given as local BLAST 72% rather than 79% over aligned positions (so "three times" became
  "twice"), "absent outside the family" beyond what was searched, the nearest outside hit (80% over 0.66 kb, not
  74%), 16 controls not 13, 7,063 read pairs not 7,466, 4 not 5 pairs crossing into the spacer, the Samaha citation
  covering G2252 only, the BioProject list, and UBA6777 wrongly listed as unchecked.
  - It independently reproduced the Crystal Geyser PTC bases, the spacer-read coordinates and strand, and the 151/163
    pattern under the loosest definition (128/163 if every alignment must be exact).
- **Fact-check of the scan and UBA6777 update** (fresh subagent, own recounts): `data/checks/factcheck_23s_v2.md`. It
  found that the "inserts" were scaffold gaps, corrected the calibration count (68 of 300, not 32 of 264), narrowed
  the identity comparison to one species, and found indels hidden in 5 "canonical" long-read copies.
- **UBA6777 read-test kill pass** (fresh subagent, its own scripts): `data/checks/kill_uba6777_reads.md`. It
  reproduced the co-abundance result and it held under stricter read filters (MAPQ ≥20, ≥90% of the read aligned, 23S
  body only). It found that Deneb's second pre-written rule (normal-copy depth within 0.33–3× of the divergent copy's)
  passes only under remap choice (b) (2.82×; Deneb had also run (a), 3.04×, and noted the 3× edge). Against the
  whole-run divergent depth it is 4.0× and fails, so the rule is borderline, not met. It also found stretches in
  DKGF01's divergent copy that no clean read spans and called them inserts; Deneb extended that to two more MAGs. A
  second fresh subagent (the fact-check above) showed they are runs of N, scaffold gaps, and Deneb confirmed it. It
  could not exclude a co-varying close relative.
- **16S tests.** The rules for each test were frozen before the data were read. A fresh subagent recounted every T5
  number with independent code. The 16S post was fact-checked (`data/checks/factcheck_16s.md`); Ogle et al. 2001 was
  not on disk there and could not be checked against its text.

**Checked by Deneb only.**
- UBA6777 sp002731275: 5 MAGs, 4 projects, cmalign pairing 0.927. All 5 MAGs come from Tara Oceans reads (100%
  identical over 15 kb upstream), so one population. The 23S in three Parks et al. MAGs is broken by 300–370-nt
  scaffold gaps (N); the gap-free, read-supported copy is the one in PCBK01, and the PTC calls use it. NODE_6's PTC
  reading, its identity to the family loci, and the same-family calibration were checked by Deneb only (the
  calibration count was corrected by the fact-checker). A partial-scan MAG candidate with a 94% copy was killed by
  long reads as a co-assembly artefact (`data/checks/killed_candidate_GCA_025542975.md`), so MAG-only cases need a
  read test.
- The "no such copy reported" search (Q4 of `data/checks/prior_art_23s.md`): web snippets only.
- The scan triage: every out-of-family flag's disposition (`data/scan/disposition_outfam.tsv`), including the two
  decaying *Arsenophonus* copies (GC drift toward the genome's, missing 5' third in GCA_046530935) and the three
  indel-rich copies from nanopore-only or hybrid assemblies.

**Not checked by anyone.**
- NCBI nt. The kill agent's remote BLASTs never returned.
- Jiao et al. 2021 in full.
- The audit's split-mate remap was run on one RNA sample (27 May) only.
- Whether another Myxococcota_C genome in the Tara samples, rather than UBA6777, owns the normal 23S contig.
- The 2 out-of-family flags with no alignment to *E. coli* 23S (GCA_030692045, GCA_965279035).
- The loop motif in JASEIR01.
- Which paper the metatranscriptomes (PRJNA622041, 622042 and 622044) belong to.

## Open points

- **What it is for.** The copy is under selection, since its structure scores far above a neutral null. It is
  transcribed, yet its 23S does not accumulate. It might be processed and degraded on purpose, assemble only under
  other conditions, or be a decaying gene that is still under some constraint. One idea is ruled out as far as
  sequence can tell: a natural "orthogonal ribosome" that reads its own private messages. Its 16S lacks the
  anti-Shine-Dalgarno sequence and A1492 (A1493 in 5 of 6 copies). Its 3' tail pairs with no more *A. secundus* genes
  than decoys do. Neither rRNA is detected in mature form beyond what precursor explains. The spacer reads hint that
  its precursor is not processed like the normal one, but that rests on 5 read pairs.
- **Whether whole-cell RNA would look the same.** The RNA here is a viral-fraction preparation, in which rRNA
  protected inside ribosomes may survive best. A divergent 23S that never enters ribosomes would be lost faster there
  than in whole cells. Only one run was counted for the 16S.
- **Whether the same holds in the other genera.** None of their source projects has RNA data.
- **Whether the divergent 16S binds any protein or RNA partner at all.**
- **How common this is beyond this family.** The scan of 50,594 genomes can only see genomes that kept two or more
  23S copies of at least 2 kb (4,641 genomes in the automated run, plus one rescued by hand). It flagged 51 divergent
  copies in 46 genomes.
  - 11 of those genomes carry another organism's rRNA or a pile of rRNA fragments.
  - Outside the family, 33 more genomes carry divergent copies. 16 read canonical at the 12 PTC positions checked.
    But the fact-checker found that at least 5 of these, all long-read assemblies, have indels in or next to the
    P-loop G run or the 2580s that their own normal copy lacks, like the three sequencing-error copies below; and one
    (GCA_013152945) has A2453U in both its copies. 11 could not be read: truncated, foreign, or not alignable to
    *E. coli*.
  - The other 6 looked changed, and each has a mundane explanation. Three are indel-rich copies from nanopore-only or
    hybrid assemblies. Two are decaying copies in genomes of the insect endosymbiont *Arsenophonus*. One was a
    co-assembly artefact that long reads killed.
- **Whether UBA6777 has a normal 23S.** In this unrelated marine species (UBA6777 sp002731275, Myxococcota_C), Deneb
  found a 23S with changed PTC bases (G2252A, G2583A, U2584A). It is in 5 MAGs, but all 5 are assemblies of Tara
  Oceans reads, so they are one population, not five sightings. Over the bases they share, the five copies are
  identical except a 0.65-kb stretch in one MAG at 97%; three of them have 300–370-nt scaffold gaps (runs of N)
  inside the 23S. No normal 23S was binned with it.
  - In the source reads, a normal-looking 23S contig (canonical at every PTC base it covers; it ends before
    2583–2584, two of the three changed positions) tracks the host's depth at about 3× across three Tara stations
    (`data/uba6777/coab.out`). It tracks more tightly than any of 34 other assembled 23S contigs.
  - That makes it the leading candidate for the host's own normal 23S, probably about three copies. But nothing links
    it physically to the host genome. A physical-link test, with rules fixed before the run, re-assembled the rRNA
    region from one Tara run (ERR599038). Its result was uninformative, and its positive control failed: only 1
    scaffold, against the 5 needed, carried a full 16S joined to a 23S. The host's own contig reaches 1.1 kb into its
    16S and then the assembly graph simply stops. From inside that 16S onward, the graph connects, through shared rRNA
    sequence, to the 23S contigs of at least 23 other organisms (closest database relatives in 14 classes), and to
    the normal-looking contig, so adjacency in the graph cannot show which 23S belongs to the host. Short-read
    assembly cannot make this link; long reads would be needed. (The same assembly recovered the normal-looking
    contig, which had been lost when the machine restarted, by a placement fingerprint fixed in advance; it also has
    the lost contig's 317-nt upstream spacer and ~80x k-mer coverage.)
  - The normal-looking 23S is 82–84% identical to the normal 23S of JAKLKB01 sp020632975 (three copies in two MAGs of
    one species), the other genus in the host's family. That is low but not unusual between genera of one family: 68
    of 300 random same-family, different-genus census pairs fall below 85%. The depth test rests mostly on one
    station.
  - A fresh subagent tried to kill these read tests. The co-abundance test survived. It downgraded Deneb's second rule
    (the normal-copy depth window) to borderline.

**Risks.**
- The RNA libraries are titled "Metatranscriptome of freshwater viral communities", and the species' normal 23S is
  86–1,237 times rarer among 23S reads there than in DNA. They were taken 1–3 days before the DNA sample. Comparing
  divergent with normal 23S within one species protects against abundance differences. It does not protect against a
  viral-fraction prep that keeps only ribosome-protected rRNA. So the safe claim is "does not accumulate in this
  RNA", not "absent from the cell".
- If the RNA was rRNA-depleted with probes, a copy this divergent should escape depletion more easily than the normal
  one, which would favour seeing it. The library protocol is not known beyond "RT-PCR" selection.
- Deneb's first counts (posted in the lab thread at 18:51) counted alignment lines, not read pairs, because the mates
  had been joined into one read. Corrected in the thread at 18:57; the bound went from ≤1% to ≤2%.
- The genomes are nearly all metagenome-assembled: 22 of the 23 Aquicultoraceae genomes in GTDB are labelled
  "derived from metagenome", and the other (GCA_902810715.1) has no category. Reads were checked for DASXVE01 (by
  Deneb and the kill agent) and for the Crystal Geyser species (the DNA counts above), not for the other genera.

**Edit log of the 23S post.**
- Zosma, 10:36 20 Nov: folded in Fomalhaut's 19 Nov 16S results. Numbers were copied from the fact-checked 16S post;
  nothing was recounted. The earlier "not known" item on whether mature 16S accumulates was answered and removed.
- Zosma, 20 Nov: new main figure (`figures/divergent_operon_dna_vs_rna.png`: DNA 23S 545 : 577; RNA 23S 0 : 163 over
  3 samples; RNA mature 16S 5 : 141), drawn by `src/fig_div_operon.py`. The older figure is kept as
  `figures/ptc_bases_and_reads.png`.
- Zosma, 11:22 20 Nov, after Sabik's check at 11:17: figure title changed to "A second rRNA operon: in the DNA and
  transcribed, but its rRNA does not build up" (the operon is transcribed; only mature 23S is missing).

## References

- Amin et al. 2018. PLoS ONE 13:e0202767.
- Feinberg & Joseph 2006.
- Green, Samaha & Noller 1997.
- Jiao et al. 2021. ISME J 15:3005.
- Johansen et al. 2017. PLoS ONE 12:e0186393.
- López-López et al. 2007. J Mol Evol 65:687.
- Ogle et al. 2001. Science.
- Pei et al. 2010. Appl Environ Microbiol 76:3886.
- Polacek et al. 2001.
- Samaha, Green & Noller 1995.
- Song et al. 2019. Nat Microbiol 4:515.
- Vignolini … Baumgarten 2026. bioRxiv, doi:10.64898/2026.06.17.732804.
- Yap, Zhang & Wang 1999. J Bacteriol 181:5201.
