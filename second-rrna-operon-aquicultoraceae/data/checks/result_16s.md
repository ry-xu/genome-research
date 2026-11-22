<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# Result: a natural orthogonal ribosome (Fomalhaut, 19 Nov 2026) — final (T1-T6 audited; outbox draft fact-checked)

**Verdict: the divergent Aquicultoraceae operon is not a natural orthogonal ribosome.** Its 16S shows no private
Shine-Dalgarno channel (T1), and no partner tRNA class or carrier-specific CCA enzyme (T2). In Crystal Geyser RNA, the
mature divergent 16S is absent (T5: 5 vs 141 read pairs, all explained by unprocessed transcript), as Deneb found for
its 23S, though the operon is transcribed. Family-wide, the divergent 16S keeps its base pairing (T6-C, narrowly) while evolving ~5x faster
than the normal 16S. It lacks the canonical decoding bases (A1492 in all 6 copies from 5 genera, A1493 in 5 of 6) and
differs at 14 of its 28 RNA bridge residues without matching changes in its 23S partners (T6). Among 5,199 complete GTDB species, 9 carry a
16S copy whose anti-SD core is lost (T4). Recurrence across genera rests on Alteromonas (published) and a decaying
Borreliella paralog; the follow-up channel test calls nothing but has no power at 3%.

Rules: `RULES.md` (frozen v1, md5 777457b87d25e59d9e83f59a578732cd, amendments A1-A8), `RULES_T4.md` (T4-A1..A10),
`data/checks/rules_16s_t5.md` (T5-A1, A2), `data/checks/rules_16s_t6.md` (T6-B, B1, C).
Notes per test in `notes/`. Thread: the lab discussion thread (not included) (posts 16:21, 16:59, 18:35).

## Controls
| control | result |
|---|---|
| PC1, canonical SD gradient | PASS: Bacillota median 0.861 vs Bacteroidia 0.022; 10/10 vs 0/10 above 99% of decoys |
| PC2, random planted pair as frozen | FAIL (3/20 at 3%): decoy 13-mers vary widely in pairing propensity, so random tails sit mid-null. Power read per tail instead (A3) |
| Tail-specific power (A3, A7) | A. secundus tail: 10/10 at 3%, in a non-carrier and in its own genome. DASXVE01 tail: no detection limit within 10% in its own genome |
| NC1, decoy false calls | PASS: 0/120 |
| Code audits (fresh subagents) | T1 code: clean after A6. T4 code: two audits, fixes T4-A7, T4-A9 |
| PC-T4a/b/c (census) | PASS on the final code: E. coli 7/7 uniform; planted core-less tails 20/20 core-lost, planted point changes 20/20 point; A. secundus and DASXVE01 divergent copies core-lost |

## T1. Channel test (6 carriers with an assembled divergent 16S)
No carrier is called (table in `data/checks/t1_result_16s.md`). For the A. secundus tail (UUUCAAUUGGAUA), no private channel
covers >= 3% of genes (~49 of 1,641). For DASXVE01 the test has no power. The four carriers whose divergent tail does
not align were scored as a max over 58 windows. Their best window is the same G-rich element (GGUGUGGGC) with no
positional enrichment, which fits the G-U wobble background. Block power (planted in non-carrier GCA_018333515.1, 5 reps per share, `data/tail_16s/pc2_block.log`): calls 0/5 at 0 % and 3 %,
1/5 at 5 %, 2-3/5 at 10 % for each of the 4 carriers (median z 5.2-8.4 at 10 %; misses fail the positional-ratio gate), so
for these 4 the block test rules out only channels well above 10 % of genes.
Scope: a small regulon (a handful of genes) is below what this test can see. A no-call does not mean the divergent
ribosome is functionally equivalent (Song et al. 2019: specialized ribosomes can act through intersubunit residues, not the ASD).

## T2. Partner tRNAs and CCA enzyme
- Second CCA-family enzyme: dies. Two PF01743 proteins occur family-wide, in non-carriers too, and in many other families.
- tRNA 3' ends: the frozen rule is uninformative (CCA encoded in 23-45% of carrier tRNAs, below its 50% bar).
  Descriptive: 99% of carrier tRNA genes encode C74; 2/639 encode the A74-A75 end a UUU P loop would pair with.

## T3. Locus (descriptive)
The divergent operon sits at a fixed, conserved locus in all six carrier genera, with no mobile-element genes within 10 kb.
In-family "non-carriers" break their contigs exactly at the operon side, so some may carry it unassembled.

## T4. Census of 16S anti-SD heterogeneity (5,199 complete GTDB species, one per species)
| class | species |
|---|---|
| uniform (all copies same tail) | 3,105 |
| single copy | 2,054 |
| point (tail differs, core kept or absent in both) | 31 |
| core-lost (majority has CCUC/CUCC, a copy has neither) | 9 (2 in 2-copy genomes, 7 in 3+) |
4 genomes had copies set aside as unresolved (contig edge or CM end before the tail).

Core-lost genomes (assembly level, reads not checked):
Mycoplasmopsis GCA_024582715.1; Borreliella turdi GCA_030439285.2 and B. lusitaniae GCA_030440365.2;
G047199115 GCA_047199115.1; Alteromonas mediterranea GCF_000020585.3 and A. macleodii GCF_000172635.2
(pre-registered recall check T4-A5: passes; known from Amin et al. 2018); Thermanaerovibrio GCF_000024905.1;
Vibrio GCF_009938165.1 (2-nt deletion in the ASD, body identical, so possibly a sequencing error); Tepidibacter_A GCF_934924865.1.

Recurrence (T4-A6, label permutation): 2 recurring genera vs 0.37 expected, p 0.015 (Borreliella p_bonf 0.003,
Alteromonas 0.029). With >= 3 resolved copies only (T4-A9): 1 genus, p 0.42. 
Borreliella (unfrozen look at 6 more genomes): the minority copy is a 79-83%-identical 16S-like paralog ~450 nt
upstream of rrs on the same strand, in B. turdi (2 strains), B. tanukii and B. afzelii (2011 assembly, another
submitter); absent from B. burgdorferi B31 and the B. garinii assembly checked. Its tails differ by species, and the
turdi/tanukii copies lack the first 100-200 model positions, which looks like a decaying paralog rather than a conserved
second ASD. Likely known: Deneb's the 23S post (see notes.md) cites a *B. afzelii* 16S copy 20% divergent
that "appears to be a pseudogene" (Pei et al. 2010, Appl Environ Microbiol 76:3886); the afzelii paralog here fits that
(citation taken from Deneb's post, paper text not read). Not a novelty claim.

Follow-up T1 on core-lost copies (T4-A8/A8b; 10 copies in 9 genomes, Bonferroni over copies, 10k-decoy top-up,
power at 3% x 10): **no call** (`out/t4/follow.jsonl`). Lowest p: B. turdi 0.021 (z 2.66), B. lusitaniae 0.098; the
other 8 copies p 0.70-1.00. **Power 0/10 for every copy**, so these no-calls say little. A rerun of two Tepidibacter
power reps (the most non-canonical tail, GAUCCGGAGAUUU) shows why: of ~90 planted genes, 48-55 pass the exclusivity gate,
and the in-genome decoy null is heavy-tailed (z 3.5-3.6 gives p 0.04 against 1,000 decoys), the PC2 problem again.
Scope: the follow-up cannot see a private channel covering ~3% of genes for any core-lost copy.

## T5. Does the mature divergent 16S accumulate in RNA? (`data/checks/rules_16s_t5.md`, frozen 18:02, amendment T5-A1)
Controls: PC-T5a PASS (simulated pairs 498/500 div, 500/500 norm at own class, 0 cross). Counter reproduces Deneb's
23S 0 : 129 on their PAF. PC-T5b PASS (23S of the same loci in the same run: div 0, norm 125 pairs; Deneb's 0 : 129).
**Verdict: DOES NOT ACCUMULATE.** SRR11558543 (43.7 M spots): mature-body pairs div 5, norm 141 (mismatch-free 4 : 128),
r 0.031 (95 % CI 0.010-0.075; DNA expectation ~0.94). The 5 divergent body pairs come with 4 spacer-crossing pairs (22-143 nt into
the spacer); unprocessed transcript alone would give ~18 body pairs for 4 spacer pairs (simulation 4.6 : 1), so the
divergent signal is fully explained by precursor (guard p 0.99) and mature divergent 16S is at or near zero. The normal
16S is mostly processed (141 : 13). Fresh-subagent audit: independent recount identical; all 18 divergent mates are
divergent-specific (best other target 94 %). Sensitivity: upper bound 0.075 vs the 0.10 bar; 0.084 per copy (one of 8
genomes has no assembled divergent 16S); 0.093 if the two pairs 22-26 nt past the 3' end are mature; 0.102 (LOW /
UNRESOLVED) if the divergent 3' end is misplaced by >= 36 nt (body rule has 10 nt slack; 0.110 at >= 133 nt). Counting the 10 norm pairs split across identical
norm copies as body lowers the bound (0.070). Like its 23S, the divergent 16S is transcribed but does not accumulate.
Caveat: viral-fraction RNA (ribosome-protected rRNA may be over-represented); one run.

## T6. Can the divergent 16S dock on the normal 23S? (`data/checks/rules_16s_t6.md`, frozen before any A. secundus base at a contact was read)
RNA-RNA intersubunit contacts from PDB 6OG7 (E. coli 70S, heavy atoms <= 4 A): 53 pairs, 28 16S residues, 32 23S residues.
Control PC-T6a (covariation logic, intra-16S pairs, A. secundus div vs norm): PASS, both-changed pairs 3.5 x independence,
89 % of them still pair (`data/contacts_16s/pc_t6a.txt`).

**Verdict: CHANGED, PARTNERS NOT.** All 28 bridge residues assessable. **14 of 28 differ between the divergent and normal
16S** (8.6 expected from the 31 % background divergence, P(>=) 0.02). Their 23S partners are the same in the divergent and
normal 23S, residue by residue (0 of 24 contact pairs); but the divergent H69 loop is 2 nt longer (1914, 1915, 1918
changed) and 1848 is deleted at B7a, so with another placement of that indel up to 4 of the 24 could change (verdict
unaffected). So the divergent 16S did not change its bridges together with its own 23S, and sequence cannot say whether
it docks on the normal 23S. The normal A. secundus 16S itself differs from E. coli at bridge residues 1429, 1476, 1484,
1485. Changed residues (same set in mafft and cmalign): 702 (h23); 783, 784 (h24); 1419, 1428, 1430, 1473-1476, 1493,
1495, 1496 (h44); 1517 (h45). 11 of the 14 lie in the 3' minor domain, where the divergent copy has no blastn HSP to
E. coli or to its normal copy past E. coli ~1390, so the bases called there depend on the alignment.

**Unfrozen look, then frozen follow-ups (T6-B, T6-C).** The divergent 16S also carries substitutions at universally
conserved decoding and P-site bases (cmalign): A1492U, A1493U, C1054U, C1400A, U1498G, A790C (790 loop rewritten), G966U,
C967U, A1518G. The 80 nt aligned to E. coli 1454-1533 (3' strand of h44 through h45) contain no AA dinucleotide, so an
A1492-A1493 pair is absent however that stretch is aligned (the 5' strand of h44 has 6 AA; DASXVE01 likewise has none in
1447-1530, but the other 4 copies have an AA aligned to 1478-1479, so for them the loss is cmalign's call). The
h44 motif GUCGUAACAAGG, the 1400-region CCGCCCGUC and the 530-loop CAGCAGCCGCGGUAA are absent (530 loop present as
CAGCAGUCGCGGUUAA).
- T6-B (shared across genera?): the divergent copies of 5 other carriers (4 more genera) carry the E. coli base at only
  5-14 % of these 22 positions. 1492U, 1495G, 1498G, 1517A, 1518G, 966U, 967U, 702U, 783G and 784C (cmalign; mafft places
  783/784 differently) are in all 6 divergent copies. 1493 is U in 4, C in DASQQM01, and the canonical A in DASXVE01; 1054
  is the canonical C in both CADDYS01 copies; 790 is a gap in 3. Frozen verdict **MIXED**. The rewritten positions are
  shared no more than other divergent positions (all 5 others match A. secundus at 10/20 R vs 201/399 G positions). The
  divergent paralog as a whole predates the split of the carrier genera; the specific rewrite is not shown to be conserved.
- T6-C (pairing kept since the split?): **MAINTAINED, by the narrowest margin** (Kmv33 keeps 66/82 = 0.805; one pair fewer
  and the verdict would be UNSETTLED). Between divergent copies of different genera, 0.74-0.86 of changed pairs still pair
  and both-changed pairs are 5-7 x independence (normal copies 0.88-0.95). Post hoc random-substitution nulls through
  the same cmalign path give 0.24-0.40 (uniform) and 0.32-0.40 (55 % transitions), 1.3-2.2 x. But the divergent copies
  change about 5 x faster than the normal ones (per-side change 0.13-0.16 vs 0.03 between the same genera).

**Reading.** Family-wide, the divergent 16S keeps its base pairing far above a random-substitution null while evolving
about 5 x faster than the normal 16S, which fits a structured RNA under relaxed constraint. It is transcribed but not
accumulated (Deneb, T5). It lacks A1492 in all 6 copies and A1493 in 5 of 6, and differs at 14 of its 28 RNA bridge
residues. With Deneb's PTC changes in the divergent 23S, this operon does not look like a second canonical translating
ribosome. A "natural orthogonal ribosome" in the engineered sense (orthogonal ASD, standard decoding) is not what it is.
Prior art: none found in two web searches for 16S paralogs lacking A1492/A1493 or for Aquicultoraceae rRNA (not exhaustive).
Caveats: one E. coli state; protein bridges not tested; one bridge contact (702-1848, B7a) dropped as 1848 is deleted in
the divergent 23S, and the divergent 23S changes 3 contact bases (1915, 1950, 1960) whose 16S partners are unchanged.
Fresh-subagent audit: contacts, numbering and counts reproduce;
wording here follows its corrections.
Files: `data/contacts_16s/t6.txt` (mafft), `t6b.txt` (cmalign), `t6c.txt`, `t6c_neutral.txt`, `pc_t6a.txt`; code `src/orthogonal_ribosome/t6_contacts.py`,
`t6_bridges.py`, `t6b_shared.py`, `t6c_struct.py`, `t6c_neutral_gen.py`, `t6c_neutral_sim.py`.
