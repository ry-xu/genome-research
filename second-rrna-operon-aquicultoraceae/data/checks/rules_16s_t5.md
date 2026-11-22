<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# T5 rules: does the mature divergent 16S accumulate in RNA? (Fomalhaut, frozen before any 16S read count is made)

**Why.** Deneb (the 23S post (see notes.md)) found the divergent 23S of *Ca.* Aquicultor secundus absent from
Crystal Geyser RNA (0 : 129 read pairs in the best run; at most ~2% of the normal 23S per copy), while the divergent operon is
transcribed (5 read pairs through the divergent 16S 3' end into the spacer). "Whether mature divergent 16S accumulates" was
left open because the reference held no 16S decoys. If the divergent 16S accumulates while its 23S does not, a divergent
small subunit could join normal large subunits: the layout of engineered orthogonal ribosomes (orthogonal 16S + wild-type
23S). If it is absent too, the operon is transcribed and both products are lost.

**Run.** SRR11558543 (Crystal Geyser metatranscriptome, PRJNA622042, 43.7 M spots, 2 x 151), chosen before any 16S count
because it has the most normal-23S reads of the three runs (Deneb's PAFs: 267 alignments vs 45 and 22). One run only (download
cap). Downloaded once with prefetch, streamed through `fastq-dump --split-spot --stdout`, no FASTQ kept.

**Reference (competitive).**
- Deneb's `expr/ref.fa` unchanged (17,380 23S targets: Crystal Geyser and family 23S loci plus census 23S decoys).
- 16S loci, each with 150 nt of flank on both sides: every RF00177 cmsearch hit (--hmmonly, E <= 1e-10, >= 1,000 model
  positions) in the 65 Aquicultoria genomes (`out/aq65.txt`).
- 16S decoys: one 16S per multi-copy census species (the first copy classed `same` in out/t4/census_v2.jsonl), no flank.
- Target classes: in the 8 Crystal Geyser *Ca.* A. secundus genomes (Deneb's list: GCA_009993245.1, GCF_001871795.1,
  GCF_002779175.1, GCF_002784845.1, GCF_002791875.1, GCA_002783345.1, GCA_002782845.1, GCA_002789655.1), a 16S locus is
  **div** if its 16S is >= 97% identical (minimap2 -c, matches / alignment columns, >= 1,000 nt aligned) to the divergent
  16S of GCA_009993245.1 (JAACTV010000069), **norm** if >= 97% to that genome's normal 16S; otherwise other. Other
  Aquicultoria loci are decoys.

**Mapping and counting.** `minimap2 -x sr -t 1 -c --secondary=no` (as Deneb). A read counts if primary, aligned >= 60 nt,
identity (PAF col 10 / col 11) >= 0.98. A read pair is counted once, at the locus of its higher-AS mate.
- **mature-body pair:** both counted mates (or the single counted mate) lie within the 16S proper, i.e. no aligned base more
  than 10 nt beyond either 16S end of the locus.
- **spacer-crossing pair:** a counted mate reaches >= 20 nt beyond the 16S 3' end (precursor). Reported, not used in the verdict.
- Classes are summed over the 8 genomes (one species; copies are near-identical, so reads split between them at random).
- r = (div mature-body pairs / div 16S length) / (norm mature-body pairs / norm 16S length), lengths of GCA_009993245.1's copies.

**Expectation.** DNA ratio per copy ~0.94 (Deneb's 23S DNA counts 545 : 577 in SRR6051694; the 16S sits in the same operons,
one copy each in the assemblies).

**Verdict (primary, identity >= 0.98).**
- Gate G1 (power): norm mature-body pairs >= 30. Otherwise: uninformative.
- **ACCUMULATES:** div pairs >= 10 and r >= 0.25.
- **DOES NOT ACCUMULATE:** the 95% Clopper-Pearson upper bound on div / (div + norm) pairs, converted to r, is < 0.10.
- Otherwise: **LOW / UNRESOLVED**, reported with r and its 95% interval.
- Secondary (reported): the same on mismatch-free reads only (NM = 0). If the primary and secondary categories differ, the
  verdict is UNRESOLVED.

**Controls, before the real run is read.**
- PC-T5a (recovery, in silico): 500 simulated pairs (2 x 151, fragment 300 +- 50, 1% substitutions, seed 7) from each of
  the GCA_009993245.1 div and norm 16S, mapped and counted by the same code against the same reference. Bar: >= 80% of each
  set counted as mature-body pairs at its own class, < 1% at the other class.
- PC-T5b (replicate): in the same real run, the 23S of the same loci counted with Deneb's rule (identity >= 0.98,
  read pairs) must give div 23S <= 2 pairs and norm 23S 90-170 pairs (Deneb: 0 : 129). A miss pauses the verdict until explained.
- Reported: mismatch-free share of norm and div pairs, and norm 16S pairs on A. secundus vs the other Aquicultoria loci.

**Caveats fixed now.** The RNA is a viral-community preparation, so rRNA protected inside ribosomes may be over-represented;
a free (non-ribosomal) divergent 16S would be lost faster than one in subunits. Accumulation would not show that divergent
small subunits translate. Close relatives at >= 98% could add reads to the norm class (it makes r smaller); the mismatch-free
secondary is there for this.

## Amendments
- **T5-A1 (18:05, reference build, before any read was mapped).** The HMM-only hit of the divergent 16S runs model
  149-1462 (both ends are the divergent parts), so the frozen locus would have put the real divergent 5' and 3' ends
  outside "the 16S proper": reads there would count as edge or as spacer-crossing (precursor). This biases against the
  divergent copy. The full CM reaches model 1533 on every divergent locus (~100 nt past the HMM end). Fix: for the 15 target
  loci (div and norm alike), the 16S span is the best full-CM cmsearch hit (no --hmmonly, E <= 1e-5) in a +-400 nt window
  around the HMM hit, and the locus is that span +-150 nt. Decoy loci unchanged.
- **PC-T5a PASS (18:07, on the T5-A1 reference).** div 498/500 pairs mature-body at div, 0 at norm; norm 500/500 at
  norm, 0 at div (bar >= 80%, < 1%). Code fix before the pass (no rule change): minimap2 appends /1 and /2 to interleaved
  same-name mates, so the counter strips them to pair mates (first sim run had counted each mate as a pair).
  Counter check: on Deneb's own SRR11558543 PAF the 23S branch gives div 0, norm 129 (Deneb's figure: 0 : 129).
- **T5-A2 (18:23, after a fresh-subagent code audit, before the real run started; no real read mapped yet).**
  1. Spacer-crossing is now >= 20 nt past the 3' end, as written (code had > 20).
  2. "Both counted mates in the 16S proper": a pair whose other counted mate lies on a different target is edge, not body.
  3. "95% Clopper-Pearson upper bound" means the upper end of the two-sided 95% interval (stated now: with div 0 it needs
     norm >= 35 to give DOES NOT ACCUMULATE).
  4. The combined rule is now in code: the line "VERDICT (T5)" is UNRESOLVED when the primary and NM=0 verdicts differ.
  5. **Precursor guard (new).** The audit showed unprocessed transcript alone yields mature-body pairs (Deneb saw 5 div
     spacer pairs). Precursor-only simulation (out/t5/prec: 2,000 pairs uniform over each GCA_009993245.1 16S +-600 nt,
     wgsim seed 11, same map/filter/count) gives body : spacer = 1114 : 244 (div, k = 4.57) and 1082 : 292 (norm, k = 3.71).
     ACCUMULATES now also needs one-sided binomial P(X >= div body | n = div body + div spacer, p0 = k/(1+k)) <= 0.01 (on
     the all-read counts); if only the guard fails, the verdict is LOW / UNRESOLVED (precursor guard). On the
     precursor-only set the guard gives p 0.52. Effect on power: needs >= 24 div body pairs even with 0 spacer pairs;
     with 5 spacer pairs ~64 div body pairs (r 0.19 at norm 300, so r >= 0.25 still binds at that depth).
  6. Reported, not in the verdict: precursor-corrected body counts (body - k x spacer) and r; r with the div 16S
     insertions (91 + 70 + 34 nt, audit) excised from the length.
  PC-T5a rerun on the final code: div 498/500, norm 500/500 at own class, 0 cross (unchanged).
- **Code note (18:27, during the real run, before its PAF was read; no rule change).** The counter now streams a pair as a run
  of consecutive same-name PAF lines (memory: ~6M lines expected). Output identical to the dict version on all 11 test
  PAFs (sims, precursor sims, Deneb's PAF, audit synthetics). Post-run check: pairs_seen = distinct read names.
- **T5 result (18:55, counts.txt read once, as frozen):** PC-T5b PASS (23S div 0, norm 125; bar <= 2, 90-170).
  pairs_seen 1,423,175 = distinct read names. Body pairs div 5, norm 141 (NM=0: 4, 128); r 0.031, 95 % CP 0.010-0.075 ->
  **DOES NOT ACCUMULATE** (both tiers). Spacer pairs div 4, norm 13; div body:spacer 5:4 is the precursor-only ratio
  (guard p 0.99; precursor-corrected div body -13). mm2 1,730 s, peak RSS 0.81 GB.
