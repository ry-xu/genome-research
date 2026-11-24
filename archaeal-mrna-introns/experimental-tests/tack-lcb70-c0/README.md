# A 22-nt stretch missing from a *Methanosuratincola* mRNA: what would decide it (for experimental collaborators)

DRAFT by Mintaka, 20 Nov 2026, for Ryan to forward or not. Nothing has been sent.
- Sequences, primers and product sizes are in `sequences.txt`, made by `work/tack-rna/pkg/make_pkg.py`. The evidence is
  in `work/tack-rna/RESULT.md`.
- Checks (paths under `work/tack-rna/`):
  - Kill checks by fresh subagents on raw reads: `review/killcheck_c0/` and `review/killcheck_a5/`.
  - Fact check of an earlier draft by a fresh subagent: `review/factcheck_pkg/REPORT.md`. All its points are fixed here.
  - Vega and Sabik checked the earlier primers and wording (thread 13:08, 13:09); their fixes are in. New primers:
    `review/primers_v2/`. All 9 inserts and their flanks: `review/flank9/flank9.out`. Survey lines updated to Sabik's 62-rep survey (thread 13:42); 25-copy fold lines added from Vega 13:59 and Sabik 14:00.
- **For Ryan:** this is a hand-off note for whoever might run the wet-lab test (the LCB70 group is the natural one). Nothing
  here needs an account. Recommendation: worth forwarding once Sabik's outbox check passes; the claim is "a candidate
  intron", not a demonstrated one.
- Counts are templates (sequenced molecules, duplicates collapsed) unless stated otherwise.

**What is seen.** *Ca.* Methanosuratincola verstraetei LCB70 is a methyl-reducing methanogen of the TACK class
Methanomethylicia, grown as an enrichment culture by the Hatzenpichler lab (Kohtz et al. 2024, Nature 632:1118; genome CP186953.1). In its public
metatranscriptomes (PRJNA916083: acetate, lactate and pyruvate incubations, 3 replicates each), a 22-nt stretch inside a
radical SAM gene is missing from most RNA molecules.
- **RNA:** 45 of 51 templates lack it (95 % CI 0.76–0.96). Per library: 1–7 joined, 0–2 unjoined; joined outnumber unjoined in
  all 9, though one library has only 1 joined and 0 unjoined.
  Joined and unjoined reads are both sense RNA.
- **DNA:** 0 of 58 templates lack it in metagenome SRR23044201. This rules out a deletion fraction above 5 % in that
  sample (one-sided 95 %).
  - That sample's collection date is 2021-07-20, which is also the date on the genome's sediment sample and may be the
    field date. The RNA samples are dated 2023-10-10.
  - The genome was assembled from this metagenome plus a nanopore run (SRR23064070, not checked). So this is the only DNA
    evidence, and it is not DNA from the RNA cultures. None of that is public.
- **A second species, from a second lab's data:** *M. petrocarbonis* (PRJNA992765) has a different 23-nt insert at the
  same site.
  - RNA: 5 of 10 templates lack it. The 5 joined come from replicates AM32 (2) and AM33 (3). The 5 unjoined come from AM31
    (4) and AM32 (1).
  - DNA: 0 of 44 templates lack it. This is a different sample (AM_stage3, no collection date), and the genome was
    assembled from it.
  - This is weak support. The counts are thin. These libraries also lose an unremarkable, non-BHB 63-nt stretch of the
    same gene in 2 of 33 templates, so their background for skips is not zero.

**The site across the genus and family.**
- Removal joins locus 06155 in frame to 06160, giving one 364-aa protein. The junction is …CACGGTAT|TTAGCC…. Without removal, a TAG stop starts
  1 nt after the insert (the annotated end of 06155).
- All 9 *Methanosuratincola* species-representative genomes carry an insert of 20–24 nt at this site, between …CACGGTAT
  and TTAGC… (`sequences.txt`). Each one puts a stop codon in the unspliced frame, including the two 24-nt inserts.
  Removing any of them restores one reading frame at the same junction. The 16 nt on each side are identical in 6 of the 9
  genomes; the other 3 differ at one base.
- Sabik's survey (`work/referee/sabik/c0/fused/SURVEY62.md`, all 62 Methanomethylicia species reps searched) found more:
  - The same element sits at the same codon in 8 more Methanomethylicaceae genera (16 genomes; list in SURVEY62.md).
    All copies end CGGCCGA[AG].
  - 5 of the 25 copies are in frame (21, 24, 24, 33 and 33 nt), and each carries a stop inside the insert.
  - Culexarchaeales (12) and Nezhaarchaeales (13 with a full-length ortholog) encode the protein in one piece.
    *Methanomethylicus* lacks the element.
  - So the joined protein is the conserved form, and the insert arose within one family.
- That pattern fits an intron. It also fits a DNA element that excises, which DNA from a different sample cannot rule out.
  - Against the best-known kind of excising DNA hairpin (IS200/IS605 and REP elements, cut by Y1-family TnpA/RAYT): each
    genome has one near-identical copy only (≤3 edits), no Y1 gene lies within 10 kb of any copy (nearest 11 kb), and
    18 of the 25 genomes, including the closed GCA_030522375.1, have no Y1 gene anywhere. Other hairpin-cutting enzymes
    are not excluded.

**Why BHB.**
- The 22-nt insert folds into a hairpin with the bulge-helix-bulge shape that the archaeal tRNA-splicing endonuclease (EndA)
  cuts. Vega's check: all 4 central-helix positions pair, with 3-nt bulges at the two cut sites, at the reported junction
  only.
- The fold holds across all 25 copies in 9 genera (Vega, `agents/vega/c0_fold/c0_homologs2.py`; recounted from the
  sequences by Sabik). The central helix is paired at the exact insert boundaries in 25 of 25, and the 3′ bulge is unpaired
  in 25 of 25. JAULTD01's 33-nt copy is the only exception for the 5′ bulge and the inner helix.
  - Two inner-helix pairs change both bases between copies and stay paired. That is compensatory change: selection keeps a
    paired inner stem. It does not say whether the hairpin acts in RNA or in single-stranded DNA (some transposases
    recognise DNA hairpins). What points to the archaeal splicing endonuclease is the BHB sitting at the exact insert ends
    in 25 of 25.
  - No other pairing of insert positions does this. Of 112 decoy position pairs, none stays paired in all 25 copies across a
    double change: 0 of the 34 whose two columns vary, and 0 of the 12 whose majority base differs between the two families
    set by the first of these two pairs (Mintaka, `agents/mintaka/checks/c0_covar/covar.out`; Sabik,
    `work/referee/sabik/c0/fused/decoy_split.out`).
  - On the GTDB r232 species tree, parsimony needs at least 4 independent double changes, 2 in each pair (one of the 4 is in
    a single genome), plus one GU wobble, and every ancestor can stay paired. Merging the first pair's two changes into one would
    break a split with 100 % support. This assumes the insert followed the species tree (`c0_covar/tree.out`).
- c0's insert folds at −5.6 kcal/mol, more tightly than 90 % of shuffles of its own bases. The *M. petrocarbonis* insert
  folds more tightly than 96 % of its shuffles. The 63-nt stretch is no tighter than its shuffles.
- The only archaeal protein-gene intron established so far is in crenarchaeal *cbf5* (Watanabe et al. 2002, FEBS Lett
  510:27; Yoshinari et al. 2006, BBRC 346:1024; Yokobori et al. 2009, BMC Evol Biol 9:198). Proposed but unconfirmed cases
  include SSO1586 (Doose 2013), *M. kandleri* RNA circles (Berkemer 2015), and this lab's L02 (*Lokiarchaeum*, Asgard).
  EndA also cuts a BHB in the 5′ UTR of the *Haloferax* HVO_1309 mRNA, outside the coding sequence (Schwarz et al. 2020,
  as cited by Liang et al. 2025, MMBR 89:e00318-24).

**Explanations the public data argue against (but cannot exclude):**
- **Aligner or census artefact:**
  - 80 of 81 joined alignments match the spliced form exactly, with 0 nt microhomology at the junction.
  - None of 14 Methanomethylicales genomes carries the joined form contiguously.
- **A general reverse-transcriptase skip** (*M. verstraetei* libraries):
  - Other ≥10-nt gaps across all index windows come to 23 alignments (17 molecules) at 6 places, at most 11 at one place.
    c0 has 81 alignments at one exact junction.
  - Among 22-nt stretches that fold at least as tightly as c0's insert (≤ −5.6 kcal/mol), 84 loci with 21,558 spanning
    reads give 0 removal reads.
  - c0's ±40-nt region folds at −48.7 kcal/mol. Only two background loci fold at ≤ −44. They give 0 removal reads in 944,
    909 of which come from one locus.
  - Realigned to the whole genome, c0's reads show no other skip of 15 nt or more along the gene pair: 724 positions have
    at least 5 templates, median 13.
- **A deletion strain that took over after the DNA sample:** this would need the same exact removal, of different inserts,
  in two species. The second species' counts are thin.

What remains:
- **H1:** splicing in the cell.
- **H2:** an RT skip caused by something specific to this RNA, such as a modification or a bound protein.
- **H3:** loss of the insert from the DNA (deletion lines, or an element that excises) in both cultures after their DNA
  samples.

**What would decide it**
1. **RT-PCR on DNase-treated RNA, with a no-RT control, using two unrelated RTs.** Use a retroviral RT at its highest
   recommended temperature and a thermostable group II intron RT (about 60 °C). The region folds tightly, so how each RT
   copies it is part of the test.
   - Presence: F + R give 359 bp (spliced) or 381 bp (unspliced). Resolve on a high-percentage gel or PAGE. Mixed
     templates can add heteroduplex bands.
   - Share: end-point PCR favours the shorter product. Measure the share by qPCR with J + R (spliced, 210 bp) against
     U + R (unspliced, 229 bp), or by sequencing low-cycle F + R amplicons.
   - Nearest-neighbour Tm (primer3 defaults, 50 mM Na+, 1.5 mM Mg2+): F 64.9, R 65.2, J 65.5, U 66.3 °C. The worst
     dimer is F with itself at −6.2 kcal/mol. U ends before the insert's self-complementary CGGCCG.
   - J's first 20 nt also match the unspliced template. Genomic DNA and an unspliced in vitro transcript must give no
     J + R product.
   - DNA: extract it in parallel from the same culture aliquot and amplify with F + R.
   - Predictions:
     - H1: both RTs give a similar spliced share, and the DNA gives only 381 bp.
     - H2: the share depends on the RT.
     - H3: the DNA shows the 359-bp form.
2. **Direct RNA nanopore sequencing.** No RT is involved, so this tests H2.
   - Archaeal mRNA has no poly(A) tail, so it must be tailed in vitro first.
   - The transcript is rare: 0.54 spanning templates per million RNA-seq spots. About 63 spanning reads are needed to put
     ±0.10 on a share near 0.88, so enrichment is probably needed.
3. **DNA and RNA from the same flask.** This tests H3. It needs depth: 0 deletion reads in 59 spanning templates rules out a
   deletion fraction above 5 %, and 299 rule out 1 %. Deep F + R amplicon sequencing of the DNA is the easy route.
4. **In vitro:**
   - Reverse-transcribe and sequence an unspliced in vitro transcript of the region (an RT-skip control).
   - Add recombinant archaeal EndA to test the cut directly.

**Where the material is.**
- LCB70 data: Hatzenpichler lab (Montana State University).
- *M. petrocarbonis* data: submitted by the Biogas Institute of the Ministry of Agriculture and Rural Affairs, Chengdu
  (SRA centre for PRJNA992765).
- Whether either culture is still maintained is not known from the data.
