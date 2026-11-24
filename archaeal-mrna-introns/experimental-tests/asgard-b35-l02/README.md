# A 31-nt stretch missing from an Asgard archaeon's mRNA: what would decide it (two pages, for experimental collaborators)

DRAFT by Vega, 20 Nov 2026, for Ryan to forward or not. Nothing has been sent.
- A fresh fact-check has been applied (`agents/vega/kills/wetlab_pkg_factcheck/REPORT.md`). Sabik checks it next.
- Sequences and primers are in `sequences.txt`, made by `work/archaeal-introns/h20/pkg/make_pkg.py`. The uniqueness
  check is `uniq.py` and the second unspliced-only primer is from `unspliced_primer2.py`, both in that folder.

**What is seen.** The organism is the Asgard archaeon *Ca.* Lokiarchaeum ossiferum B35 (Schleper lab enrichment
culture; closed genome CP104013.1). In its public RNA-seq, 31 nt inside a RuBisCO-family gene are missing from
90-99 % of RNA reads, depending on run and counting method. In DNA the stretch is always present.
- **RNA, main count:** 385 of 411 molecules lack the stretch, pooled from SRR27313660 and SRR32148620 (first 10 M spots
  of each).
- **RNA, recount without an aligner:** 285 of 316 and 342 of 349 reads in those two runs.
- **RNA, replicate library** SRR32148621 (same submission as SRR32148620): 0.99.
- **DNA:** 0 of 46 fragments (34 distinct molecules) lack the stretch in SRR22273278. That culture DNA was deposited
  in 2022; its BioSample collection date is 2019-04-21 and the extraction date is not stated.

The removal joins gene UYP47557 in frame to UYP47558. A 31-nt removal can do this because UYP47558 sits in the +1
frame relative to UYP47557, and the stop of UYP47557 lies inside the stretch. The result is one uninterrupted protein,
the form that about 66 of 107 relatives encode without a break. It looks like an mRNA intron, but every RNA-seq
library used a reverse transcriptase (RT).

Public data argue against (but cannot exclude) these explanations:
- Skipping at the PCR stage: duplicate families show 0-1 mixed against 8-13 expected, depending on how families are keyed.
- A general polymerase skip at DNA folds: 0 of the 3,735 loci that fold at least as stably as L02 lose ≥ 10 %.
- A general RT skip at strong RNA hairpins: no locus with ΔG ≤ −7 kcal/mol loses more than about 1 %.
- Aligner artefacts: ruled out by direct k-mer counts.
- RNA from another enrichment member or strain: the RNA reads are 99.6–99.9 % identical to B35's closed genome, though
  the RNA and DNA come from different samples.

Two explanations other than splicing remain:
- **H2:** an RT skip caused by something unique to this sequence or its RNA, such as a modification or bound protein.
- **H3:** a B35 line that lost the 31 nt after the DNA sample and took over the culture.

The alternative to both is **H1: splicing in the cell.**

**Sequences** (mRNA sense = genome plus strand):
- Removed: `AATTAAGAATAGGATTTCCCTATCTTAGATT`. Its edges share a 2-nt direct repeat (AA), so the cut point is ambiguous
  by 2 nt; the spliced sequence is not.
- Spliced junction: `...gctcaataacatgatc|aaaccgtgtactgggc...`.
- Unspliced gene: ends ...EKRPLLNNMIN-stop.
- Spliced protein: ...RPLLNNMI|KPCTGHTAAVAGDLAYQVALGGVDVVK...

**Experiments and expected outcomes.**
- The main L02 pair gives 189 nt unspliced and 158 nt spliced.
- Read the ratio by amplicon sequencing; short-amplicon bias favours 158, so calibrate it with the mixes below.
- RT priming is gene-specific, because archaeal mRNA has no poly(A).

| Experiment (same nucleic-acid prep throughout) | H1 spliced in the cell | H2 RT skip at L02 | H3 deletion line |
|---|---|---|---|
| 1. RT-PCR with two unrelated RTs (e.g. an MMLV-type such as SuperScript IV, and a group II intron RT such as TGIRT-III or MarathonRT), flanking primers, amplicon-seq | 158 dominant with both RTs | 158 with one or both RTs. A difference between RTs points to H2; equal shares do not exclude it | 158 dominant with both |
| 2. The same two RTs on in vitro-transcribed **unspliced** RNA (401-nt template in `sequences.txt`, T7 promoter added; DNase-treated), spiked into carrier RNA | 189 only (158 < 1 %) | 158 at a high share, like the libraries (if sequence-driven; a skip that needs cellular modifications or bound proteins gives 189 here) | 189 only |
| 3. No-RT control on each RNA | no product if DNase is complete; any residual-DNA product is 189 | as H1 | residual-DNA product is **158** |
| 4. PCR on **DNA from the same culture batch** as the RNA, ideally un-DNased total nucleic acid from the same extraction | 189 only | 189 only | 158 at a share near the RNA's (~0.9) unless expression differs |
| 5. RT-free RNA: an inside-the-31-nt RNase protection or northern probe, or nanopore direct RNA (the sequenced strand is RNA even when the kit's optional RT step is used) | spliced dominates; inside-probe weak on mRNA | unspliced dominates; inside-probe strong | spliced; only 4 tells H3 from H1 |
| 6. The excised 31-nt piece (small-RNA northern with ends defined by splinted ligation; RNase R, then sequence any circle junction) | may be found; its absence does not argue against H1 | absent | absent |
| 7. Targeted MS (PRM, heavy standards): spliced junction peptide, and the unspliced C-terminal peptide | junction peptide if the fused protein is made; its absence does not refute H1 | unspliced peptide only | junction peptide |

**Notes on the experiments.**
- Experiments 2 and 4 are the cheapest decisive pair: 2 tests a sequence-driven H2 and 4 tests H3.
- 2 caveat: IVT RNA lacks cellular modifications and bound proteins. An RT skip that depends on them would still give
  189 there, so a clean 2 leaves that form of H2 open. Experiment 5 closes it.
- 5 feasibility: one 1.3-kb mRNA of a slow-growing minority organism will give very few direct-RNA reads without rRNA
  depletion. RNase protection is the more realistic option.
- 6: BHB/RtcB-type splicing usually leaves a circle. The culture's own 30-nt tRNA-SeC intron shows no circle in
  RNA-seq, so even the positive control may be negative here.
- 6: a ~31-nt band can also be a degradation fragment of the precursor, which is why its ends must be defined.
- 6: RT template switching makes false circle junctions.
- 7 peptides: the tryptic junction peptide is the 35-mer above, because trypsin does not cut K-P. Lys-C does cut K-P
  and gives RPLLNNMIK, which only the spliced form can make (the unspliced form ends RPLLNNMIN).

**Controls a wet lab will want.**
- **IVT spliced:unspliced mixes** (0, 1, 10, 50, 90 and 100 %) taken through the whole RT-PCR and amplicon-seq. They
  calibrate length bias, heteroduplex bands and template switching. Use the 370-nt spliced template for these.
- **Form-specific primer checks** on each IVT RNA.
- **RNA and reaction controls:** a reference-gene RT-PCR for RNA integrity, NTCs, and two culture batches.
- **Size standards for 189 and 158:** B35 genomic DNA and a synthetic deleted-allele fragment.
- **Two RT temperatures.**
- **Second case: S859.** A 46-nt RNA-only removal that joins two hypothetical genes, with the same unknown status as
  L02. Expect 167 nt if spliced.
- **Not-removed control: L01.** A 51-nt stretch in an HdrA-like gene that RNA does not remove (0 of 24). Expect the
  210-nt product only.
- **tRNA-SeC 30-nt intron.** Removed in RNA, so it is the splicing positive control for 5.
  - The left primer sits in the 5′ leader, so RT-PCR with this pair sees only 5′-unprocessed precursors. That makes it
    useless as a positive control for experiment 1.
  - The right primer's 3′-terminal 10 nt sit on the mature tRNA's 3′ end, not in the trailer.
  - Possible fix: all 11 spliced reads reach the trailer and none reach the leader. A left primer in the mature 5′ half
    might therefore make this pair an experiment 1 positive control. It has not been designed or checked.

**Primers** (`sequences.txt`).
- **Flanking pairs:** primer3 designs at Tm ~60 °C. Every primer's full sequence and 3′-terminal 15 nt occur once in
  the B35 genome.
  - They have not been checked against the enrichment partners' genomes; do this before ordering.
  - Coordinates in `sequences.txt` are B35 genome positions, 0-based and half-open.
  - The RT primer is the right primer for L02 and tRNA-SeC, and the left primer for S859 and L01, which lie on the
    minus strand.
- **Spliced-only forward** (`TTGCTCAATAACATGATCAAACCG`, Tm 58.4 °C): absent from the genome. With L02 R0 it gives
  98 nt.
  - Because of the AA repeat, only its 3′-terminal 4 nt mismatch the unspliced template. Show discrimination on the
    IVT RNAs.
- **Unspliced-only forward:** with L02 R0, either primer can be paired with the spliced-only primer for a qPCR ratio.
  - `CCTTTGCTCAATAACATGATCAATTAAGA` (preferred; Tm 59.1 °C, 132 nt) spans the 5′ edge, and its last 6 nt mismatch the
    spliced form.
  - `AGAATAGGATTTCCCTATCTTAGA` (Tm 53.8 °C, 106 nt) lies wholly inside the stretch.

**Material.** The culture is the Schleper lab's (University of Vienna; Rodrigues-Oliveira et al. 2023). Experiments
1-4 need fresh RNA plus DNA from one batch, or RNA from the existing batches together with their DNA. Who to approach
is Ryan's call.

**Optional, no lab needed:** `REQUEST_logan.md` lists three Logan Search queries Ryan can run in a browser to look for
the joined form in other public runs.

**Where the evidence is.**
- `outbox/2026-11-14-asgard-protein-intron/post.md`: the full draft, with risks and prior art.
- `work/archaeal-introns/confirm/RESULT.md`: the artefact tests and the independent kill attempt.
- `work/archaeal-introns/h20/RESULT.md`: why no other cultured Asgard's RNA can test L02.
