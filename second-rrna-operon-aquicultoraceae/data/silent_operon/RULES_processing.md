<!-- Copied from the lab record (work/silent-operon/, 22 Nov 2026). Paths such as out/, src/, recheck/ or frozen-core/ refer to the original working folder; the small tables from out/ are in this folder (see ../README.md), the rest is not included. -->
# A_processing rules (written 22 Nov 2026, before any result was looked at)

Question: can the divergent operon mature, and how is it expressed?
Genomes: the 15 Aquicultoraceae carriers (`orthogonal-ribosome/out/aq_status.tsv`), fasta from
`frozen-core/results/lsu/bac/cca/aq/fna/`. Controls: the normal operon of the same genomes; *E. coli* K-12 (7 operons,
GCF_000005845.2); *Thermobispora bispora* DSM 43833 (GCF_000092645.1), whose two operon types differ ~6% in 16S and are
both transcriptionally active (Wang, Zhang & Ramanan 1997, J Bacteriol 179:3270) = a known divergent working operon.

Locating: 16S by cmsearch RF00177, 23S by RF02541, 5S by RF00001, tRNAs by ARAGORN 1.2.41. Mature ends = the genomic
positions aligned to the first and last model positions of a hit that is complete at that end. An operon is "div" or
"norm" by which locus class (frozen-core `loci_class.tsv` / t0 DIV list) or, failing that, by best blastn identity to the
DASXVE01 div vs norm 16S/23S. Only operons whose relevant flank is assembled (>=200 nt beyond the mature end, on the
same contig) enter a test; the rest are reported as "not assessable".

## T1 RNase III processing stems
Leader L = 200 nt upstream of the mature 5' end; trailer T = 200 nt downstream of the mature 3' end (16S and 23S
separately). Stem = RNAduplex(L, T) (ViennaRNA 2.8, 37 C). Recorded: dG, number of base pairs, longest helix (pairs with
at most 1-nt bulges/1x1 loops), and the distance of the helix from each mature end.
Null: 200 pairs of dinucleotide-shuffled L and T (each shuffled separately); p = fraction of nulls with dG <= observed.
- **Stem forms** = p <= 0.05 AND >= 15 bp AND the helix lies within 150 nt of both mature ends.
- Calibration: the rule must call stems in >= 6/7 E. coli 16S and 23S and in both T. bispora operon types, else T1 is
  "uncalibrated" and reported descriptively only.
- Calls:
  - div stems absent in most assessable div operons while present in the same genomes' normal operons -> "div rRNA
    cannot be excised by RNase III; maturation block explains non-accumulation" (a lead; recheck by fresh agent).
  - div stems present at rates like the normal -> "excision stems intact; non-accumulation is downstream of RNase III
    (later maturation, assembly or turnover)".
  - mixed (16S stem present, 23S absent or the reverse) -> reported per rRNA; a broken 23S stem with an intact 16S stem
    fits "16S made, 23S not".
- Descriptive: whether the div leader/trailer pair with each other in the same register across genera (helix positions
  relative to mature ends conserved) -> structure kept by selection.

## T2 Promoter and antitermination
Region = 600 nt upstream of the mature 16S 5' end (or to contig end). Recorded:
- Best sigma70/HrdB-type promoter: -35 TTGACA and -10 TATAAT, spacer 16-18, scored as matches (max 12); also the best
  extended -10 (TGnTATAAT) and an AT-rich UP element (>= 70% A+T in the 20 nt upstream of -35). Discriminator = the 6-8
  nt between -10 and the presumed +1, reported as GC count.
- boxA: best match to the rrn boxA consensus TGCTCTTTAACA (and the shorter core GCTCTTTAA) within the 600 nt; called
  present at <= 2 mismatches over 12 nt. Also searched in the 16S-23S spacer (E. coli has one there).
- boxB/nut hairpin: descriptive only (stem-loop 5' of boxA).
- Calls: div intact = boxA present and a promoter scoring >= the lowest-scoring normal operon's promoter in the same
  genome. "Changed" = present but with a different sequence/score than norm. "Absent" = no boxA <= 2 mm and no
  promoter >= 8/12.
- Conservation: one operon per genus (6 genera). Mean pairwise identity (global, blastn-free Needleman-Wunsch over the
  300 nt upstream of 16S) of div upstreams across genera vs norm upstreams across genera. Div upstream conserved at
  least as much as norm upstream, with the boxA/promoter in the same position -> "selection on expression". Div
  upstream at background identity (<= 55%, ~random-aligned) -> no evidence.
- Controls: E. coli (all 7 operons should show boxA <= 2 mm; promoter >= 9/12), T. bispora both types.

## T3 Intergenic content
- 16S-23S spacer: ARAGORN tRNAs (aa, anticodon, ARAGORN score, encoded CCA or not); spacer length.
- 23S downstream 400 nt: 5S by cmsearch RF00001 (complete = model 1..end; score vs the normal 5S in the same genome);
  plus any tRNA.
- Calls: "intact" = spacer tRNA with a canonical anticodon for its aa and an ARAGORN score within 10% of the normal
  spacer tRNA of the same aa (if any), and a complete 5S scoring >= 80% of the normal 5S. Missing or truncated = "absent".
- Interpretation: intact tRNAs/5S in div -> the operon is built to be expressed as a full rrn transcription unit
  (selection on the whole unit). Missing 5S -> a div 50S would need the normal 5S in trans.

Each test gets a RESULT.md section (one sentence, numbers, mechanism). One headline line per test to
`lab/commons/2026-11-22-silent-operon.md`. Light compute: 1 thread.

## Amendment A1 (19:05, after seeing control and normal-operon T1 rows only; no div T1 row had been viewed)
The frozen duplex rule failed calibration on *T. bispora* 16S (p = 0.075-0.19 in all 4 copies; 23S passed 3/3; E. coli
passed). Reason: a 200 x 200 RNAduplex returns a ~100-bp hybrid whose dG is dominated by long weak pairing, so the
frozen T1 call is "uncalibrated" and is reported descriptively only. Added test T1b (decided now, before div rows):
- Fold the whole precursor segment, leader 200 + rRNA + trailer 200 (RNAfold MFE, 37 C, 1 thread).
- Stem forms = >= 15 leader-trailer base pairs, containing a run of >= 8 pairs with at most 1-nt bulges/1x1 loops, and
  the leader-trailer pairs come within 150 nt of both mature ends.
- Same calibration bar as T1 (>= 6/7 E. coli for each rRNA, all T. bispora copies), same calls as T1.

## Amendment A2 (post hoc, written after T1 and T1b 16S results were seen)
T1b also misses the strict calibration bar (T. bispora 16S 3/4; one copy folds with 0 leader-trailer pairs). Both T1
calls are therefore reported descriptively. Post-hoc descriptive statistic, labelled as such in RESULT.md: the longest
helix (pairs with <= 1-nt bulges / 1x1 loops) in the T1 leader-trailer duplex, with its own dinucleotide-shuffle null
(200 shuffles; p_helix = fraction of shuffles with a helix >= observed). Not a frozen call.
