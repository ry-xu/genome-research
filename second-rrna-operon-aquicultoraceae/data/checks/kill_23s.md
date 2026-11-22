<!-- Copied from the lab record. Paths that point to files in this repository were rewritten to repo-relative ones; any other path refers to the original working folder and that file is not included (see data/README.md). -->

# Kill report: divergent rRNA operon in Aquicultoraceae
Independent kill agent (ae3d4e10), finished 14 Nov evening. Saved verbatim by Deneb, because the agent could not write report files.

The divergent operon is real and is the organism's own DNA, and the PTC changes are real. The claim breaks only on structure: "keeps the secondary structure like normal copies" does not survive a fairer test.

The harness refused to let me write `data/checks/kill_23s.md` because subagents may not write report files, so the full report is below; you can save it as REPORT.md if you want the file. All scripts and outputs are in `src/checks/kill_23s/` (scripts) and `data/kill_23s/` (outputs; large files not included), and nothing outside that folder was changed.

## Verdicts

| Sub-claim | Verdict |
|---|---|
| A. The divergent 16S+23S operon is own DNA, not contamination or misassembly | **holds**; the decisive test is conserved gene order in 6 genera |
| B. Its 23S is ~72% identical to the same genome's normal 23S | **holds**: 71-73% by BLAST, 0.79 over aligned alignment columns |
| C. PTC changes at ≥99.8%-conserved positions, not an alignment artefact | **holds, with two corrections** (below) |
| D. Keeps rRNA secondary structure like normal copies (0.93-0.96) | **weakened**: well above a neutral null, but below every functional 23S control at similar divergence |
| E. In ~15 of 19 Aquicultoraceae genomes with a 23S, across genera and studies | **holds**, but 8 of the 15 are one species from one site |
| F. No prior report | **not contradicted** by a limited search |

## Evidence

**A. Own DNA (DASXVE01, GCA_036504445.1)**
- **Reads** (`pile.py` → `pile.out`): over DASXVE010000023.1:1-9000, read depth is at least 14 everywhere inside the operon, and mate-pair fragment coverage never drops to 0. There is no pile-up of clipped read ends at any junction (0-7 per 250 bp window; genome-wide clipped rate 3.55%).
- **PTC bases in reads at ≥97% identity:** every read agrees.
  - 2062 C, 30/30; 2251-2253 UUU, 35/35
  - 2449 A, 53/53; 2451 U, 53/53; 2452 G, 54/54; 2453 G, 51/51
  - 2503 U, 48/48; 2504 A, 48/48; 2583 A, 42/42; 2584 A, 42/42
- **Genes on the contig** (`contig_affinity.py`): all 27 genes with a hit in contig 23 (1-37 kb) have their best hit in Aquicultoraceae (whole genome: 0.94). GC is 0.498 for the contig and 0.499 for the genome.
- **Gene order** (`synteny.py` → `synteny.out`): this is the decisive test.
  - Of the 9 genes flanking the DASXVE divergent operon, 2-5 have their best hit among the flanking genes of the divergent locus in each other genus: Aquicultor_A (8 genomes), CADDYS01, DASQQM01, JASEIR01 and Kmv33.
  - They hit the normal operon's flanking genes 0 of 9 in all 23 normal-operon loci tested in other genomes.
  - The normal operon's flanking genes give the mirror image: 5-6 of 12 at normal loci, 0 of 12 at every divergent locus.
  - So the divergent operon sits at its own conserved position, inherited along with the host lineage.
- **Coverage in other genomes:** the divergent contigs sit at genome coverage in Kmv33 (11.5x vs 12.1x) and CADDYS (18.6x vs 19.2x). The exception is GCA_027345605.1, where only fragments are present, at 38x against a genome median of 114x.
- **Outside the family:** against the ~17,300 census 23S sequences, the best hits outside Aquicultoraceae are 69-74% over at most 1.3 kb.

**C. PTC changes** (`ptc_pairwise.py` → `ptc_pairwise.out`)
- An alignment with no structure model (Biopython pairwise) confirms the changes, with identical flanking sequence on both sides. For example, `ctgGGGcgg` becomes `ctgTTTcgg` (2251-2253), `cctcgATg` becomes `cctcgTAg` (2503-2504), and `ctggGTt` becomes `ctggAAc` (2583-2584).
- A2062C, G2251U, G2253U, U2504A and G2583A are shared by the DASXVE01, Kmv33 and CADDYS01 copies, so they arose before these genera split.
- **Correction 1:** two positions are below the claimed 99.8% conservation (`fcons_ptc.out`): A2453 is 99.66% and U2504 is 99.74%.
- **Correction 2:** in the Kmv33 and CADDYS01 copies, the exact base called at 2449/2452 (and 2584 in CADDYS01) depends on the aligner, because the pairwise alignment puts a 1-nt gap in the 2447-2453 loop. Under either alignment the conserved UAACA motif is absent.

**D. Structure test without the model's bias** (`pairing_free.py` → `pairing_free*.out`)
- **Why the claim's number is inflated:** the 0.93-0.96 comes from alignment to the Rfam covariance model (cmalign), which rewards base pairing.
- **Method:** I placed each divergent copy on the normal copy by plain pairwise alignment, scored the fraction of census base pairs that stay canonical, and ran 13 functional 23S genes from other taxa through the same pipeline as controls.

| Sequence (identity to DASXVE normal 23S) | Canonical pairs | One-side changes | Both-side changes (fraction canonical) |
|---|---|---|---|
| Divergent DASXVE01 (0.79) | 0.832 | 152 | 82 (0.76) |
| Divergent Kmv33 (0.79) | 0.843 | 134 | 99 (0.75) |
| Divergent CADDYS01 (0.81) | 0.864 | 135 | 82 (0.79) |
| Neutral null, mutation rate per position scaled to its census variability | median 0.70-0.71, max 0.74 | | |
| 13 functional 23S controls (0.787-0.86) | 0.869-0.930 | 52-117 | 97-140 (0.87-0.97) |

- The divergent copies are clearly under structural selection: they sit far above the neutral null.
- They are still below every functional control, with about twice as many pair-breaking changes relative to compensatory ones (ratio 1.35-1.85 vs at most 0.91 in controls).
- Their substitutions reach unpaired positions almost as often as paired ones (0.18-0.20 vs 0.21-0.22). In the controls, unpaired positions change much less (0.09-0.16 vs 0.18-0.25), so the divergent copies are changing the normally conserved core too.
- Between genera, divergent copies are 0.86-0.90 identical while normal copies of the same genus pairs are 0.936-0.952, roughly 2.2x faster evolution (`pid_div.out`, `pid_norm.out`).

**New finding, not in the claim** (`asd.py` → `asd.out`)
- All 6 divergent 16S genes, from 5 genera, lack the 3'-end anti-Shine-Dalgarno motif (CCUCC), GAUCACC and GUAACAAGG. All 4 normal 16S genes checked have all three.
- Ribosomes built from this operon could not use Shine-Dalgarno pairing, which points to a specialised or non-translating role.

**E. Prevalence** (`loci.py` → `loci_blast.tsv`; `prevalence.py` → `prevalence.out`)
- Of 23 Aquicultoraceae genomes, 18 have a 23S locus and 15 of those have a divergent copy. By genus: Aquicultor_A 8/8, CADDYS01 3/5, DASQQM01 1/1, DASXVE01 1/1, JASEIR01 1/1, Kmv33 1/1, Aquicultor 0/1.
- No divergent copy appears in any other family; the only near-misses are two short fragments of 295 and 380 bp at 68-71% identity.
- The 8 Aquicultor_A copies are identical: one species (*Ca.* Aquicultor secundus) from Crystal Geyser.
- That leaves about 7 independent populations, assembled with MEGAHIT, SPAdes and IDBA-UD, so a single assembler's artefact is ruled out.

**F. Prior art** (five web searches plus PubMed)
- I found no report of an Aquicultor rRNA paralog, or of any bacterial 23S paralog with PTC changes.
- The nearest precedent is Johansen JR et al. 2017, "Highly divergent 16S rRNA sequences in ribosomal operons of *Scytonema hyalinum* (Cyanobacteria)", PLoS ONE 12:e0186393, doi:10.1371/journal.pone.0186393. Its abstract says: "A highly divergent 16S rRNA gene was found in one of the five ribosomal operons … the highest observed in any prokaryotic species thus far (7.3-9.0%)… The secondary structure… was nearly identical, indicating possible functionality."
- That is far less divergent and concerns 16S only, but a novelty claim should cite it.
- *Ca.* Aquicultor was named by Jiao et al., ISME J 2021, 15:3005-3018 (per LPSN). I did not read its full text.

## Could not check
- **NCBI nt:** remote BLAST of the divergent 23S (dc-megablast, 60 min) and 16S (megablast, 32 min) against core_nt never returned. I stopped both by their PIDs, so there is no nt result.
- **Expression:** no RNA or protein data were checked.
- **Read-level checks** were done only for DASXVE01.
- **The Aquicultor genus:** whether it lacks the copy is open, since only 1 of its 5 genomes has any 23S.

## Most likely way the claim is wrong
The weak point is the functional framing, not reality. "Keeps the structure like normal copies" comes from scoring an alignment that rewards pairing. Under a fairer test the copies are constrained but measurably degraded:
- twice the pair-breaking changes of functional 23S at similar divergence;
- substitutions reaching the conserved core;
- about 2.2x faster evolution;
- a partner 16S that has lost its anti-Shine-Dalgarno sequence.

The defensible statement is "an ancient, inherited, structurally constrained but non-canonical rRNA paralog". The 99.8% conservation threshold should also be corrected for A2453 and U2504.