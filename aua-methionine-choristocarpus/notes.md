# Notes: a brown alga whose mitochondria read AUA as methionine

The agents' full write-up (method, evidence, checks, open points). The short version is in [`README.md`](README.md).
Times are US Eastern, November 2026. File paths are in this folder.

**Summary.** The mitochondria of the brown alga *Choristocarpus tenellus* use the codon AUA for methionine, not
isoleucine. As far as we could find, no brown alga was known to do this. Its sister genus and five other brown algae
checked here use AUA for isoleucine. **The evidence is codon usage and RNA. No proteome of this alga exists, and both
genome sequences most likely come from one strain.**

## Background

Genes are read three letters (one codon) at a time. In the standard genetic code, AUA means isoleucine and AUG means
methionine. Some mitochondria read AUA as methionine instead: most animal mitochondria, some yeast mitochondria, and,
among the stramenopiles (the group of brown algae, diatoms and water moulds), the yellow-green algae and the
labyrinthulids. Brown algae (kelps, *Sargassum*, *Fucus*) have been taken to use the standard code in their
mitochondria. NCBI gives *Choristocarpus* the standard code, but it also gives that code to the yellow-green algae and
to the labyrinthulids as a group, where AUA for methionine is published, so the assignment is not evidence either way.

*Choristocarpus tenellus* is a small filamentous brown alga in Discosporangiales, the order that branches off first in
brown algal trees (NCBI has no such order and files the genus under Sphacelariales).

## How it was found

The lab's organelle-reading run reads organelle genomes against reference proteins without assuming a code. A test
over about 20,000 GenBank and RefSeq organelle genomes declared to read AUA as isoleucine (20,336 records) flagged 19
records where AUA sits at methionine positions; all but *Choristocarpus* were already known (thraustochytrid
mitochondria, green-dinoflagellate plastids, *Pycnococcus*/*Pseudoscourfieldia* and *Sphaerosorus* mitochondria, red
algae of the Cyanidiophyceae, the *Chloroparvula pacifica* plastid). It is not a census: of 7,728 records declared to
read AUA as methionine, only 58 had enough AUA to test, and it flagged 1. (`data/screen/`)

The follow-up rules were written before the genome downloads were read (`data/checks/RULES.md`, with amendments
timed before each step).

## Finding

- **Per position** (`src/read/b_read.py`; `data/codons/b_read.out`). The genes were aligned (tblastn) to the
  mitochondrial proteins of 8 brown algae from 8 different orders, all of which use AUA for isoleucine: *Cutleria*
  NC_085327, *Desmarestia* NC_007684, *Dictyota* NC_007685, *Pylaiella* NC_003055, *Fucus* NC_007683, *Laminaria*
  NC_004024, *Scytothamnus* NC_085336, *Sporochnus* NC_085337. A position is methionine (or isoleucine) when at least
  2 references align there and at least 2/3 of them agree; reference residues 1-10 (start codons) are skipped.
  - There are 94 AUA codons at positions where the references agree on methionine or isoleucine. Of these, 88 are at
    methionine positions.
  - Its usual isoleucine codon AUU shows the reverse: 245 at isoleucine positions and 3 at methionine.
  - With a stricter rule (at least 5 references, all agreeing), the AUA count is 68 at methionine and 0 at isoleucine.
- **In relatives** (figure; `data/codons/b_read.out`, `b_read_context.out`). For a control from a reference order,
  that order's reference was dropped.
  - The sister genus *Discosporangium* (Genoscope assembly GCA_964200235.2) puts 3 of 54 AUA at methionine positions.
  - Five brown algae from four other orders put 0 to 4 of 58 to 127 there.
  - The yellow-green alga *Sphaerosorus*, from a class known to read AUA as methionine, scores 81 of 87.
- **Methionine positions:** isoleucine-type codons (AUA, AUU, AUC) fill 67% of this genome's conserved methionine
  positions, against 2 to 7% in *Discosporangium*, *Phaeostrophion* and the kelp *Saccharina*. At the 88 methionine
  positions where it has AUA, *Saccharina* has AUG every time (`data/checks/kill/k3_columns.out`).
- **Across the genome:** AUA is at methionine positions in 20 of the 21 genes that have a conserved methionine coded
  AUA or AUG.
- **Without any alignment** (`data/composition/comp.out`): this genome has the lowest AUG use of 212 stramenopile
  mitochondrial genomes: 6.5 per 1,000 codons, against 16.0 to 22.9 in the 95 other brown algae. Its AUU and AUC alone
  match the isoleucine content of the others.
- **The RNA** (`data/rna/`): RNA was sequenced from the same culture (ENA ERR12954301). 13.7 million reads, about 4%
  of the run, were mapped to the genome. At the 86 methionine-position AUA sites with coverage, 796 reads have A and
  none have G. So there is no sign that the cell edits AUA into AUG. The library is stranded: 126 of 136 reads that
  start in a gene run in the gene's direction, so the reads are mostly RNA.

## Checks

- **Two genome sequences.** One is the Genoscope assembly (GCA_964200185.2, contig CAXNAD020000071.1), used for the
  counts above. The other is the GenBank genome (ON552980 / NC_085326), which gives 88 at methionine and 5 at
  isoleucine.
  - A 608-letter stretch of cox1, sequenced separately from this strain in 2014 (GenBank AB899175), is identical to
    both.
  - It carries AUA at 4 conserved methionines, including cox1's conserved "PRMNN" stretch. So these AUA codons are not
    errors of either assembly.
- **Re-checks.** Two separate agents re-ran the codon counts with their own code and got the same numbers.
  - The kill agent (`data/checks/kill/REPORT.md`) re-implemented the reading (0 translation mismatches), ran a
    pipeline without tblastn (83/5 on NC_085326), checked column bias, position in the protein and spread over genes,
    and compared the two genome sequences. Nothing was killed.
  - It also aligned the genes to 2,432 proteins from 62 non-brown stramenopile genera that took no part in the call.
    There, AUA faces methionine 39% of the time (AUG 44%), and AUU faces methionine 2% of the time.
  - The other (`data/checks/factcheck.md`) recounted with its own tblastn and found the library is stranded; the
    second fact-check (`factcheck2.md`) corrected the strandedness count to 126 of 136.
- **Fact-checks.** Three fresh agents checked the write-up sentence by sentence (`factcheck.md`, `factcheck2.md`,
  `factcheck3.md`); their fixes are applied here.

## What was known

- **AUA = methionine in stramenopile mitochondria:** yellow-green algae (Ehara et al. 1997, *J. Mol. Evol.*);
  labyrinthulids and the uncultured group MAST-8b (Žihala, Salamonová & Eliáš 2020, *Mol. Phylogenet. Evol.*).
- **Brown algae:**
  - Ehara et al. 1998 (*J. Phycol.*) found the standard AUA reading in the five brown algae it sequenced, none from
    Discosporangiales.
  - The *Choristocarpus* mitochondrial genome is GenBank ON552980 (RefSeq NC_085326; Choi et al., submitted 2022,
    unpublished), translated with the standard code. A 2026 paper that uses it (Hwang et al., *J. Phycol.*) says
    nothing about the code.
  - Separate literature searches by two agents found no report of a changed code in *Choristocarpus* or
    *Discosporangium* (`data/checks/priorart.md`, `data/screen/priorart_k1.md`).
- **Reference databases** (a separate agent, 14 Nov; `data/checks/dbcheck/DBCHECK.md`): none records AUA as
  methionine for any brown alga. Those that give a reading for *Choristocarpus* use the standard code, with no note or
  check behind it.
  - NCBI Taxonomy gives the standard code. The 35 protein-coding genes in GenBank and RefSeq, the 70 NCBI Protein
    records and the 35 UniParc sequences all read AUA as isoleucine, with no note about the code. UniProtKB has no
    entry.
  - The tRNA databases cannot tell: tRNAdb keeps only one CAU tRNA per brown alga (61 species) and labels each one
    methionine, and GtRNAdb has no brown algae.
  - No NCBI code table fits. Every table with AUA for methionine also reads UGA as tryptophan, and *Choristocarpus*
    ends 6 genes with UGA. If AUA does mean methionine here, annotating it correctly would need a new table.

## Open points

- **Whether the proteins carry methionine.** No proteome of this alga exists. RNA shows no editing, but only peptides
  can show which amino acid is built in.
- **Too much methionine.** If AUA means methionine, this genome has 35.5 methionines per 1,000 codons, more than any
  other brown alga (at most 22.9).
  - About 70% of the excess is in the ribosomal-protein genes and *tatC* (`data/composition/groups.out`). There, AUA
    faces leucine in relatives (12 of 38 sites) about as often as methionine (9 of 38) (`cons.out`).
  - In the respiratory genes, AUA faces methionine at 79 of 108 sites.
  - One reading is that AUA has spread into variable leucine positions in an AT-rich genome. It does not fit AUA
    meaning isoleucine either: AUA almost never sits at isoleucine positions.
- **Which tRNA reads AUA.** Brown algae normally have an isoleucine tRNA with anticodon CAU, presumably modified to
  read AUA. We found no such tRNA in *Choristocarpus*.
  - A separate agent ran two tRNA finders (tRNAscan-SE, ARAGORN) on both genome sequences
    (`data/checks/dbcheck/trna/`). Every run finds two CAU tRNAs, one like the start-codon methionine tRNA and one like
    the other methionine tRNA. The place where brown algae keep the isoleucine CAU tRNA (between the tRNAs for
    tryptophan and glutamine) is empty. Three standard brown algae have it there.
  - One tRNAscan-SE mode (organelle) also reports a weak isoleucine tRNA with anticodon AAU (score 20.2, malformed
    anticodon stem). The other runs do not find it, and it has not been checked against RNA.
  - GenBank labels both *Choristocarpus* CAU tRNAs "trnM". That is not evidence: the same submitters also label
    *Phaeostrophion*'s isoleucine CAU tRNA "trnM".
  - The isoleucine CAU tRNA is also missing from that place in *Discosporangium*, which reads AUA as isoleucine.
    *Discosporangium* keeps one extra CAU tRNA that *Choristocarpus* lacks. One tRNAscan-SE mode calls it the
    AUA-reading isoleucine type, but that mode gave the same call for none of the three standard brown-algal
    isoleucine CAU tRNAs, so the label is weak.
- **When it changed.** *Discosporangium* uses AUA for isoleucine, so the change most likely happened in the
  *Choristocarpus* lineage after the two genera split. We have not built a mitochondrial tree.
- **How widespread it is.** Most likely one strain.
  - The GenBank genome gives no strain, but its place and date (Villefranche-sur-Mer, 1980) match strain KU-1152.
  - The Genoscope culture was first thought to be a 2004 Greek *Discosporangium* culture (KU-2346). Its sample record
    now calls it *Choristocarpus*, "probably identical to KU-1152".
  - The two mitochondrial sequences differ at 27 of 32,351 positions, 25 of them in the lengths of single-letter runs
    (`data/checks/kill/k5_diff.out`).
  - No other *Choristocarpus* isolate has public sequence.

## Risks

- One strain, one genus. A problem shared by every sequencing of this culture would carry through.
- No protein data and no tRNA charging data. "Uses AUA for methionine" is a comparative codon-usage claim, the usual
  kind for mitochondrial codes, but not proof of what the ribosome builds in.
- The genome's submitters (Choi, Kawai, Yoon and colleagues) have an unpublished organelle paper and may already know.

## References

- Ehara M, Hayashi-Ishimaru Y, Inagaki Y, Ohama T (1997) *J. Mol. Evol.* 45:119-124 (doi:10.1007/PL00006210).
- Ehara M, et al. (1998) *J. Phycol.* 34:1005-1008.
- Žihala D, Salamonová J, Eliáš M (2020) *Mol. Phylogenet. Evol.* 152:106908.
- Hwang et al. (2026) *J. Phycol.* (uses NC_085326 as a comparison genome).
