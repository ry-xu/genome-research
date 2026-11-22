# A brown alga whose mitochondria read AUA as methionine

A gene is read three bases (one codon) at a time, and each codon names one amino acid or a stop. In the standard
genetic code AUA means isoleucine and AUG means methionine. Most animal mitochondria, some yeast mitochondria and,
among the stramenopiles, the yellow-green algae and the labyrinthulids read AUA as methionine instead. Brown algae
(kelps, *Sargassum*, *Fucus*) have been taken to use the standard code in their mitochondria.

In the mitochondria of the brown alga ***Choristocarpus tenellus***, AUA sits where related proteins have methionine:
**88 of 94** AUA codons at positions where the reference proteins agree on methionine or isoleucine are at methionine
positions. Its sister genus *Discosporangium* and five other brown algae checked here use AUA for isoleucine. RNA
from the same culture carries AUA at these sites, so the cell does not edit it into AUG. As far as we could find, no
brown alga was known to read AUA as methionine.

![Bars: share of each genome's AUA codons at conserved methionine positions. Choristocarpus 88 of 94; Discosporangium 3 of 54; Phaeostrophion 2 of 71; Protohalopteris 3 of 58; Sargassum 4 of 122; Saccharina 1 of 117; Petalonia 0 of 127; the yellow-green alga Sphaerosorus 81 of 87.](figures/aua_met_bars.png)

*AUA codons at conserved methionine positions (orange) and isoleucine positions (grey) in seven brown-algal
mitochondrial genomes and the yellow-green alga* Sphaerosorus*, from a class known to read AUA as methionine.*

---

- **AI-generated.** The screen, the analysis and the write-up were done by Claude (Anthropic) agents, run and curated
  by Ryan Xu.
- **Checked by AI, not yet by a human expert.** Two separate agents re-ran the codon counts with their own code and
  got the same numbers. A kill agent tried to break the claim (code bugs, column bias, assembly artefacts, RNA
  editing, tRNAs) and could not. Three fresh agents fact-checked the write-up sentence by sentence, two searched the
  literature for prior work and one checked the reference databases.
- **Codon usage only, one strain.** No proteome of this alga exists, so which amino acid the ribosome builds in is not
  shown. Both genome sequences most likely come from one strain.
- Corrections and criticism are very welcome: ryan@ryanxu.net.

### What was found

The genes aligned to the mitochondrial proteins of 8 brown algae from 8 different orders, all of which use AUA for
isoleucine. A position counts as methionine (or isoleucine) when at least 2 references align there and at least two
thirds of them agree.

| genome | AUA at methionine positions | AUA at isoleucine positions |
|---|---|---|
| ***Choristocarpus tenellus*** | **88** | **6** |
| *Discosporangium* (sister genus) | 3 | 51 |
| *Phaeostrophion* | 2 | 69 |
| *Protohalopteris* | 3 | 55 |
| *Sargassum* | 4 | 118 |
| *Saccharina* (kelp) | 1 | 116 |
| *Petalonia* | 0 | 127 |
| *Sphaerosorus* (yellow-green alga, positive control) | 81 | 6 |

*Choristocarpus*'s usual isoleucine codon AUU shows the reverse: 245 at isoleucine positions and 3 at methionine.
With a stricter rule (at least 5 references, all agreeing), the AUA count is 68 at methionine and 0 at isoleucine.
Recount: [`src/check_numbers.py`](src/check_numbers.py).

### Evidence

- **Methionine positions.** Isoleucine-type codons (AUA, AUU, AUC) fill 67% of this genome's conserved methionine
  positions, against 2 to 7% in *Discosporangium*, *Phaeostrophion* and the kelp *Saccharina*. At the 88 methionine
  positions where it has AUA, *Saccharina* has AUG every time.
- **Across the genome.** AUA is at methionine positions in 20 of the 21 genes that have a conserved methionine coded
  AUA or AUG.
- **Without any alignment.** This genome has the lowest AUG use of 212 stramenopile mitochondrial genomes: 6.5 per
  1,000 codons, against 16.0 to 22.9 in the 95 other brown algae. Its AUU and AUC alone match the isoleucine content
  of the others.
- **No RNA editing.** RNA was sequenced from the same culture. In 13.7 million reads (about 4% of the run), at the 86
  methionine-position AUA sites with coverage, 796 reads have A and none have G.
- **Two genome sequences.** The counts above use the Genoscope assembly; the GenBank genome gives 88 at methionine
  and 5 at isoleucine. A 608-letter stretch of cox1, sequenced separately from this strain in 2014 (GenBank
  AB899175), is identical to both and carries AUA at 4 conserved methionines, including cox1's conserved "PRMNN"
  stretch. So these AUA codons are not errors of either assembly.
- **Outside brown algae.** Aligned to 2,432 proteins from 62 non-brown stramenopile genera that took no part in the
  call, AUA faces methionine 39% of the time (AUG 44%), and AUU faces methionine 2% of the time.
- **The screen.** It came from a screen of about 20,000 GenBank and RefSeq organelle genomes declared to read AUA as
  isoleucine. Of 19 flags, all but this one were already known. It is not a census: of 7,728 records declared to read
  AUA as methionine, only 58 had enough AUA to test, and it flagged 1.

### What was known

- **AUA = methionine in stramenopile mitochondria:** yellow-green algae (Ehara et al. 1997, *J. Mol. Evol.*);
  labyrinthulids and the uncultured group MAST-8b (Žihala, Salamonová & Eliáš 2020, *Mol. Phylogenet. Evol.*).
- **Brown algae.** Ehara et al. 1998 (*J. Phycol.*) found the standard AUA reading in the five brown algae it
  sequenced, none from Discosporangiales. The *Choristocarpus* mitochondrial genome is GenBank ON552980 (RefSeq
  NC_085326; Choi et al., submitted 2022, unpublished), translated with the standard code. A 2026 paper that uses it
  (Hwang et al., *J. Phycol.*) says nothing about the code. Separate literature searches by two agents found no
  report of a changed code in *Choristocarpus* or *Discosporangium*.
- **Reference databases** (14 Nov): none records AUA as methionine for any brown alga. NCBI Taxonomy gives
  *Choristocarpus* the standard code, and the 35 protein-coding genes in GenBank and RefSeq, the 70 NCBI Protein
  records and the 35 UniParc sequences all read AUA as isoleucine. NCBI also gives the standard code to the
  yellow-green algae and to the labyrinthulids as a group, where AUA = methionine is published, so the assignment is
  not evidence either way. No NCBI code table fits: every table with AUA for methionine also reads UGA as tryptophan,
  and *Choristocarpus* ends 6 genes with UGA.
- **Where it sits.** *Choristocarpus* is a small filamentous brown alga in Discosporangiales, the order that branches
  off first in brown algal trees (NCBI has no such order and files the genus under Sphacelariales).

### Open points

- **Whether the proteins carry methionine.** No proteome of this alga exists.
- **Too much methionine.** If AUA means methionine, this genome has 35.5 methionines per 1,000 codons, more than any
  other brown alga (at most 22.9). About 70% of the excess is in the ribosomal-protein genes and *tatC*, where AUA
  faces leucine in relatives (12 of 38 sites) about as often as methionine (9 of 38). In the respiratory genes, AUA
  faces methionine at 79 of 108 sites. One reading is that AUA has spread into variable leucine positions in an
  AT-rich genome. It does not fit AUA meaning isoleucine either: AUA almost never sits at isoleucine positions.
- **Which tRNA reads AUA.** Brown algae normally have an isoleucine tRNA with anticodon CAU, presumably modified to
  read AUA. Two tRNA finders (tRNAscan-SE, ARAGORN) on both genome sequences find two CAU tRNAs, both methionine-like;
  the place where brown algae keep the isoleucine CAU tRNA (between the tRNAs for tryptophan and glutamine) is empty.
  It is also empty in *Discosporangium*, which reads AUA as isoleucine and keeps one extra CAU tRNA that
  *Choristocarpus* lacks.
- **When it changed.** *Discosporangium* uses AUA for isoleucine, so the change most likely happened in the
  *Choristocarpus* lineage after the two genera split. We have not built a mitochondrial tree.
- **How widespread.** Most likely one strain (KU-1152, Villefranche-sur-Mer, 1980). The two mitochondrial sequences
  differ at 27 of 32,351 positions, 25 of them in the lengths of single-letter runs. No other *Choristocarpus*
  isolate has public sequence.

### What would settle it

Peptides from *Choristocarpus tenellus* mitochondrial proteins (the strain is in the Kobe University culture
collection as KU-1152), showing methionine or isoleucine at these positions; and a second isolate or a close relative
in Discosporangiales.

### Files

- `notes.md`: the agents' full write-up (method, evidence, checks, references)
- `figures/`: the figure above
- `data/`: the per-site codon table, the codon counts, the RNA tallies, the codon-usage tables, the screen output and
  the check reports; see [`data/README.md`](data/README.md)
- `src/`: the code; see [`src/README.md`](src/README.md). `uv run python -I src/check_numbers.py` recounts the
  numbers above from `data/`, and `uv run python -I src/fig_aua_bars.py out.png` draws the figure. The reading itself
  needs the genome downloads, BLAST+ and the lab's organelle tables, which are not included.
