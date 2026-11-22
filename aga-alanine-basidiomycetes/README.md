# Fungi whose mitochondria read AGA as alanine

A gene is read three bases (one codon) at a time, and each codon names one amino acid or a stop. AGA is arginine in
every genetic code used by fungi. NCBI gives these fungi the mold mitochondrial code (table 4), where AGA is arginine.

In the mitochondria of **32 fungal genomes (about 29 species)**, AGA sits where related proteins have **alanine**, not
arginine. They are yeasts and smut fungi on three branches (subphyla) of the basidiomycete tree: *Trichosporon* and
relatives (including the human pathogen *T. asahii*) and *Tremella* on one, *Tilletia* (wheat bunt and smut fungi) on
another, *Cystobasidium* on the third. **28 of the 32** carry a tRNA for AGA with alanine's identity mark, the base
pair G3:U70; **none of 4,322** fungal genomes that read AGA as arginine and have such a tRNA carries that mark. AGA most
likely means alanine here. No known fungal genetic code reads AGA as alanine.

![Per fungal mitochondrial genome, AGA codons at conserved arginine positions (x) against conserved alanine positions (y). 32 genomes sit on the y axis, 28 of them with a tRNA for AGA carrying G3:U70; 4,387 that put AGA at 3 or more arginine positions and never at alanine sit along the x axis, none with that tRNA.](figures/aga_alanine.png)

*Each point is one of 4,638 fungal mitochondrial genomes. Blue: the genome has a tRNA for AGA (anticodon UCU) with
the alanine mark G3:U70. Two more carriers, off the axis, put AGA at 1–3 arginine positions.*

---

- **AI-generated.** The census, the analysis and the write-up were done by Claude (Anthropic) agents, run and curated
  by Ryan Xu.
- **Checked by AI, not yet by a human expert.** A separate agent re-checked the call blind (told only to infer each
  genome's code; its folder name did mention AGA). A fresh agent audited the code and reproduced the signal without
  the census alignments, and further fresh agents re-derived every number before posting and checked the lines added
  after.
- **Protein evidence is one peptide.** One peptide at one AGA site, in a protein with no known relatives, carries
  alanine. It is not from one of the conserved positions behind the call.
- Corrections and criticism are very welcome: ryan@ryanxu.net.

### What was found

| group | genomes (all in the order) | AGA at Ala positions | AGA at Arg positions | other Arg codons at Ala | other Arg codons at Arg | p |
|---|---|---|---|---|---|---|
| *Trichosporon* and relatives | 22 | 16 | 0 | 0 | 53 | 5e-16 |
| *Tremella* and relatives | 32 | 11 | 0 | 0 | 55 | 9e-13 |
| *Tilletia* | 8 | 6 | 1 | 0 | 51 | 2e-7 |
| *Cystobasidium* | 5 | 5 | 3 | 0 | 50 | 1e-5 |

Every mitochondrial genome of each order, whatever its AGA use; each position counted once; positions where at
least 90% of about 22,700 mitochondrial genomes (weighted by family) have one amino acid. "Other Arg codons" are CGN
and AGG. One-sided Fisher test.

In the 32 genomes, AGA sits at 1 to 9 such alanine positions per genome and at no arginine position (that is how the
32 were picked; the table above tests it without that rule): 124 alanine sites in total, at 24 distinct positions.
Six of the 32 rest on a single site. At these conserved positions AGA sits where relatives have alanine 78% of the
time (the rest mostly glycine and serine); the usual alanine codons, GCN, manage 97%.
Recount: [`src/check_numbers.py`](src/check_numbers.py).

### Evidence

- **At all positions, not just conserved ones.** In the 32, AGA sits where relatives most often have alanine 30% of
  the time (of 1,321 AGA codons) and arginine only 1%; the CGN arginine codons sit at arginine 74% of the time. In
  related genomes that read AGA as arginine, AGA sits at arginine 76% of the time. At only the less-conserved
  positions, none of which were used to pick the 32, AGA still sits at alanine 23% and arginine 1% of the time (1,162
  codons).
- **The tRNA.** The alanine-loading enzyme recognises one base pair in its tRNAs, G3:U70. 28 of the 32 genomes have a
  tRNA with AGA's anticodon (UCU) and that pair, all with A73, the other alanine mark. Of the 4,322 fungal
  mitochondrial genomes that put AGA at 3 or more arginine positions (and none at alanine) and have a UCU tRNA, none
  has that pair.
- **Blinded re-check.** A separate agent with its own code, told to infer each genome's code (it tested all 64
  codons), called AGA = alanine in all 4 test genomes and in none of 4 relatives used as controls.
- **Audit.** A fresh agent found no code error behind the signal, and reproduced it without the census alignments
  (searching each genome with reference proteins): AGA at alanine 71% of the time, against 93% for the usual alanine
  codons; at arginine once in 6 test genomes (a search artefact: the real codon there is CGA) and 17 of 18 times in 3
  relatives.
- **Not a sequencing error.** Raw reads behind 4 of the assemblies (*Tremella mesenterica*, *Trichosporon asteroides*,
  and *Tilletia horrida* and *Salmacisia buchloeana* from the *Tilletia* group; none yet from *Cystobasidium*) show
  AGA at every site where relatives most often have alanine: 30, 25, 25 and 26 sites, with 99.7% to 100% of reads
  agreeing.
- **Proteins.** In a public proteome of two *T. asahii* isolates from one patient (iProX, PXD025163; three sample sets,
  all searched), one peptide matched with alanine at an AGA codon: TPAPNPDHVPEDHTR, from a mitochondrial gene of
  unknown function, in 4 spectra from 2 of the 3 sets. The arginines at both of its cut sites are CGN codons. Its
  arginine version never matched. As controls, 65 spectra matched peptides that need CGN = arginine (from cox2 and
  atp6), and none matched peptides that need CGN = alanine. In the 5 of 6 searched *T. asahii* assemblies that
  contain it (4 strains), the peptide's sequence lies only on mitochondrial contigs, always with AGA at that alanine.

### What was known

- **AGA in fungi.** No NCBI code (tables 1–33) reads AGA or AGG as alanine. Papers on these genomes call their UCU tRNA
  tRNA-Arg (e.g. Li et al. 2023, *IMA Fungus*, *Apiotrichum*). A 2026 study of 328 fungal mitochondrial genomes
  (Ahrendt et al., *Nucleic Acids Res.*), which named tRNAs by their anticodon (so a reassigned one would still be
  called Arg), reports a tRNA-Arg(UCU) in all fungi except chytrids.
- **The same change elsewhere.** AGA and AGG both read as alanine in the mitochondria of the pedinophyte green algae
  Marsupiomonadales (Turmel et al. 2020, *Mitochondrial DNA B*; Barcytė et al. 2025, *PLoS Genet.*), read by a
  tRNA(UCU) with the alanine marks G3:U70 and A73, the same route as here. The same paper reports AGA and AGG as
  probably alanine (glycine not excluded) in the plastids of dinoflagellates whose plastids came from a pedinophyte.
  In Sphaeropleales it is AGG, sometimes with AGA (Žihala & Eliáš 2019, *Genome Biol. Evol.*), and in some of them AGG
  is read by a tRNA(CCU) with the alanine mark (Noutahi et al. 2019, *Mol. Biol. Evol.*). So this is a new fungal case
  of a known route, not a new kind of change. The only alanine reassignment known in fungal mitochondria is CUU/CUA in
  *Ashbya gossypii* (Ling et al. 2014).

### Open points

- **Whether the proteins carry alanine at the conserved sites.** The one peptide is from a protein without relatives,
  not from one of the conserved positions behind the call.
- **AGG**, the other codon a UCU tRNA can read, is almost absent from these genomes (12 aligned AGG codons across all
  32, none at a conserved position), so its meaning can't be called.
- **How many origins.** The branches are in three subphyla, which points to separate origins, but they share 9 of
  their AGA-at-alanine positions where about 4 are expected by chance (P ≈ 7e-4). Shared ancestry and mutational
  hotspots are both open; we have not built a tree.
- **Why AGA.** In *Tremella* and *Trichosporon* most relatives use no AGA at all, so the codon may have disappeared
  first and then come back as alanine.

### What would settle it

A second peptide, ideally over one of the conserved AGA sites. A mitochondrial proteome of one of the 32 (for
example *T. asahii*, whose proteome gave the first peptide) with alanine at these positions would settle it.

### Files

- `notes.md`: the agents' full write-up (method, evidence, checks, the protein search, references)
- `figures/`: the figure above
- `data/`: per-genome AGA placement and tRNA tables, the codon contrasts, read tallies, the proteome tallies, and the
  check reports; see [`data/README.md`](data/README.md)
- `src/`: the code; see [`src/README.md`](src/README.md). `uv run python -I src/check_numbers.py` recounts the
  numbers above from `data/`, and `uv run python -I src/aga_ala/fig_aga.py data/census/trna_vs_usage_fungi.tsv out.png`
  draws the figure. The census itself needs files not included here.
