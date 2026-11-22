# Genome research

Research findings by AI research agents, steered by Ryan Xu.

## Results

Newest first. Each folder has a short README, the agents' full write-up (`notes.md`), the data tables and the code.

- **[A fungus whose mitochondria read UGA as arginine](uga-arginine-pneumocystis/)** (November 2026). In the
  mitochondria of *Pneumocystis jirovecii*, the fungus that causes pneumocystis pneumonia, UGA (tryptophan in the
  code these fungi are filed under) sits at 23 strongly conserved arginine positions and at none of 37 tryptophan
  positions. RNA from a patient carries UGA there (7,109 of 7,298 reads), so it is not edited into an arginine codon.
  UGA most likely means arginine here. None of NCBI's genetic codes reads UGA as arginine; as far as we can find,
  this is the first such case in fungi. There is no protein evidence yet.
- **[A marine flagellate that reads the stop codon TAG as glutamine](bilabrum-tag-gln/)** (November 2026). In the
  nuclear genes of *Bilabrum* cf. *latius*, a small marine flagellate, the stop codon TAG sits at 166 aligned positions
  in 76 conserved proteins, most often where the human protein has glutamine, and ends none of 149 genes checked. RNA
  from a related *Bilabrum* strain carries the TAG too. As far as we can find, no change to the nuclear genetic code was
  known before in stramenopiles, the group that includes diatoms and kelp. The genome is from the Wellcome Sanger
  Institute's Aquatic Symbiosis Genomics project.
- **[Fungi whose mitochondria read AGA as alanine](aga-alanine-basidiomycetes/)** (November 2026). In the
  mitochondria of about 30 basidiomycete fungi, among them *Trichosporon*, *Tremella*, *Tilletia* and *Cystobasidium*,
  the arginine codon AGA sits where related proteins have alanine and never where they have arginine. 28 of the 32
  genomes carry a tRNA for AGA with alanine's identity mark G3:U70, which none of 4,322 fungal genomes that read AGA as
  arginine has. AGA most likely means alanine here; no known fungal genetic code reads it that way. Some green algae
  use the same route (Barcytė et al. 2025), and the protein evidence so far is one peptide.
- **[A brown alga whose mitochondria read AUA as methionine](aua-methionine-choristocarpus/)** (November 2026). In
  the mitochondria of the brown alga *Choristocarpus tenellus*, 88 of 94 AUA codons at conserved positions sit where
  related proteins have methionine, not isoleucine as in the standard code. Its sister genus *Discosporangium* and five
  other brown algae use AUA for isoleucine, and RNA from the same culture shows no editing of AUA into AUG. Most animal
  mitochondria read AUA this way, but as far as we can find no brown alga was known to. The evidence is codon usage
  from what is most likely one strain; no proteome exists.
- **[A red alga whose mitochondria read AGG as methionine](agg-methionine-porphyridium/)** (November 2026). In the
  mitochondria of the single-celled red alga *Porphyridium aerugineum*, AGG, arginine in its code, sits where related
  proteins have methionine: 7 times at strongly conserved methionine positions and never at conserved arginine ones.
  In two isolates sequenced by two labs, AGG swaps with the methionine codon ATG 17 times and with an arginine codon
  never. AGG most likely means methionine here. No red alga was known to do this; a calcareous sponge's mitochondria
  are reported to read AGA and AGG as methionine (Lavrov et al. 2016, as restated in Lavrov & Pett 2016 and Lukeš et
  al. 2025). Which tRNA reads AGG is unknown, and there is no protein data yet.
- **[Beetle symbionts that read TGA as tryptophan and deleted RF2](tga-tryptophan-trixagus/)** (November 2026). The
  bacterial symbionts (*Candidatus* Shikimatogenerans) of three *Trixagus* beetles read the stop codon TGA as
  tryptophan, and each has deleted RF2, the protein that stops translation at TGA. On a tree of the symbionts the three
  do not group together, so RF2 was deleted, and TGA spread into tryptophan codons, at least twice, as Parks et al.
  2026 found in Eggerthellaceae gut bacteria. GenBank annotates the genomes with the standard code, so 19-25 genes in
  each are filed as pseudogenes. The genomes were sequenced by Wierz, Malm, Kaltenpoth and Engl (Max Planck Institute
  for Chemical Ecology), whose paper on them is listed as unpublished.
- **[Bacteria that keep a second rRNA operon whose rRNA barely accumulates](second-rrna-operon-aquicultoraceae/)**
  (November 2026). Most genomes of the Aquicultoraceae, a family of uncultured groundwater bacteria, carry a second
  rRNA operon only about 79% identical to the normal one, with changed catalytic-centre bases. It has been inherited
  since before the family's genera split, and selection keeps its folded shape. Yet in the one species with RNA data,
  0 of 280 23S read pairs fall on its copy, against about one to one in the DNA. We found no bacterial rRNA operon
  reported that is kept like this but does not accumulate.
- **[Bacteria that changed the ribosome's ricin-target loop](ricin-loop-zymomonas/)** (November 2026). The
  sarcin-ricin loop of the ribosome, where ricin and α-sarcin strike, reads GAGA and is described as universally
  conserved. It reads GAAA instead in 85 of 15,739 bacterial genomes we checked, from 13 phyla, among them the
  ethanol bacterium *Zymomonas mobilis*, in every copy in all 18 of its complete genomes. The cheese bacterium
  *Propionibacterium freudenreichii* reads GCAA.
- **[Yeasts that read CUU and CUA as proline](proline-codon-yeasts/)** (November 2026). In the mitochondria of four
  yeasts, the codons CUU and CUA (leucine in the standard code) sit where related proteins have proline, 286 of 290
  times at strongly conserved positions. Proline is spelled CCN in every genetic code known so far. The change most
  likely happened twice. The evidence is DNA only, so RNA editing is not yet ruled out.
- **[Gut bacteria that skip the C of CGG](cgg-skip-tanb77/)** (November 2026). Almost all bacteria of TANB77, an
  order of gut bacteria, have no tRNA that reads the arginine codon CGG. At most CGGs in their genes the C is passed
  over and the next three bases are read as glycine, so the gene continues one base over. As far as we can find, no
  bacterium or archaeon was known to shift reading frame at the same codon throughout its genome.
- **[Bacteria that stop their genes at TTA](tta-myxobacteria/)** (November 2026). A branch of myxobacteria stops its
  genes at TTA, normally read as the amino acid leucine: the first bacteria known to use a new stop codon. They lost
  the tRNA that reads TTA, and RF2, the protein that recognizes stop codons, carries two changes that no other GTDB
  species representative has.
- **[Bacteria that no longer use TAG as a stop codon](tag-bacteria/)** (November 2026). In 12 bacterial genomes, genes
  run straight through TAG, tens to hundreds of times per genome, because these bacteria dropped RF1, the protein that
  stops at TAG. It happened 6 separate times, and TAG is most often read as glutamate.

Not yet reviewed by human experts. Corrections and criticism are very welcome: ryan@ryanxu.net.
