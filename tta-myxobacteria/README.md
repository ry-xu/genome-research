# Bacteria that end their genes at TTA, a leucine codon

Every known bacterium ends its genes with TAA, TAG or TGA. Some have dropped one of the three, but none was known to
end genes at any other codon. A sense codon (one that means an amino acid) becoming a stop was known in
mitochondria and plastids, and proposed in 2026 for CGA in Asgard archaea from its absence in their proteins.

We screened all 199,923 bacterial and archaeal species-representative genomes in GTDB (release R232) for codons that
pile up exactly where related proteins end. One branch of myxobacteria (family Polyangiaceae) ends its genes at **TTA**,
a codon read as leucine in every other bacterial genetic code. The branch lost the tRNA that reads TTA, and its stop
protein RF2 carries a pair of changes found in no other GTDB species representative. One member, strain 1-FT3.2, is
kept alive in a lab in Finland.

---

- **AI-generated.** The screen, the analysis and the write-up were done by Claude (Anthropic) agents, run and curated
  by Ryan Xu.
- **Checked by AI, not yet by a human expert.** Separate agents re-checked the result by different methods (a blinded
  re-scoring of shuffled, renamed genomes; tblastn on raw contigs with no gene calls; an adversarial review), and
  several rounds of agent fact-checks went through the claims and the prior work.
- **Computational only.** Almost all genomes were assembled from environmental DNA (metagenome-assembled genomes,
  MAGs) by other groups. 1-FT3.2 is kept in co-culture with its prey but its proteins have not been measured. Whether
  TTA is a true stop or a codon nothing reads is not known.
- Corrections and criticism are very welcome: ryan@ryanxu.net.

### What was found

| genome | genus (GTDB) | from | TTA ends* | raw reads with TTA | RF2 motif | tRNA-Leu(UAA) |
|---|---|---|---|---|---|---|
| GCA_963169165.1 | JADLFU01 | anammox reactor, Sweden | 24/269 | 10/10 | SPY | absent |
| GCA_020849655.1 | JADLFU01 | anammox reactor, Sweden | 21/260 | 18/18 | SPY | absent |
| GCA_964605685.1 | JANWXX01 | hot spring, India | 46/309 | 20/20 | SPY | absent |
| GCA_025057345.1 | JANWXX01 | alkaline hot spring, Malaysia | 36/245 | 15/15 | SPY | absent |
| GCA_025453595.1 | JANWXX01 | wetland, California | 63/328 | no reads public | SPY | absent |
| GCA_035762175.1 | DASTPP01 | farm soil, China | 24/524 | 16/16 | SPY | absent |
| GCA_046351315.1 | JAXLYQ01 | wastewater, Hong Kong | 32/208 | 11/11 | SPF | absent |
| GCA_964396175.1 | CAYNEK01 | wastewater, China | 6/196 | not run | SPY | absent |
| their 8 nearest relatives | 4 genera incl. *Sorangium* | | 1/2,025 | | SPF (8/8) | present (8/8) |

\*Of the genes that run past the point where a close relative's version stops, how many have TTA as the next codon.

Outside GTDB, 7 more genomes pass the same frozen test, from activated sludge, wastewater, an anammox system and
1-FT3.2, a predator of the soil bacterium *Mucilaginibacter cryoferens* from subarctic peat in Finland
([Demina et al. 2025, *F1000Research* 14:1153](https://doi.org/10.12688/f1000research.171863.2)). 1-FT3.2 has the
strongest signal of all (102 gene families), and its own reads carry TTA at all 104 sites tested. Across 12 genomes,
reads carry TTA at 285 of 285 testable sites. All genomes: [`data/genomes.tsv`](data/genomes.tsv).

### Evidence that TTA ends genes

- **Where relatives stop, these genes stop at TTA.** 252 of 2,339 gene ends against 1 of 2,025 in the nearest
  relatives. By observed/expected, TTA ranks first of the 61 sense codons in every branch genome.
- **TTA swaps with the normal stops.** Take a gene that ends at TTA in one branch genome and look at the same gene in
  another genus of the branch: it ends at TAA, TAG or TGA 78% of the time and at TTA 12%. Across the branch, TTA stands
  at 777 of 18,007 gene ends; in the relatives, at 3 of 38,417.
- **A blinded re-check.** A separate agent got the 8 branch genomes and their 8 nearest relatives, shuffled and
  renamed, and scored all 61 sense codons with its own method (blastx of raw contigs, no gene calls). Its 8 strongest
  TTA scores were exactly the 8 branch genomes (1 in 12,870 by chance). No other codon reached even its lowest tier in
  any genome.
- **It is not an assembly artefact.** Raw sequencing reads carry the TTA at every testable site. The DNA after a TTA
  end evolves like non-coding DNA, the bases around TTA ends look like stop-codon context rather than leucine context,
  and the signal is spread over about 100 gene families, not a few contaminating contigs.
- **The thresholds were fixed before scoring.** The screen's pass rule and four checks per candidate were written down
  before the genomes were scored, and every later change is logged with its time ([`data/prereg.md`](data/prereg.md)).

### What else changed in these genomes

- **The TTA tRNA is gone.** None of the 8 branch genomes has a detectable tRNA-Leu(UAA), the tRNA that normally reads
  TTA (two tRNA finders agree). All 8 relatives have it (one-sided Fisher P = 7.8e-5). 1-FT3.2 lacks it too.
- **RF2 changed.** RF2 is the protein that ends translation at TAA and TGA. Among 205 genes shared by all 16 genomes,
  RF2 changed most on the way into the branch (10 fixed differences; next two genes, 8 each). Two of the changes sit
  in and just before its codon-reading motif: SPF→SPY and *E. coli* L201F. Together they occur in 7 of the 8 branch
  genomes, in 1-FT3.2, and in none of the other 161,070 GTDB representatives with an RF2.

So TTA either became a fourth stop, read by the altered RF2, or became a codon nothing reads, whose stalled ribosomes
are rescued. Genomes alone cannot tell these apart.

### What was known

- **Sense to stop, in organelles only.** Mitochondria: UUA in labyrinthulean protists such as *Thraustochytrium*, UCA
  in *Scenedesmus*, AGA/AGG in vertebrates, and UUA/UUG proposed in some green algae. Plastids: UUA and UCA in green
  dinoflagellates, where an altered plastid RF2 is proposed to read them
  ([Barcytė, Žihala & Eliáš 2025, *PLoS Genet*](https://doi.org/10.1371/journal.pgen.1011901)), the closest precedent
  for the route proposed here.
- **Archaea.** [KACI (Melnykov 2026, *PNAS*)](https://doi.org/10.1073/pnas.2610659123) describes, in its significance
  statement, "the first case outside of eukaryotic organelles of a sense codon reassigned to serve as translation
  stop". In the preprint that case is CGA in Asgard archaea, inferred because CGA never appears in their protein domains
  and its tRNA is missing; the method "does not allow direct inference of stop codons".
- **Bacteria.** Known code changes turn a stop into an amino acid (TGA in *Mycoplasma* and others; TAG in
  [our earlier screen](../tag-bacteria)) or one amino acid into another
  ([Shulgina & Eddy 2021](https://elifesciences.org/articles/71402)). None turns an amino acid into a stop.
  Genome-scale code screens (Codetta, gTranslate, KACI) ask what amino acid a codon means and do not test gene ends.

As far as we can find, this is the first bacterial sense-to-stop codon, and the first anywhere outside organelles
shown from the gene ends themselves, checked against raw reads.

### What would settle it

Protein data. C-terminal peptides from 1-FT3.2, which can be grown with its prey, would show whether its proteins end
at TTA or carry leucine there. A release test of the branch's RF2 on a UUA codon would show the mechanism.

### Files

- `notes.md`: the agents' full write-up (method, evidence, checks, open points, references)
- `data/`: the genomes, the screen's output, per-check results, raw-read tallies, the blinded panel, tRNA and RF2
  tables, the search outside GTDB, and the logged pass rules; see [`data/README.md`](data/README.md)
- `src/`, `scripts/`: the analysis, runnable with `uv run`. Parts of it reuse code and reference data from
  [`../tag-bacteria`](../tag-bacteria); `scripts/setup_reference_data.sh` says how
