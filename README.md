# Genome research

Genome findings by Claude research agents, run and curated by Ryan Xu.

## Results

Newest first. Each folder has a short README, the agents' full write-up (`notes.md`), the data tables and the code.

- **[Bacteria that stop their genes at TTA](tta-myxobacteria/)** (November 2026). A branch of myxobacteria stops its
  genes at TTA, normally read as the amino acid leucine: the first bacteria known to use a new stop codon. They lost
  the tRNA that reads TTA, and RF2, the protein that recognizes stop codons, carries two changes that no other GTDB
  species representative has.
- **[Bacteria that no longer use TAG as a stop codon](tag-bacteria/)** (November 2026). In 12 bacterial genomes, genes
  run straight through TAG, tens to hundreds of times per genome, because these bacteria dropped RF1, the protein that
  stops at TAG. It happened 6 separate times, and TAG is most often read as glutamate.

---

- **AI-generated.** The screens, the analyses and the write-ups were done by Claude (Anthropic) agents, run and curated
  by Ryan Xu.
- **Checked by AI, not yet by a human expert.** Separate agents re-checked each result by different methods, and more
  agents fact-checked the claims and searched for prior work.
- **Computational only.** The genomes are public, mostly assembled from environmental or gut DNA by other groups. No
  protein data exist yet for any of these organisms.
- Corrections and criticism are very welcome: ryan@ryanxu.net.

Code is MIT and text and data are CC BY 4.0 (`LICENSE`). `tta-myxobacteria/` reuses code and reference data from
`tag-bacteria/`; its `scripts/setup_reference_data.sh` says how.
