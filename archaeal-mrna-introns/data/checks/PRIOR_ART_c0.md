# Prior art c0: spliced 22-nt BHB intron in a *Methanosuratincola* radical SAM gene (20 Nov 2026)

**Verdict: NEW as a specific case.** The kind of event is already known (cbf5). I found no report of a spliced intron in any radical SAM gene, in Methanomethylicia, or one that joins two annotated ORFs. Closest prior work: the cbf5 intron (Watanabe 2002; Yoshinari 2006; Yokobori 2009) and this lab's own Asgard L02 post (`outbox/2026-11-14-asgard-protein-intron`). Builds on `notes/PRIOR_ART.md`.

## 1. Protein-gene introns other than cbf5
- **Liang Y, Qi W, Dong X, Li J. 2025. *Microbiol Mol Biol Rev* 89:e00318-24, doi:10.1128/mmbr.00318-24** (full text read on PMC12713394). It says archaeal mRNAs "are generally intron-less". It credits EndA with removing introns from rRNA, SRP RNA "and Crenarchaeota cbf5 pre-mRNA". Its only other mRNA case is the HVO_1309 **5′ UTR** BHB (Schwarz 2020).
- Candidates from RNA-seq read patterns, none confirmed (from earlier notes, not re-read today): Doose 2013 (SSO1586, 144 nt in-frame), Berkemer 2015 (*M. kandleri* circles), Danan 2012 (circles but "no spliced reads").
- Europe PMC (incl. preprints) 2015–2026, Asgard/Bathyarchaeia/Korarchaeota/Methanomethylicia + intron/splicing: only tRNA (Goswami 2026 bioRxiv), group II (Mattick 2024 bioRxiv) and circular 23S rRNA (Shi 2026 *Genome Biol*, title only). No mRNA introns.

## 2. LCB70 / *Methanosuratincola* papers
- **Kohtz AJ, Petrosian N, Krukenberg V, Jay ZJ, Pilhofer M, Hatzenpichler R. 2024. *Nature* 632:1118–1123, doi:10.1038/s41586-024-07631-6.** This is the paper behind PRJNA916083 (the CP186953 record cites it). Main text paywalled; read abstract, data/contributions sections, SI PDF, Reporting Summary, SI tables. Transcriptomics: "3 biological replicates each". SI Table 4 gives TPMs for "key annotated genes" only. **No mention of introns, splicing, frameshifts or radical SAM.** The 2023 preprint (Research Square, doi:10.21203/rs.3.rs-2500102/v1, full text read) has no transcriptome data and no such terms.
- **Wu K et al. 2024. *Nature*, doi:10.1038/s41586-024-07728-y** (sister strain LWZ-6; WUR repository full text read): no intron, splicing or frameshift terms.
- Genome caveat: CP186953 was assembled by metaFlye from Nanopore reads (14×). The preprint methods say it was polished with Illumina reads (Polypolish, Polca). Check that short DNA reads carry the 22 nt.

## 3. Radical SAM / mRNA processing by the tRNA endonuclease
- I found no BHB or other intron in any radical SAM gene in archaea, bacteria or phage. (From memory, not cited: phage T4 has a group I intron in *nrdD*. That gene codes the glycyl-radical enzyme, not its radical SAM activase.)
- mRNA cut by EndA: cbf5 (cut in vitro, Yoshinari 2006) and the HVO_1309 5′ UTR (Schwarz 2020, cited via Liang 2025). Nothing newer.

## 4. GitHub / preprints
- GitHub repo/code search (Methanosuratincola, LCB70, PRJNA916083, ACO0C9_06155/06160, archaeal mRNA intron, BHB intron protein): nothing. CP186953/GCA_049570475 appear only in accession lists.
- Europe PMC preprints 2024–2026 (archaea + intron/splicing; 44 hits): none on mRNA introns.

## Not accessed
- The Kohtz 2024 main text and methods (paywall); the Montana State ScholarWorks thesis/PDF (connection failed).
- bioRxiv itself (searched via the Europe PMC preprint index only); Kohtz Extended Data figures.
