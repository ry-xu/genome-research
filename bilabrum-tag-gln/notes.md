# Notes: a marine flagellate that reads the stop codon TAG as glutamine

The agents' full write-up (screen, evidence, checks, risks). The short version is in [`README.md`](README.md). Times
are US Eastern, November 2026. File paths are in this folder; the agents' own logs, with paths in the lab's working
tree, are in `data/checks/`.

**Summary.** *Bilabrum* cf. *latius*, a small heterotrophic marine flagellate, has reassigned the stop codon TAG to
the amino acid glutamine in its nuclear genes. No change to the nuclear genetic code was known before in
stramenopiles, the large group that includes diatoms, kelp and the potato-blight pathogen. **The genome is the
Wellcome Sanger Institute's (Aquatic Symbiosis Genomics), unpublished; there is no RNA or protein from this strain.**

## Background

Genes end at a stop codon: TAA, TAG or TGA. A few eukaryotes read one or more of them as an amino acid instead, mostly
ciliates and a handful of other lineages, nearly all single-celled. TAG (often with TAA) read as glutamine has arisen
many times.

## The screen

I screened 1,785 protist genomes and transcriptomes for genes that keep matching known proteins past an in-frame stop
(1,927 units scored, 1,539 with enough aligned codons to call, 1,406 of them non-duplicate). The call rules were set
on 67 genomes and transcriptomes of known code (43 positives, 24 negatives; `data/screen/controls.tsv`) and frozen
at 21:17 on 11 Nov, before the screen ran (`src/screen/call.py`). *Bilabrum* cf. *latius*, whose genome the
Wellcome Sanger Institute released in March 2025 (GCA_965178125.1, 43 chromosomes), was the clearest new case. In the
screen: TAG at 190 inner columns in 121 proteins, top meaning glutamine; TAA and TGA 0. Stage 2 (is the codon on the organism's own genes?) passed on every test against the control codon CAA
(`data/screen/GCA_965178125.1_TAG.txt`).

The genome is a cobiont: it was sequenced with the dinoflagellate *Dinothrix rugata*. Chromosome 12 carries two rDNA
copies matching *Bilabrum latius* 18S (99.6-99.8%), the host's 18S is on one unplaced contig, and all 123 TAG-bearing
hits of stage 2 lie on chromosomes (`data/ssu/`).

## Finding

- **TAG sits inside conserved genes.** I aligned 76 conserved human proteins to the genome with the standard code
  (tblastn; `data/review/`). TAG appears at 166 aligned positions; TAA at 1 and TGA at 4. The same test finds no TAG in
  *Cafeteria*, a relative in the same order (Bicosoecida). At TAG's rate, point errors would give about 38 TAA and 40
  TGA (`data/review/bil.analysis.txt`).
- **It means glutamine.** At the 91 cleanest sites, the human protein has glutamine at 44, then lysine at 11 and no
  other amino acid at more than 5. That is the pattern of CAG, the usual glutamine codon (best-fitting codon profile:
  glutamine). The 91 sites are in 36 proteins on 30 chromosomes, among them enolase, PCNA, BiP, Sec61A, RNA
  polymerase II subunits, PRP8 and SF3B1.
- **It ends none of 149 genes checked** (TGA ends 136, TAA 13), against 11 of 145 in *Cafeteria* and 35-39 in each
  of three other bicosoecids (149-172 genes each; `data/ends/terminal_codons_relatives.tsv`). *Cafeteria* uses TAG
  least of the relatives, so it is the conservative comparison.
- **The machinery fits.** The genome has a tRNA with anticodon CUA, which pairs with UAG (`data/trna/`); identity
  models score it as glutamine- or histidine-like, not tyrosine-like. No nuclear tRNA reads TAA, which fits TAA
  staying a stop. Its eRF1, the protein that recognises stop codons, has valine at a serine that 2,714 of 2,776 eRF1s
  in UniProt reference proteomes keep and none has as valine (domain 1 TAVNIKS for human TASNIKS;
  `data/erf1/erf1_aln.fa`). The eRF1 gene itself has an in-frame TAG.
- **TAG is in the RNA.** No RNA-seq exists for this genome, but a related *Bilabrum* contaminates algal cultures
  sequenced by the 1KP project (One Thousand Plant Transcriptomes), as Vierock et al. 2022 found. Its reads, from three
  of the four cultures it contaminates (about 96% identical to this genome), carry TAG at 41 of 60 covered sites and
  CAG at 18 (`data/rna/all3.table.txt`). The three cultures carry one strain, so the runs are pooled
  (`data/rna/strain_compare.txt`).
- **TAG changes like a silent codon.** Where this genome has CAG, the 1KP strain has TAG at 12.3 per 1,000 sites,
  close to the rate of a silent leucine change (CTG to TTG, 14.4). Comparable changes that alter the protein (His to
  Tyr, Arg to Trp, Pro to Ser) run at 0.8 per 1,000 (`data/rna/ALL3.pairs.txt`). A stop codon in a conserved gene
  would be weeded out like those; TAG is not.

The figure (`figures/enolase_tag.png`, `src/fig_bilabrum_enoa.py`) shows enolase: a TAG opposite glutamine 166 of
human enolase, where *Cafeteria* has CAG (all 32 RNA reads there carry TAG), and a second TAG opposite a lysine (all
36 reads TAG).

## Known before

- **Nuclear stop codons read as sense** in ciliates, diplomonads and their relatives, some green algae, the
  trypanosomatid *Blastocrithidia*, a parasitic dinoflagellate (*Amoebophrya*), oxymonads and a few others.
- **Stramenopiles.** Swart et al. 2016 screened the MMETSP transcriptomes, which include *Cafeteria*, and found no
  change outside ciliates. Some Labyrinthulea, and the Arctic *Plocamiomonas* (Gastineau et al. 2026, *Sci. Rep.*),
  change their code in mitochondria only. *Blastocystis* uses the standard code but completes some TAA stops by adding
  the poly(A) tail (Klimeš et al. 2014). NCBI lists the standard code for every *Bilabrum* taxon.
- **This TAG has been seen.** Vierock et al. 2022 (*Sci. Adv.* 8:eadd7729) found *Bilabrum* channelrhodopsins in 1KP
  transcriptomes. Some transcripts "appeared to contain premature stop codons", so they built the protein from a
  stop-free haplotype in the reads. One such stop is a TAG where that haplotype has CAG.

## Checks

- **Fixed in advance:** the call rules were frozen before the full screen.
- **Independent re-check:** a separate agent re-derived the TAG counts with tblastn and a 76-protein human panel; same
  verdict (`data/review/`).
- **Fact-check** (second fresh agent; `data/checks/factcheck.md`): checked every sentence against the project files
  and recomputed the RNA numbers. It found a misattributed paper (Vierock, not Govorunova), an overstated tRNA line
  and a wrong genus count; all were fixed.
- **Pre-post check** (third fresh agent, 12 Nov; `data/checks/prepost_check.md`): verdict POST. It re-derived every
  body number from the raw files, re-ran the literature search (Europe PMC incl. bioRxiv, the Vierock supplement,
  NCBI code tables and taxonomy) and found nothing earlier, and checked every codon in the figure against the genomes.
- **Novelty:** no earlier report of a nuclear code change in any stramenopile was found. The only earlier sighting is
  Vierock 2022, who treated the TAG as a premature stop.

## Risks and open points

- **Not our genome to publish.** It is a Sanger Aquatic Symbiosis Genomics cobiont, sequenced with the dinoflagellate
  *Dinothrix rugata*, and unpublished. The Sanger team may report the code themselves.
- **RNA from a related strain.** The RNA comes from a strain about 96% identical (95.6-95.8%), not from this genome.
  58 of the 60 covered sites come from one 1KP run.
- **The "silent rate" rests mostly on one run.** 21 of the 23 CAG/TAG differences are from that run. The two smaller
  runs show TAG at about half the leucine rate (5.7 vs 10.9, 8.1 vs 15.8 per 1,000).
- **The tRNA is suggestive, not decisive.** By sequence it is closest to tRNA-Tyr. The identity models do not
  separate glutamine from histidine in this genome.
- **Whether TAG still ends a few genes**: no RNA or protein from this strain.
- **How far the change reaches.** Genomes of three other bicosoecid genera (*Cafeteria*, *Caecitellus*,
  *Halocafeteria*) and an unnamed relative use the standard code, so the change arose after *Bilabrum* split from
  them; unsequenced close relatives may share it, and no other *Bilabrum* genome exists.

## The rest of the screen

Every other call or watch on a lineage with no known code failed the follow-up checks: contaminants in a culture,
mostly ciliates (in red-alga and *Symbiodinium* transcriptomes); near-copies of the host's own transcripts carrying
errors or C-to-U changes (*Aduncisulcus*, *Ochromonas*); mitochondrial or plastid genes; probably bacteria that read TGA
as glycine; and known codes under uninformative names. One other genome-backed case, *Maristentor dinoferus* (TAA=E,
TAG=Q, TGA=W), is a new combination, but inside heterotrich ciliates, where stop changes are already known. A
follow-up over lineages that exist only as SRA reads found nothing new; its two stramenopiles (*Pirsonia*,
*Mutimo*) use TAG as a stop.

## Who would care

The Sanger Tree of Life Aquatic Symbiosis Genomics team; Vierock and Hegemann (Berlin) and Rozenberg and Béjà
(Technion), whose *Bilabrum* channelrhodopsin carries the TAG; the Arndt group (Cologne), who described *Bilabrum*.
A courtesy note to the Sanger team before posting was recommended, since it is their unpublished genome.
