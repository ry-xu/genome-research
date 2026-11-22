# Choristocarpus AUA = Met, WGS test: result (mine only; 09:10 14 Nov; before kill and prior art)
Rules: RULES.md (+ amendments A, B). Scripts a_contigs.py, b_read.py; outputs a_contigs.out, b_read.out.

| target | ATA at M / I | ATG at M / I | ATT at M / I | ATC at M / I | call |
|---|---|---|---|---|---|
| C. tenellus WGS (Genoscope GCA_964200185.2, contig CAXNAD020000071.1) | 88 / 6 | 27 / 2 | 3 / 245 | 2 / 13 | ATA=Met |
| C. tenellus NC_085326 (Choi et al.; positive control) | 88 / 5 | 27 / 3 | 3 / 242 | 2 / 13 | ATA=Met |
| Discosporangium mesarthrocarpum WGS (GCA_964200235.2, CAXNAE020000396.1) | 3 / 51 | 113 / 4 | 6 / 190 | 0 / 10 | none |
| Phaeostrophion ON552993 (neg) | 2 / 69 | 98 / 12 | 5 / 126 | 0 / 28 | none |
| Protohalopteris MZ156064 (neg) | 3 / 55 | 93 / 3 | 3 / 166 | 0 / 17 | none |
| Saccharina NC_013473 (neg, Laminariales ref dropped) | 1 / 116 | 137 / 6 | 2 / 233 | 0 / 25 | none |
| Petalonia NC_025227 (neg, Ectocarpales ref dropped) | 0 / 127 | 135 / 2 | 3 / 214 | 1 / 39 | none |
| Sargassum NC_024613 (neg, Fucales ref dropped) | 4 / 118 | 131 / 8 | 1 / 202 | 1 / 32 | none |

- In Choristocarpus ATA is the main Met codon (88 of 115 ATA+ATG at Met positions) and has almost left Ile positions (6 vs 245 ATT);
  in every other brown alga tested ATA is a common Ile codon (51-127 at Ile) and rare at Met (0-4).
- Spread: Choristocarpus WGS has ATA at Met columns in 20 of 21 genes with any ATA/ATG Met site (Discosporangium 3 of 22,
  Phaeostrophion 2 of 22).
- Independence: the Genoscope assembly is a separate sequencing and assembly (99.9% identical over 32.3 kb), strain KU-1152
  (Kobe). The Choi record lists Kawai (Kobe) as an author and no strain, so it is probably the same culture: two submitters,
  one organism. Not two origins.
- Discosporangium, the only other genus in Discosporangiales, reads ATA as Ile: the change is on the Choristocarpus branch.
- Plastid of the same alga reads Ile (K1: 3/146), so this is mito-specific.

## Prior art (subagent, 09:24; checks/chor_wgs/priorart.md)
- NOT FOUND for Choristocarpus or Discosporangium. The Choi et al. mito paper is still unpublished. Choi 2024 (Curr Biol,
  plastids), Hwang 2026 (J Phycol, uses NC_085326) and Phaeoexplorer (Denoeud 2024 bioRxiv) don't mention the code.
  Codetta's benchmark predates the record (public Feb 2024). NCBI taxonomy gives mito code 1 for the genus.
- Precedents for AUA = Met in stramenopile mitochondria: Xanthophyceae (Ehara 1997; a related PX-clade class; the sister class is Schizocladiophyceae, corrected 14 Nov),
  Labyrinthulea and MAST8b (Žihala 2020). Ehara 1998 sampled only later-branching brown algae (no Discosporangiales) and
  judged a shared ancestral AUA = Met "highly unlikely". So the known xanthophyte state and this one would be separate origins
  unless the early-branching brown algae carry it. Not tested.
- Strain: both C. tenellus data sets come from culture KU-1152 (Villefranche 1980), first labelled Discosporangium and
  later re-identified as Choristocarpus (BioSample). So: one culture, two sequencings. The Discosporangium used here is a
  different strain (MT17-79, Malta), and its mito is only 76% identical, so it is a different taxon.

## Context run (09:26; post hoc, same b_read.py; b_read_context.out)
| mito | ATA at M / I | ATG at M | call |
|---|---|---|---|
| Sphaerosorus PV826132 (Xanthophyceae; AUA=Met published, Ehara 1997 lineage) | 81 / 6 | 16 | ATA=Met |
| Schizocladia NC_053869 (Schizocladiophyceae, sister class of brown algae) | 3 / 36 | 115 | none |
| Nematochrysopsis QB113265 (Chrysoparadoxophyceae) | 7 / 56 | 80 | none |
- The reader recovers the published xanthophyte code blind, with the same pattern as Choristocarpus (ATA has replaced most
  ATG at Met). Schizocladia, Discosporangium and the other brown algae read Ile, so on the accepted tree the Choristocarpus
  change is separate from the xanthophyte one.
- Contaminant check (tblastn of 8 proteins, length-weighted aa identity to the C. tenellus WGS mito): Sphaerosorus 63.9,
  Schizocladia 73.2, Fucus 72.8, Dictyota 73.1, NC_085326 99.6. So it is not a xanthophyte mito. Against the Discosporangium
  mito: C. tenellus 78.0, Fucus 77.2, Schizocladia 76.6. Choristocarpus is a long branch, about equally far from all
  of them. Identity does not place it inside Discosporangiales; a tree would. Not done.

## RNA-editing kill test (09:48; Amendment C/D/D2; c_rna.py, c_rna.out, c_rna_sites.tsv)
- ERR12954301 read 1, first 400 MB (4,466,748 reads): only 484 reads map to the mito (a poly-A library, probably).
- By the fixed rule: NO CALL. Only 6 ATA-at-Met sites reach depth >= 10; the rule needs 10. Those 6 have median G share 0.0
  and max 0.0.
- Pooled, descriptive only (not the rule): at the 71 covered ATA-at-Met sites, reads show A 270, G 0, T 1. At ATG-at-Met control
  sites: G 94, T 2. At ATA-at-Ile: A 11. If half the AUA codons were edited to AUG, we would expect about 135 G.
  So in these reads the Met-position AUA codons are transcribed and not edited.

## Kill check (10:12; kill/REPORT.md): SURVIVES, nothing killed
- An independent reimplementation reproduces the counts exactly. An annotation-based pipeline without tblastn gives 83/5.
- Not N-terminal. 20 of 21 genes. Conserved Met columns (Saccharina ATG at 88/88 of them).
- Non-brown stramenopiles (62 genera): at Choristocarpus ATA, Met share 0.390 and Ile 0.134, like ATG (0.437/0.093). ATT is 0.020/0.505.
- The substitution null (ATA share at Met columns should equal its Ile-column share, 0.023) is rejected: 88 of 93, P = 1e-137.
- RNA: in the kill run's 13.7M reads, ATA@M third base is A 796, G 0. My own 400 MB run gives A 270, G 0. No editing.
- tRNAs: Choristocarpus and Discosporangium both lack the canonical trnI(CAU) slot. Discosporangium keeps a divergent extra
  CAU tRNA ("X") that Choristocarpus lacks. That fits the claim, but charging is unknown.
- Corrections: "independent assembly" means two sequencings of one strain (KU-1152 = RCC7273); the trnI(CAU) loss predates
  the split.
STATUS: open lead, survives kill, prior art NOT FOUND. One strain, one genus, one origin (a separate origin from the xanthophytes).

## Second isolate? (10:17 14 Nov, Lyra)
- ENA read runs under Choristocarpus (tax_tree 116064): 3, all strain KU-1152 (Genoscope WGS ERR12682627/8, RNA ERR12954301).
- Non-WGS nuccore mito records: ON552980 / NC_085326 (genome) and AB899175 (cox1 partial, Kawai et al. 2015), which is also
  KU-MACC:KU-1152. So no public sequence from a second isolate exists; a second isolate needs new material.
