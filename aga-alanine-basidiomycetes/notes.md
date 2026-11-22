# Notes: fungi whose mitochondria read AGA as alanine

The agents' full write-up (method, evidence, checks, open points). The short version is in [`README.md`](README.md).
Times are US Eastern, November 2026. File paths are in this folder.

**Summary.** In the mitochondria of about 30 basidiomycete fungi (32 genomes, about 29 species), the codon AGA most
likely means alanine, not arginine. No known fungal genetic code reads AGA as alanine. **The protein evidence is one
peptide at one AGA site, in a protein with no known relatives.**

## Background

Genes are read three letters (one codon) at a time. In the standard genetic code, 61 codons each name one of 20 amino
acids and 3 are stops. Some lineages, many of them mitochondria, read a few codons differently; NCBI lists 26
alternatives to the standard code. AGA is arginine in every code used by fungi, and NCBI Taxonomy gives these fungi the
mold mitochondrial code (table 4, AGA = arginine).

## Method

- **Census.** About 22,700 mitochondrial genomes (18,015 of the 18,031 in RefSeq, plus one per additional species from
  GenBank) were translated with their declared codes, every protein was aligned to profile models, and the codon at
  every match column was recorded. A column is *strict* when at least 90% of the aligned genomes, each weighted by
  1/(genomes in its family), have one amino acid there (`src/codes_census/`).
- **Codon contrast** over all 4,638 fungal genomes in the census (`src/codes_census/codon_contrast.py`): per genome,
  AGA and its arginine synonyms (CGN, AGG) at strict alanine and strict arginine columns, one-sided Fisher test.
- **tRNA scan** (ARAGORN via pyaragorn) of all 4,638 genomes (`src/codes_census/org_trnas.py`), joined to AGA use by
  `src/aga_ala/trna_vs_usage.py`. A genome is in the alanine class if AGA sits at at least 1 strict alanine column and
  0 strict arginine columns, and in the arginine class if at 3 or more arginine columns and 0 alanine columns.

## Finding

The 32 genomes of the alanine class are yeasts and smut fungi on three branches (subphyla) of the basidiomycete tree:
*Trichosporon* and relatives (including the human pathogen *T. asahii*) and *Tremella* on one (Tremellomycetes,
Agaricomycotina), *Tilletia* (wheat bunt and smut fungi) and *Salmacisia* on another (Tilletiales,
Ustilaginomycotina), *Cystobasidium* on the third (Cystobasidiales, Pucciniomycotina).

| lineage | genomes | AGA at strict Ala | AGA at strict Arg | G3:U70 tRNA(UCU) |
|---|---|---|---|---|
| Trichosporonales | 18 | 67 | 0 | 17 |
| Tremellales | 2 | 13 | 0 | 1 (+1 at ARAGORN threshold 0.8) |
| Tilletiales | 7 | 34 | 0 | 6 |
| Cystobasidiales | 4 | 9 | 0 | 4 |
| Sporidiobolales | 1 | 1 | 0 | 0 |

(`data/census/trna_vs_usage_fungi.tsv`.) Two more genomes have AGA at both: *Tilletia* BK074393.1 (5 Ala / 1 Arg) and
*Cystobasidium* MT366950.1 (3 / 3); both carry the G3:U70 tRNA. MT366950.1 is the same genome as BK072865.1, one of the
32 (rotated, 99.6% identical); its 3 AGA-at-Arg codons are CGA in BK072865.1, so they are assembly differences.

- **Per genome:** AGA at 1 to 9 strict alanine positions and at none where they have arginine (that is how the 32 were
  picked). 124 such alanine sites in total, at 24 distinct positions. Six of the 32 rest on a single site. At strict
  columns AGA sits at alanine 78% of the time (124 of 159; the rest G 15, S 7, P 5, L 4, V 2, T 1), while GCN in the
  same 32 genomes sits at alanine 97.5% of the time (1,841 of 1,889) (audit and pre-post check,
  `data/checks/`).
- **Against its synonyms, selection-free** (`src/aga_ala/lineage_nosel.py`, `data/census/lineage_nosel.txt`): every
  census genome in each order whatever its AGA use, each strict column counted once per order. AGA at Ala / AGA at Arg
  / synonyms at Ala / synonyms at Arg: Trichosporonales 16 / 0 / 0 / 53 (22 genomes, p 5.2e-16), Tremellales 11 / 0 /
  0 / 55 (32, 9.3e-13), Tilletiales 6 / 1 / 0 / 51 (8, 1.7e-7), Cystobasidiales 5 / 3 / 0 / 50 (5, 1.2e-5). The
  pre-post checker's independent count is identical. Sporidiobolales (1 / 0 / 0 / 51, 21 genomes, p 0.019) rests on
  one site and was dropped from the post. The pooled per-class contrasts (`data/census/aga_ala_*.tsv`) give far smaller
  p-values, but they pool near-identical genomes; the selection-free version is the one to quote.
- **All aligned columns, no strict filter** (`src/aga_ala/all_columns.py`, `data/census/all_columns.txt`). In the 32
  genomes, the column consensus under AGA is A 30%, S 15%, L 8%, G 8%, T 7% and R only 1.0% (13 of 1,321 sites);
  under GCN it is A 53%, S 17%, G 6%, T 5%; under CGN it is R 74%, K 14%. In the 5 same-class genomes that put AGA at 3
  or more strict arginine columns (4 Exobasidiomycetes and *Solicoccozyma*), AGA sits under R 76% (59 sites), like CGN
  (71%). At only the non-strict columns, none of which was used to pick the genomes, AGA sits at A 23.1% and R 1.1%
  (1,162 codons; pre-post check, `src/prepost_check/`).
- **tRNA:** G3:U70 tRNA(UCU) in 28 of the 32, all with A73, the alanine discriminator. Of the 4,387 genomes with AGA at
  3 or more strict arginine columns and 0 alanine columns, 4,322 have a UCU tRNA (22 have no ARAGORN output, 65 no UCU
  tRNA), and none of those 4,322 has the G3:U70 form (`data/trna/fungi_ucu.aragorn.tsv`). Of the 4 of the 32 without
  one: *Hannaella* NC_039747.1 has one at threshold 0.8 (G3:U70, A73); in *Salmacisia* and *Pascua* the tRNA is
  found by blastn or has its anticodon mis-called; *Sporobolomyces* QB252587.1 rests on 1 site. *Auricularia*
  NC_058230.1's UCU tRNA is identical over 74 nt to *Tremella fuciformis*'s, likely contamination, so it is not counted
  as a carrier. Each carrier's tRNA(UCU) compared to its own other tRNAs: `data/trna/trna_origin_basidio.tsv`.
- **AGG**, the other codon a UCU tRNA can read, is almost absent: 12 aligned AGG codons across the 32, none at a strict
  column (`data/census/agg_ala.txt`). Its meaning can't be called.
- **Editing can't make it:** AGA to any GCN needs two base changes.

## Checks

- **Blinded re-check** (separate agent, its own code, Codetta v2.0 plus an ortholog-column tally against 11 table-4
  relatives; told only to infer each code against table 4, all 64 codons tested; its folder name did mention AGA):
  AGA = alanine in all 4 test genomes (*Trichosporon inkin*, *Tremella mesenterica*, *Tilletia indica*, *Apiotrichum
  gamsii*), no AGA = alanine call in the 4 controls. It also found AUA = Met in *Ceraceosorus* unprompted.
  (`data/checks/blind_recheck/`, `src/blind_recheck/`)
- **Audit** (fresh agent, `data/checks/audit/report.md`, `src/audit/`): no code bug makes or inflates the signal. Its
  own reproduction without the census or profile models (tblastn of reference proteins, code 4) gives AGA at alanine
  71% vs GCN 93%; AGA at arginine once in 6 test genomes (an artefact: the hit ran into an intron, the real codon is
  CGA) and 17 of 18 times in 3 controls; pooled p 4e-104. AGA is the only codon in the census whose column majority is
  not its standard residue.
- **Raw reads** (each assembly's source run, mapping quality 0 allowed; `data/reads/`,
  `src/codes_census/reads_codons.py`). At sites where relatives most often have alanine: *T. mesenterica* BK074411.1
  (SRR11469946) 30 sites, 164 of 164 reads; *Trichosporon asteroides* BK074449.1 (SRR1521503) 25 sites, 99.8%;
  *Tilletia horrida* BK074391.1 (SRR28015936) 25 sites, 99.7%; *Salmacisia buchloeana* BK074121.1 (SRR24299923) 26
  sites, 100%. No *Cystobasidium* reads were tested.
- **Pre-post check** (fresh agent, 11 Nov; `data/checks/aga-ala-prepost-check.md`): re-derived the census, tRNA, read
  and peptide numbers with its own code (`src/prepost_check/`). Verdict FIX, THEN POST; its fixes are applied.
- **Check of the lines added after** (fresh agent, 12 Nov; `data/checks/postcheck-additions-2026-11-12.md`): every
  number reproduces.
- **Tweet check** (21 Nov; `data/checks/tweet_check.md`).

## Protein test

PXD025163 (iProX): a TMT proteome of two *T. asahii* strains isolated 15 years apart from one patient (Q Exactive;
three 4-plex sets A, B, C of 16 fractions each). Searched with Sage against the CBS 2479 proteins (UniParc
UP000002748, 8,310) plus every stop-to-stop ORF of at least 50 codons of three *T. asahii* mitochondrial genomes, each
in three versions: table 4 (AGA = Arg), AGA = Ala, and a control with CGN = Ala (`src/prot/build_db.py`,
`run_sage.sh`, `tally_aga.py`). Peptide-spectrum matches at spectrum q ≤ 0.01 that tell the versions apart
(`data/proteome/set*/tally.txt`):

| Needs | Set A | Set B | Set C |
|---|---|---|---|
| CGN = Arg (positive control) | 20 (cox2) | 24 (21 cox2, 3 atp6) | 21 (20 cox2, 1 atp6) |
| AGA = Ala | 2 | 2 (+1 decoy) | 0 |
| AGA = Arg | 0 | 0 | 0 |
| CGN = Ala (negative control) | 0 | 0 | 0 |

- The peptide is TPAPNPDHVPEDHTR, from A1Q1_00013 (EJT45009.1, a 558-aa hypothetical protein of the CM001777.1
  mitochondrial genome, no relative in the profiles). Context EEPR|TPAPNPDHVPEDHTR|APKK: the flanking arginine codons
  are CGN, the peptide's alanine is AGA. Set A: fractions A14 and A13, hyperscore 58.0 and 41.3, spectrum q 8.8e-5
  and 8.3e-4; b1, b2, b3 and b5 matched in both. Set B: fraction B15, hyperscore 57.9 and 56.1, spectrum q 3e-5 and
  4e-4. (`data/proteome/setA/`, `setB/`)
- The arginine version TPRPNPDHVPEDHTR never matched.
- Genome control (`src/prot/pep_in_genome.py`, `data/proteome/setA/peptide_in_6_assemblies.tsv`): in 5 of 6 *T.
  asahii* assemblies (4 strains) the peptide region lies only on mitochondrial contigs, always with AGA at the alanine.
- One decoy (the reversed peptide SYYAIIK at q 0.0096, just inside the cut) also fell in the AGA = Ala class in set B,
  so that class's own error rate is not tiny; the peptide's spectra score far above it.
- A peptide over atp6's AGA site shows up only below the 1% cut-off, in both its alanine and arginine versions; it
  counts for neither.
- Sets B and C are more runs of the same two strains: they replicate the measurement, not the finding in other
  species.

## Known before

- No NCBI code (tables 1–33) reads AGA or AGG as alanine; NCBI Taxonomy gives these fungi table 4. Papers on these
  genomes use table 4 and call the UCU tRNA tRNA-Arg (e.g. Li et al. 2023, *IMA Fungus*, *Apiotrichum*). Ahrendt et
  al. 2026 (*Nucleic Acids Res.*; 328 fungal mitochondrial genomes, tRNAs named from the anticodon) report a
  tRNA-Arg(UCU) in all fungi except chytrids.
- AGA and AGG read as alanine in the mitochondria of the pedinophyte green algae Marsupiomonadales, through a
  tRNA(UCU) with G3:U70 and A73 (Turmel et al. 2020, *Mitochondrial DNA B*, first report; Barcytė et al. 2025, *PLoS
  Genet.*, order-wide). Barcytė et al. also report AGA and AGG as probably alanine (glycine not excluded) in the
  plastids of dinoflagellates whose plastids came from a pedinophyte. AGG (sometimes AGA) reads as alanine in
  Sphaeropleales (Žihala & Eliáš 2019, *Genome Biol. Evol.*), in some through a tRNA(CCU) with the alanine mark
  (Noutahi et al. 2019, *Mol. Biol. Evol.*). The only alanine reassignment known in fungal mitochondria is CUU/CUA in
  *Ashbya gossypii* (Ling et al. 2014). So this is a new fungal case of a known route.

## Open points

- **Proteins at the conserved sites.** One peptide at one AGA site, in an unconserved ORF; TMT spectra can carry
  co-isolated peptides. A second peptide, ideally at a conserved site, would make it solid.
- **How many origins.** The three groups sit in three subphyla, but they share 9 AGA-at-alanine columns against about
  4 expected by chance (P ≈ 7e-4; audit), and those columns are no less conserved than other alanine columns. Shared
  ancestry and codon hotspots are both open; no tree was built.
- **Why AGA.** Most *Tremella* and *Trichosporon* relatives use no AGA at all, so the codon may have been lost and then
  reused as alanine.
- **Novelty rests on searches.** The green algal precedent uses the same route, so a fungal case could have been
  noticed and not published.
