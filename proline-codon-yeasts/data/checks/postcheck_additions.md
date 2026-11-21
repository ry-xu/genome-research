# Post-check of the passages added after the 2026-11-11 pre-post checks

Independent fact-check, 2026-11-12, written at 12:20. Scope: the three passages Clover added or changed after the ~21:00
pre-post checks. Each sentence was checked against the project's saved result files. I recomputed only where no saved
file gives the number. I did not edit the drafts or any project file and downloaded nothing.

Scratch scripts are in `/private/tmp/claude-501/-Users-ryan-lab/c532536a-a44b-4f0d-84b9-39c5fcb71590/scratchpad/pc/`
(session scratch; may be cleaned). They are:

- `iso_mm2.py`: the checker's `isolates_mm2.py` plus a per-isolate printout;
- `swaps_all.py`: every UGA and Arg codon of NC_020331 in the 18 isolates;
- `nad1col.py`: unweighted and weighted Arg share at one census column;
- `ds02_maincols.py` and `acc_maincols.py`: DS02's CUN columns against the main census consensus;
- `ds02_link.py` and `ds02_sites.py`: DS02 reads against the assembly and BK074532;
- `ds02_vs_bk.paf` and `ds02_vs_wick.paf`: minimap2 asm20 alignments.

Labels:
- **OK**: the saved files support the claim as written.
- **LOOSE**: true in substance, but the wording overstates, understates or could mislead.
- **UNSUPPORTED**: no saved file supports it, and I could not check it offline.
- **WRONG**: the saved files contradict it.

---

## 1. `outbox/2026-11-11-mito-aga-alanine/post.md`

### 1a. "**Proteins:**" bullet (line 13)

| Sentence or claim | Source file(s) | What I found | Label |
|---|---|---|---|
| "a public proteome of two *T. asahii* isolates from one patient (iProX, PXD025163; three sample sets, all searched)" | `results/aga_ala/prot/set{A,B,C}/results.json`, `set{B,C}/stream.log`, README 250–252 | Sets A, B and C were each searched, 16 mzML files apiece. "Two isolates from one patient" was already checked in the pre-post check (OmicsDI). That the dataset has only these three sets comes from the README's reading of the ProteomeXchange record. No file list is saved (`src/prot/stream_tmt2.sh` builds the URLs from a pattern). Low risk. | OK |
| "one peptide matched with alanine at an AGA codon: TPAPNPDHVPEDHTR" | `set{A,B,C}/tally.txt`, `setA/ala_peptide_psms_A13A14.tsv`, `setB/ala_peptide_psms.tsv`, `setB/mito_psms_q01.tsv` | The AGA = Ala class at q ≤ 0.01 has 2 targets in set A, 2 targets plus 1 decoy in set B, and 0 in set C. Every target is TPAPNPDHVPEDHTR. | OK |
| "from a mitochondrial gene of unknown function" | README 264–265 (A1Q1_00013 / EJT45009.1, hypothetical, CM001777.1 mitogenome) | As stated. | OK |
| "in 4 spectra from 2 of the 3 sets (none in the third)" | as above; `setC/target_peptides_anyq.tsv` (header only) | The 4 spectra are A14 scan 3271, A13 scan 3307, B15 scan 3203 and B15 scan 3184. Set C has no hit for either version at any q. Note: within each set the two spectra are one elution event (set A: adjacent fractions at RT 15.11/15.23 min; set B: one fraction, charge 3 and 4 at RT 15.15/15.10). So this is 2 independent observations, one per set. The sentence does not claim more. | OK |
| "the arginines on either side of it are CGN codons" | `results/aga_ala/prot/sites.tsv`, `db.fasta` (ORF CM001777.1_+1_24946) | The R before the peptide is CGT (site 371). The peptide's Ala is AGA (374). The R that ends the peptide is CGA (386) and is the peptide's own last residue; the residue after the peptide is A (…HTR\|APKK). So both cleavage-site arginines are CGN, as intended, but "either side of it" puts the second one outside the peptide. | LOOSE (minor) |
| "Its arginine version never matched." | `set{A,B,C}/target_peptides_anyq.tsv` | TPRPNPDHVPEDHTR is absent at any q in all three sets. The only related hit is the missed-cleavage TPRPNPDHVPEDHTRAPK in set A, at q 0.228 with isotope error +1. | OK |
| "65 spectra matched peptides that need CGN = arginine (from the cox2 and atp6 proteins)" | `set{A,B,C}/tally.txt` | 20 + 24 + 21 = 65, with no decoys. That is 61 cox2 (ORF CM001777.1_+2_24401) and 4 atp6 (AVSLGVR, ORF CM001777.1_+1_23092: 3 in B, 1 in C). | OK |
| "none matched peptides that need CGN = alanine" | same | The C class is 0 in all three sets. The RC and R classes are 0 too. | OK |
| "The peptide's sequence exists only on mitochondrial contigs of 5 *T. asahii* assemblies, always with AGA at that alanine." | `setA/peptide_in_6_assemblies.tsv`, README 271–273 | 6 assemblies were searched. 5 have the core PRPNPDHVPEDH (6 hits), all CCT **AGA** CCA. CBS 17435 has no copy at all. The 5 assemblies are 4 strains (CBS 2479 twice: CP116791.1 and CM197668.1). Two of the contigs are not labelled mitochondrial and are mitochondrial by size and content: M6011 ctg_44 (32,648 bp), and Z1 Scaffold11 (62,836 bp, two copies 31,420 nt apart, so a mitogenome dimer). The sentence reads as if it covers all assemblies; it leaves out the one that lacks the sequence. | LOOSE |
| "One peptide at one site, in a protein with no known relatives, is thin." | — | Judgement; consistent with the data. | OK |

### 1b. Note bullet "**Proteome sets B and C**" (line 37)

| Sentence or claim | Source file(s) | What I found | Label |
|---|---|---|---|
| "searched 2026-11-11 22:00 and 2026-11-12 00:16, after the pre-post check" | `setB/stream.log` (ends Wed 22:01:06), `setC/stream.log` (ends Thu 00:16:52) | These are the completion times. | OK |
| "set B matched the same peptide in 2 more spectra (one fraction; hyperscore 57.9 and 56.1, spectrum q 3e-5 and 4e-4" | `setB/ala_peptide_psms.tsv` | B15 scans 3203 and 3184: hyperscore 57.85 and 56.12, spectrum q 3.16e-5 and 3.51e-4, charge 3 and 4. | OK |
| "its TMT reporter ions sit in the same 4 channels as cox2's, as in set A" | `setB/reporters.txt`, `setA/reporters.txt` | True: 126, 129N, 129C and 131 dominate for the peptide and for cox2 in both sets. But this is a 4-plex, so those are the only labelled channels; every TMT-labelled peptide, nuclear or mitochondrial, shows the same 4. The channel pattern shows the spectrum is a real labelled one. It says nothing about mitochondrial origin, which the parallel with cox2 implies (README 276–277 says it outright: "So it comes from the same samples as a known mitochondrial protein"). | LOOSE |
| "set C matched neither version" | `setC/target_peptides_anyq.tsv` (empty), `setC/tally.txt` | As stated. | OK |
| "No new AGA = Ala or AGA = Arg peptide in either" | `set{B,C}/tally.txt` | The only AGA = Ala target is the same peptide. AGA = Arg (RC class) is 0. | OK |
| "CGN = Arg controls 24 and 21 spectra" | same | Set B has 24 (21 cox2, 3 atp6); set C has 21 (20 cox2, 1 atp6). | OK |
| "One decoy (a reversed peptide at q 0.0096, just inside the cut) also fell in the AGA = Ala class in set B" | `setB/mito_psms_q01.tsv` | The decoy is rev SYYAIIK, ORF CM001777.1_+0_8661, hyperscore 33.45, q 0.009577. | OK |
| "the peptide's spectra themselves score far above it" | the PSM files above | In q, all 4 spectra are ≤ 3.9e-4, against the decoy's 0.0096, so at least 25× lower. In hyperscore the set B spectra are 57.9 and 56.1 against 33.4. The weakest spectrum overall (A13, hyperscore 41.3) is closer, but the sentence is about set B. | OK |
| "Sets B and C are more runs of the same two isolates, so this is a replication of the measurement, not of the biology in other species" | README 250–252, 299 | Consistent. | OK |
| "Proteins line updated; nothing else in the post changed." | file times | No earlier copy of the text exists to diff against. `image.png` is unchanged (11-11 13:39). | OK |

### Replacement text, AGA post

1. **(LOOSE, recommended)** Line 13. Replace
   "The peptide's sequence exists only on mitochondrial contigs of 5 *T. asahii* assemblies, always with AGA at that alanine."
   with
   "In the 5 of 6 searched *T. asahii* genome assemblies that contain it (4 strains), the peptide's sequence lies only on mitochondrial contigs, always with AGA at that alanine."
2. **(LOOSE, optional)** Line 13. Replace
   "the arginines on either side of it are CGN codons."
   with
   "the arginines at both of its cut sites (the one before it and its own last residue) are CGN codons."
3. **(LOOSE, recommended)** Line 37. Replace
   "(one fraction; hyperscore 57.9 and 56.1, spectrum q 3e-5 and 4e-4; its TMT reporter ions sit in the same 4 channels as cox2's, as in set A)"
   with
   "(one fraction; hyperscore 57.9 and 56.1, spectrum q 3e-5 and 4e-4)".
   If you keep a reporter clause, make it "its TMT reporter ions fill the 4 labelled channels, as every labelled peptide's do in this 4-plex".

---

## 2. `outbox/2026-11-11-mito-proline-codons/post.md`

### 2a. "**Two more sequencings:**" bullet (line 9)

| Sentence or claim | Source file(s) | What I found | Label |
|---|---|---|---|
| *W. siamensis* PYCC 7069 reads "carry CUU/CUA at all 82 sites (1,432 of 1,450 reads)" | `results/cun_pro/audit/iso_census.out` | 82 census sites, all covered. CTA 70 sites with 1,234 of 1,248 reads; CTT 12 sites with 198 of 202. Total 1,432 of 1,450 (98.8%). The lowest site is 8 of 9 (cox1 Vir:247). | OK |
| "PacBio reads of a second strain, DS02 (2017; …)" | README 138–140; `results/cun_pro/audit/factcheck.md` 166; pre-post check line 100 | Both earlier checks read the ENA record (PRJNA318722, public 2017). The data back "a second strain": DS02's contigs differ from BK074532 (NRRL Y-17858) at about 3% of aligned bases (see 2b) and at 5 codon sites where the reads and the assembly agree with each other. | OK |
| "filed as *W. infanticola* with its reads" | README 138; earlier checks as above | Not in a saved file. Two earlier checks read the ENA record and treat it as *W. infanticola*. | OK |
| "and as *W. sorbophila* with its … genome assembly" | `data/wgsmito/picks_b5.tsv`, `b5.records.tsv` | GCA_002251995.2, taxid 45607, "Wickerhamiella sorbophila". | OK |
| "its 2018 genome assembly" | — | The year is in no saved file. `picks_b5.tsv` keeps only the FTP path (ASM225199v2), and `b5.gb.gz` has synthetic headers dated 1980. | UNSUPPORTED (minor) |
| "the ATCC lists *Candida infanticola* as a synonym of *W. sorbophila*" | — | This is in no saved file and appears only in README 147 and the Note. I could not check it offline, and it is the one clause that explains the naming conflict. | UNSUPPORTED |
| "carry CUU or CUA as the commonest triplet at all 65 sites with 3 or more reads" | `reads_isolates/SRR4787334_pb.BK074532.tsv`, `results/cun_pro/audit/BK074532.1.census_ctn_p.txt` | My own count: 78 sites; 67 have reads (1,942 reads, 1,540 = 79.3% showing the assembly codon). Of those, 65 have ≥ 3 reads, and all 65 have CTT or CTA as the strict commonest triplet. There are no ties, and the unlisted reads cannot overturn any of them. 61 of the 65 have the assembly codon as majority. | OK |
| "That assembly's mitochondrial contigs" (i.e. DS02's own assembly) | my `ds02_link.py` on the reads, GCA_002251995.2 contigs and BK074532 | There are 5 sites where the DS02 PacBio majority differs from BK074532: Vir_cox1:221 TTA, nad4:399 TTG, cox3:56 ATA, cox3:114 TTG and nad5:274 CTT. The assembly matches the reads at all 5, against 0 that match BK074532, and agrees at 543 further sites. So the reads and the assembly are the same strain. | OK |
| "which a later scan of nuclear genome assemblies reached on its own" | README 1416–1432; `src/wgsmito/run_batch.sh` | The b5 scan covers all whole-genome (WGS) assemblies not in b1–b4 (997 genomes) and takes their mitochondrial contigs, so "on its own" holds. These are whole-genome assemblies, not nuclear ones; the scan found mitochondrial contigs in them. | LOOSE |
| "put CUU or CUA 78 times at conserved proline positions and never at conserved leucine" | `results/profile_scan/contrast_wgs_b5.tsv` | CTT: a = 17 at P, 0 at L, n_strict 17, p 4.1e-8. CTA: a = 61, 0 at L, n_strict 61, p 8.1e-36. 17 + 61 = 78, at 78 distinct model columns, so the about 5-fold contig redundancy (308 kb of contigs for a 56 kb genome) did not inflate the count. Note: the b5 contrast scan takes column consensus from the b5 batch alone. Against the main-census consensus used everywhere else in the post (mito + gbmito + gbmitoorf, family-weighted, 0.9; checker's `bg.cols.npz`), DS02 has **85** (CTA 66, CTT 19) at strict proline columns, 0 at strict leucine, and 3 at other strict columns (2 Phe, 1 Met). My method gives 78 for BK074532, matching the post. So 78 is conservative, and its match to BK074532's 78 is a coincidence of definitions. "Never at conserved leucine" is true under both. The pre-post check called it "true but not evidence", and line 5 already carries that context. | OK |

### 2b. Note bullet "**DS02 naming**" (line 37)

| Sentence or claim | Source file(s) | What I found | Label |
|---|---|---|---|
| "DS02's 2018 assembly (GCA_002251995.2, KRIBB, PRJNA383455) is filed as *Wickerhamiella sorbophila*" | `data/wgsmito/picks_b5.tsv`, `b5.records.tsv` | The accession and species name are saved. The year 2018, KRIBB and PRJNA383455 are in no saved file (README only). | UNSUPPORTED (the year, submitter and BioProject only) |
| "its PacBio reads (PRJNA318722) as *W. infanticola*" | README 138; earlier checks | As in 2a. | OK |
| "ATCC 201079 lists *Candida infanticola* as a synonym of *W. sorbophila*" | — | No saved source. | UNSUPPORTED |
| "My WGS mito-contig scan (batch b5) picked the assembly up independently: CTT 17 and CTA 61 at strict proline columns, 0 at leucine" | `results/profile_scan/contrast_wgs_b5.tsv`, README 1416–1432 | As in 2a. "Strict" here means by the b5 batch's own consensus; by the main census it is 85 (CTA 66, CTT 19). | OK |
| "Its mito contigs are 94–97% identical to BK074532 over 7–8 kb blocks" | my `ds02_vs_bk.paf` (minimap2 asm20, DS02 contigs onto BK074532) | All 15 contigs align only to BK074532, over about 280–300 kb of contig sequence (the contigs overlap heavily). Outside gaps (1 − de), identity is 96.8% overall and 95.3–97.7% in every block ≥ 5 kb; blocks run to 26 kb. Counting gaps (matches/block), the 7–8 kb blocks are 92.8–94.5%. The long blocks drop to 67–78% because of kb-scale insertions and deletions (e.g. 1,192 and 1,315 bp deletions, 969 and 812 bp insertions). "94–97% over 7–8 kb blocks" is roughly right for those particular blocks, but it hides both the long alignments and the large indels, and the two identity measures are mixed. | LOOSE |
| "so a different strain of the same or a sister species" | `ds02_vs_wick.paf`; the 20 *Wickerhamiella* in `data/gb_mito/gb_024.gb.gz` | BK074532 against its nearest other species in the census is 92.0% outside gaps for *W. verensis* BK074542 (27 kb) and 90.1% for *W. siamensis*. DS02 at about 97% is much closer than any other species, which fits a second strain of *W. infanticola* or a very close sister. | OK |
| "The post's Two more sequencings line now says this; nothing else changed." | file times | Not gradable: there is no earlier copy to diff. Observation only: `image.png` was rewritten at 12:00 with `post.md`. It shows the same 49-genome heat map with no DS02 row. If it was only re-rendered, ignore this; if its content changed, the Note should say so. | — |

### Replacement text, proline post

1. **(UNSUPPORTED, required)** Line 9. Either save the ATCC 201079 record (for example to `results/cun_pro/lit/`) and keep the clause, or cut it. In place of
   "; the ATCC lists *Candida infanticola* as a synonym of *W. sorbophila*)"
   a data-backed alternative is
   "; its mitochondrial contigs match the *W. infanticola* genome at about 97% identity, against 92% between *W. infanticola* and its closest other species)".
   Or simply end the parenthesis after "assembly)".
2. **(UNSUPPORTED, minor)** Line 9. Replace "with its 2018 genome assembly" with "with its genome assembly", unless you save the NCBI assembly report (it has the date and submitter).
3. **(LOOSE, required)** Line 9. Replace "which a later scan of nuclear genome assemblies reached on its own" with "which a later scan of whole-genome assemblies reached on its own".
4. **(LOOSE)** Line 37. Replace
   "Its mito contigs are 94–97% identical to BK074532 over 7–8 kb blocks, so a different strain of the same or a sister species."
   with
   "Its mito contigs match BK074532 at about 97% identity outside gaps (95–98% per block, blocks up to 26 kb, with some kb-scale insertions and deletions); BK074532 and its closest other species, *W. verensis*, are 92%. So DS02 is a different strain of the same or a sister species."
5. **(UNSUPPORTED)** Line 37. Either save the records for "2018", "KRIBB", "PRJNA383455" and "ATCC 201079 lists …", or drop them. Without the ATCC sentence the Note still stands on the accession, the species names and the identity figures.
6. *Optional:* if you want the count on the same definition as the rest of the post, say "85 times at conserved proline positions (78 by that scan's own batch consensus)". Keeping 78 is fine; it is the conservative figure.

---

## 3. `outbox/2026-11-11-mito-uga-arginine/post.md`, "**UGA and arginine codons swap like synonyms:**" (line 10)

| Sentence or claim | Source file(s) | What I found | Label |
|---|---|---|---|
| "inside *P. jirovecii*, UGA and regular arginine codons trade places at two sites" | `results/tga_arg/isolates.tsv`; my `swaps_all.py` (minimap2 asm20, NC_020331 against the 18 isolates in `data/tga_arg/others.fa`, JX499143.1 excluded) | I checked every UGA (42) and Arg codon (39) in the reference's 15 CDSs, all of them read in all 18 isolates. Reference UGA with an isolate Arg codon, or reference Arg codon with an isolate UGA, happens at exactly 2 sites: nad2:282 and nad1:33. So "two sites" is exhaustive, not a lower bound. | OK |
| "At one, the reference genome has UGA and 12 of 18 isolates have AGA." | `results/tga_arg/isolates.tsv`; `iso_mm2.py` re-run | nad2:282 is TGA in the reference. AGA in MW530527, MW530526, MW530525, MT319762, MH010446, MH010445, MH010444, MH010442, MH010441, MH010439, MH010438 and JX855938 (12). TGA in the other 6. | OK |
| "At the other, a position where 99% of genomes have arginine" | RNA tables (`results/tga_arg/rna/pjir_*.tsv`: column mito_Fun_nad1:40, share 0.99, strict 1); my `nad1col.py` on mito + gbmito + gbmitoorf census | The saved 0.99 is family-weighted (0.9948; the checker's `bg.cols.npz` agrees). Unweighted it is 4,185 of 4,195 sites = 99.76%, and 4,179 of 4,189 genomes all-Arg. The others are 7 K, 1 M, 1 L and 1 N. So "99% of genomes" holds under either reading. (Loading the census peaked at 0.62 GB resident.) | OK |
| "the reference genome has CGA" | RNA tables, codon column; `iso_mm2.py` | CGA in NC_020331 (and JX499143.1). | OK |
| "but 6 of 18 isolates have UGA" | `iso_mm2.py` re-run (the checker's `isolates_mm2.py`, no saved output until now) | TGA in MW530528, MW530527, MH010443, JX855938, JX855937 and JX855936. CGA in the other 12. This matches the Note ("CGA in 12 of 18 isolates and UGA in 6"). | OK |
| "and the patient's RNA reads UGA in 41 of 41 reads" | `results/tga_arg/rna/pjir_ERR135872–81.tsv`, nad1 33 row | Reads per run: 5, 5, 4, 6, 2, 3, 4, 7, 3, 2 = 41, all TGA, 0 CGA. | OK |

Side note, not a fix: the two swaps do not split the isolates the same way. MW530527 and JX855938 have AGA at nad2:282 and UGA at nad1:33; MH010440 and MH010437 have UGA at nad2:282 and CGA at nad1:33. That fits independent changes, and nothing in the bullet depends on it.

No replacement text needed.

---

## Side issues outside the three passages (README; not in the drafts)

- README 268: "posterior error probability 8.8e-5 and 8.3e-4". Those are the **spectrum_q** values in `setA/target_peptides_anyq.tsv`. Sage's `posterior_error` column is log10: −3.06 and −2.04, i.e. PEP 8.7e-4 and 9.2e-3.
- README 281: "Same database, config and tally". Set A was searched with `isotope_errors` [−1, 3]; sets B and C used [0, 0] (`set*/results.json`). Everything else is identical.
- README 298: "AVSLGVRLFANM(ox)VAGHTLLK at q 0.028" is in no saved file. `setA/target_peptides_anyq.tsv` has AVSLGVRLFANMVAGHTLLK (no oxidation) at q 0.396.
- README 276–277: "So it comes from the same samples as a known mitochondrial protein" is tautological in a 4-plex (see 1b).
- README 147 and 1428–1430: the ATCC 201079, KRIBB 2018 and PRJNA383455 statements have no saved record.

---

## Verdicts

- **AGA = alanine post: FIX FIRST.** Two quick wording fixes: the assembly sentence (5 of 6 assemblies, 4 strains) and the reporter-ion clause in the Note. Every number checks out.
- **CUU/CUA = proline post: FIX FIRST.** The ATCC synonym clause is unsupported by any saved file; cut it or save the record. "Nuclear genome assemblies" should be "whole-genome assemblies". Drop "2018" or save its source. In the Note, restate the identity figure. Every count checks out.
- **UGA = arginine post: POSTABLE AS IS.** Every claim in the bullet reproduces, including the 6-of-18 isolate count (re-run) and "99%" (both weighted and unweighted). "Two sites" is exhaustive.
