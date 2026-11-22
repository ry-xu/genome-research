# Pre-post check: *Bilabrum* reads TAG as glutamine (stop-sense-euk) (independent check, 2026-11-12)

Claim checked, from the outbox draft `outbox/2026-11-11-bilabrum-tag-gln/post.md` (updated 03:13) and its `image.png`:
- *Bilabrum* cf. *latius* (GCA_965178125.1, a Sanger Aquatic Symbiosis Genomics cobiont of *Dinothrix rugata*) reads the stop codon TAG as glutamine in its nuclear genes.
- "No change to the nuclear genetic code was known before in stramenopiles."

I read the post, `evidence.md`, the first fact-check (`factcheck/report.md`), the independent review outputs (`review/`), the 1KP RNA outputs (`rna_map/`), the screen README, `summary_all.txt`, `labelled_all.tsv`, `stage2/`, `triage_full.md` and `notes/known-nuclear-stop-codes.md`. I re-derived the post-body numbers from the genome FASTA, the tblastn tables and the three PAF files with my own scripts (`scratchpad/prepost/scripts/`: `rederive_sites.py`, `ends_check.py`, `erf1_s60.py`, `rna_codons.py`, `rna_rates.py`, `site_codons.py`, `near_tag_identity.py`). I ran the literature search myself. I edited nothing else. Written 12:16 EST; details re-verified against source files 12:22 EST.

## Verdict: POST. A few optional wording and figure fixes; nothing in the body is wrong.

- **Every number in the post body re-derives** from raw files (section 2). Where my calls differ, the cause is a different read-calling rule, and the difference is small (RNA: 12.5 / 15.6 / 1.0 per 1,000 against the post's 12.3 / 14.4 / 0.8).
- **Novelty holds** (section 1).
  - No paper, preprint, genome note or NCBI table reports a nuclear code change in any stramenopile, bicosoecid or *Bilabrum*.
  - The only earlier sighting is Vierock et al. 2022's supplement, which treats the TAG as a "premature stop codon". The post already says this.
  - NCBI's tables still end at 33, and every *Bilabrum* taxon is listed with code 1.
- **The alternative explanations are ruled out** as far as genome plus related-strain RNA can rule them out (section 2, end): host or cobiont contamination, assembly or base errors, pseudogenes, and RNA editing in the related strain.
  - Two things stay open, and the post's "Not known" says both: RNA or protein from this strain itself, and whether TAG still ends a few genes.
- **The figure is correct**: every codon, the boxed site and the "32 reads" all check. Two optional label fixes (section 4).
- **Pending, not errors**:
  - the SRA/Logan line in the Note;
  - *Maristentor*;
  - the *Pirsonia*/*Mutimo* gap runs.
  The Note's "Screen finished" and "the screen has finished" (L35, L43) should be squared with those lanes when the final numbers come in.

## 1. Novelty

I searched Europe PMC (REST; it covers bioRxiv and Research Square), the web, Wellcome Open Research, NCBI and the Vierock supplement. I searched by mechanism as well as by organism. Downloads are in `scratchpad/prepost/lit_epmc/`, `dl_vierock_supp/` and `dl_ncbi_gc/`.

- **"Bilabrum", all fields:** 8 records in total. None is about the genetic code.
  - Schoenle 2020 *Eur J Protistol* (the genus description) and Baričević 2023 *Eur J Protistol* (Adriatic record).
  - Channelrhodopsin papers:
    - Vierock 2022 *Sci Adv* (PMC9733931);
    - Govorunova 2022 *mBio* (PMC9765531);
    - Tajima/Zhang 2023 bioRxiv 10.1101/2023.07.22.550149;
    - *Nat Commun* 2024 39181878 (PMC11344840);
    - *Biophys Physicobiol* 2023 (PMC10865875).
  - *ISME J* 2024 stramenopile review (PMC11412368).
  - I grepped the five open full texts for "stop codon", "genetic code", "UAG/TAG" and "premature". The only hit is a mitochondrial-code reference in the ISME review.
- **Vierock et al. 2022, supplementary text S1** (`sciadv.add7729_sm.pdf`, the "Origin of B1ChR1 and B1ChR2" section). It says: "The original transcript 1KP_LLXJ_2042068 and some of other B1ChR2 variants appeared to contain premature stop codons… The haplotype with no premature stop codons was chosen".
  - This matches the post's quote and paraphrase (L16) word for word. It treats the TAG as an error to avoid, not a code. **Does not pre-empt.**
  - Two details the post doesn't use: Vierock call the contaminant "a non-clonal *Bilabrum* sp.", and found it in **four** 1KP cultures (LLXJ, JIWJ, ALZF and JQFK *Picochlorum atomus*). See the optional fix to L9.
- **Bicosoecids**, 2025-2026 (26 records), and (bicosoecid OR *Cafeteria*) AND (genetic code OR stop codon) (113 records). Nothing reports a code change. What exists:
  - the *Halocafeteria* genome preprint 10.64898/2026.06.25.734631;
  - *Cafeteria* spatial proteomics 2026;
  - *Hirugamonas* gen. nov. 2025 (morphology);
  - virophage papers.
- **Stramenopiles**: (stramenopile/heterokont/Gyrista/Bigyra/Placidozoa/Opalozoa) AND (genetic code OR codon reassignment), 154 records. Every code change found is **organellar**:
  - Labyrinthulea mitochondria (*Mol Phylogenet Evol* 2020, with a 2026 corrigendum);
  - *Plocamiomonas psychrophila* mitochondria, code 4, "a rare feature among Stramenopiles" (*Sci Rep* 2026, 10.1038/s41598-026-58736-z; new since the ledger);
  - eustigmatophyte organelles (*GBE* 2026).
  The 2026 *Pirsonia diadema* paper (*Nat Microbiol*) says nothing on the code. **None is nuclear.**
- **Stop-to-sense or UAG-glutamine, 2024-2026, any eukaryote** (90 and 268 records). Hits:
  - Maturová 2026 bioRxiv (10.64898/2026.11.04.764322, ciliates only);
  - McGowan 2024 *PLoS Genet* (phyllopharyngean UAG = Leu or Gln);
  - Záhonová 2025 *mBio* (*Blastocrithidia*);
  - *Palmarella* 2024;
  - the *Amoebophrya* genome 2025;
  - heterolobosean *Dactylomonas* 2025.
  All are in lineages the post's L14 already covers ("…and a few others"). **None is a stramenopile.**
- **Genome-wide code screens.** Codetta's 250,000-genome run (Shulgina & Eddy 2021) covered bacteria and archaea only. The ledger lists KACI (Melnykov 2026 *PNAS*), gTranslate 2026 and Fakih 2026 as prokaryote-only (`notes/known-nuclear-stop-codes.md` L123, judged from titles and abstracts). EukPhylo (*mBio* 2025) is a phylogenomics toolkit, not a code screen.
  - Swart 2016 is the only broad eukaryote screen, and it found no change outside ciliates (ledger, full text read). It included *Cafeteria*: MMETSP0942 and MMETSP1104 are *Cafeteria* in `data/units_all.tsv`.
- **Sanger ASG genome note for *Dinothrix rugata* or its cobionts:** none. Europe PMC finds no genome note naming *Dinothrix* or *Galeidinium*, and Wellcome Open Research has no article matching "Dinothrix" or "Bilabrum". The web search found only general ASG project pages.
- **NCBI genetic codes:**
  - `gc.prt` and the live `wprintgc.cgi` page both end at table 33 (Cephalodiscidae mitochondrial). There is no table 34, and no stramenopile nuclear table. Table 23, *Thraustochytrium* mitochondrial, is the only stramenopile entry.
  - Taxonomy efetch: *Bilabrum* (2721177), *B. latius* (2604063), *B.* cf. *latius* (3115832), *B.* sp. (3503367) and unclassified *Bilabrum* (3503366) all have GCId 1 and MGCId 1. This confirms the post's "NCBI lists the standard code for every *Bilabrum* taxon".
- **Not seen by me or by anyone in the project:** Lukeš et al. 2025 *Curr Biol* Table S1 (the publisher returned 403). The ledger says its Fig 4 has no stramenopile row. The risk is low, because a stramenopile nuclear code would have its own primary paper, and none exists.

**Bottom line:** "No change to the nuclear genetic code was known before in stramenopiles" holds today. The post's "Known before" list (L14-L16) is accurate.

## 2. Correctness: every number in the post body

| post claim (line) | source file | re-derived check |
|---|---|---|
| 1,785 genomes and transcriptomes screened (L3) | `results/summary_all.txt`, README 00:59 | holds: 1,927 scored, 1,785 non-duplicate (1,539 powered / 1,406 non-duplicate, as in the Note) |
| Rules set on 67 known-code units, frozen before the screen (L3) | README L60-64, `scan_controls.log`, `scan_all.log`, `scan_bico.log` | holds: 43 positive + 24 negative = 67; rules frozen 21:17; full scan started 21:20-21; bicosoecid scan 21:55-21:56. *Spironucleus* was relabelled at 21:16, after its score was seen but before the freeze (README says so) |
| GCA_965178125.1, Sanger, March 2025, 43 chromosomes (L3) | `factcheck/dl/web1/asm_bil.json` | holds: released 2025-03-04, submitter Wellcome Sanger Institute, 43 chromosomes |
| TAG at 166 aligned positions; TAA 1, TGA 4 (L5) | `review/bil_human76.tsv` + genome FASTA; my `rederive_sites.py` | holds exactly. The codon at every `*` was re-read from the genome; all TAG sites are on OZ chromosomes |
| No TAG in *Cafeteria* (L5) | `review/caf_human76.tsv` + *Cafeteria* FASTA | holds: TAG 0 (TGA 13, 7 aligned) |
| 91 cleanest sites: Q 44, K 11, no other > 5 (L6) | same; criteria ≥10 columns from the HSP end, no query gap within ±6, ≥8/25 identical | holds exactly: Q 44, K 11, E 5, S 5, R 5, N 4, … |
| "the pattern of CAG" (L6) | `out/rna_codons.tsv` (all deduplicated panel codons, best HSP) | holds. Human residue opposite all aligned codons: CAG → Q 40.9%, E 9.0%, K 8.8% (n = 2,633); CAA → Q 35.2%, K 13.0%, E 8.1% (n = 247); TAG → Q 33.7%, K 13.9%, R/E 6.6% (n = 166). Same rank order, Q then K/E |
| Ends none of 149 genes (TGA 136, TAA 13); *Cafeteria* 11 of 145 (L7) | `review/terminal_codons_relatives.tsv`; `new-stops-euk/.../scan/GCA_965178125.1.ends.tsv`; `ends_check.py` | holds for the table it came from (`scan/`, Nov 11 00:25): 149 loci, 136/13/0. The screen's later `scan2/` table (10:44) gives 163 loci, 147/16/0. TAG is 0 in both. *Cafeteria* 11/145 holds; one-sided Fisher p = 3.4e-4 |
| eRF1 valine where 2,714 of 2,776 have serine; none valine (L8) | `erf1_aln.fa`; `scratchpad/erf1/erf1_ref.tsv` (3,189 UniProt entries); `erf1_s60.py` | V = 0 holds under every method. The *Bilabrum* motif reads TA**V**NIKS (human TASNIKS, S60). The exact S count depends on how position 60 is located: 2,708 of 2,752 (motif regex) or 2,792 of 2,875 (local alignment to human 1-140), against the post's 2,714 of 2,776 (the 2,776 of 3,189 entries that carry the motif, per `evidence.md` L20). All give 97-98% serine, so the wording is fine |
| tRNA with anticodon CUA scores glutamine- or histidine-like, not tyrosine-like (L8) | `trna.out`, `review/iso.tbl` | holds: OZ237757.1 Sup-CTA, score 70.6. Isotype models: Gln 73.1, His 73.2, Glu 51.4, **Tyr 51.3**. A second Sup-CTA on OZ237737.1 scores 25.5 and is marked pseudo; the post doesn't need it |
| 1KP: TAG at 41 of 60 covered sites, CAG at 18 (L9) | the three PAFs; my `rna_codons.py` / `rna_rates.py` (primary, MAPQ ≥ 10, all three bases aligned; ≥ 2 reads, strict majority) | holds exactly: 41 TAG, 18 CAG, 1 CAA. Reads: TAG 646, CAG 363. Opposite human Q: 29 covered, 22 TAG |
| CAG→TAG 12.3 per 1,000 vs CTG→TTG 14.4 (L9) | same | mine: 24/1,923 = 12.5 vs 40/2,563 = 15.6 (two-sided Fisher p = 0.45). Per run: LLXJ 12.2 vs 15.4; ALZF 5.6 vs 10.7; JIWJ 8.9 vs 16.2. The Note's "about half in the small runs" holds |
| Missense (H→Y, R→W, P→S) 0.8 per 1,000 (L9) | same | mine: 4/4,119 = 1.0 (post 3/3,945 = 0.8). All first-position C→T missense: 5/7,554 = 0.7. CAG→TAG vs the missense trio: one-sided Fisher p = 5e-9 pooled; 0.018 (ALZF) and 0.0005 (JIWJ) in the two runs without LLXJ |
| About 96% identical; what looks like one strain (L9) | `rna_map/*.summary.txt`, `strain_compare.txt` | holds: read identity 95.55-95.84%. Runs agree at 97.9-99.0% of shared codons, and 88-90% of one run's differences from the genome recur in another (LLXJ→ALZF 1,483 of 1,693; LLXJ→JIWJ 2,257 of 2,522; ALZF→JIWJ 1,235 of 1,375) |
| Figure: all 32 reads carry TAG at enolase Q166 (L11, image) | `out/rna_codons.tsv`, row OZ237754.1:82154 (0-based) | holds: TAG 14 (LLXJ) + 1 (ALZF) + 17 (JIWJ) = 32, and no other codon |

**Is each alternative ruled out?**

- **Host dinoflagellate or another cobiont, including a ciliate. Ruled out.**
  - The host 18S is only on unplaced contig CBCXTV010000052.1. The *Bilabrum latius* 18S is on chromosome 12 (99.8%).
  - All TAG sites sit on OZ chromosomes, and 28 of the 30 contigs with a clean TAG also carry TAG-free HSPs. GC3 matches (0.771 with a clean TAG, 0.779 without).
  - TAA is not read (1 seen, 37.6 expected), so a UAR = Q ciliate is out.
  - A UAG-only ciliate (phyllopharyngeans) is out on the chromosome evidence. So is the ciliate MAG, which is a peritrich, a UAR = E code.
  - *Bilabrum*'s own eRF1 gene carries a TAG (opposite human L69, OZ237769.1:456067; `evidence.md` L18).
  - The 1KP reads that cover TAG sites come from the same organism as the rest: codon agreement with the genome is 87.7% within 30 codons of a clean TAG and 88.9% elsewhere (`near_tag_identity.py`).
- **Assembly or base-call errors. Ruled out.** Under a uniform point-error model, 37.6 TAA and 39.6 TGA would be expected; 1 and 4 are seen. An independent culture's mRNA carries TAG at 41 of 60 sites.
- **Pseudogenes. Ruled out**, by the same TAA/TGA deficit and and by the genes involved. The 91 clean TAG sites fall in 36 of the 76 conserved proteins, among them eIF4A (12 sites), RAN (7), enolase, actin, ARF1, EF-1α, VCP, PRPF8, proteasome subunits and RNA polymerase II. Around each site the 25-column window keeps a median of 13 residues identical to human (range 8-22).
- **RNA editing. Ruled out in the related strain.** It reads TAG at 646 of about 1,020 reads at clean sites, so its mRNA is not edited at those sites. Not tested in this strain; the post's L19 says "no RNA or protein from this strain".
- **Context-dependent or ambiguous code (TAG both sense and stop). Not excluded.**
  - The first in-frame stop within 20 codons of an aligned gene end is TAG at 1-2 *Bilabrum* loci (vs 19 in *Cafeteria*). The review's gene-end test found TAG as the first stop at 1 of 60 loci.
  - The post says this at L19 ("Whether TAG still ends a few genes").

## 3. What should be cut or held

- **Nothing needs cutting.** The headline, every bullet and the Known-before / Done / Not-known split all survive.
- **No hold.** The pending lanes (Logan SRA, the *Pirsonia* and *Mutimo* gap runs, *Maristentor* RNA) can't change "known before". At most they add a second case, which would be a new post.
- **Keep the courtesy note to Sanger ASG** before posting (Note L43). It is their unpublished genome, and no genome note exists yet.

## 4. Image

- **Correct.**
  - *Bilabrum* zoom TGC GCC TTT **TAG** GAG TTC ATG = OZ237754.1:82146-82166 (human ENOA 163-169, L A M **Q** E F M).
  - *Cafeteria* zoom CTC GCG TTT **CAG** GAA TTC ATG = VLTN01000067.1:94779-94799, the best *Cafeteria* enolase HSP. The second *Cafeteria* enolase copy, VLTN01000019.1, also has CAG at Q166.
  - The boxed position is human Q166 in both rows. The caption's "all 32 RNA reads" checks.
- **Reads at a glance:** yes. Title, two rows and one orange box carry the message, and the style matches the earlier TAG-thread barcode.
- **Two label issues (optional, but worth the minute):**
  1. **The unexplained second orange bar.** The *Bilabrum* barcode has a second orange bar near the left end of the gene. It is a TAG opposite human **K64** (OZ237754.1:81837), a clean site, and all **36** RNA reads there carry TAG (23 LLXJ + 13 JIWJ). Without a word, a reader may take it for a second stop. Suggested caption ending: "…carry TAG there. A second TAG, at a lysine (left), is also in all 36 reads."
  2. **"stops at TAG" beside a CAG box.** The *Cafeteria* subtitle sits next to a boxed CAG, and a quick reader can link the two. Suggested label: "uses TAG only as a stop", or "standard code: TAG = stop".
- **Better lead image?** No. The gene barcode is the right lead. An RNA-rate panel (CAG→TAG vs silent vs missense) is the strongest *argument*, but it needs reading, and the post states it in L9.

## 5. Wording fixes

**Required:** none.

**Optional:**

1. **L9: the cultures and the strain.**
   - Current: "Its reads, from three cultures that carry what looks like one strain (about 96% identical to this genome),"
   - Replacement: "Its reads, from three of the four cultures it contaminates (about 96% identical to this genome),"
   - Why: Vierock et al. found it in four cultures and called it "non-clonal". "What looks like one strain" is defensible, but it isn't needed.
2. **L9: the silent-rate comparison.** It rests mostly on the LLXJ run; the two other runs give about half the leucine rate, still about 10 times the missense rate.
   - Current: "about the rate of a silent leucine change (CTG to TTG, 14.4)"
   - Replacement: "close to the rate of a silent leucine change (CTG to TTG, 14.4)"
   - Or add after "run at 0.8 per 1,000": "The two smaller runs alone give the same picture, more weakly."
3. **L7: the gene-end table.** The numbers come from the earlier `scan/` ends table. The screen's own `scan2/` table gives 163 genes (TGA 147, TAA 16, TAG 0).
   - Either is fine; just keep the Note's file pointer consistent.
   - Or strengthen the line, since *Cafeteria* is the weakest of the four relatives.
   - Current: "against 11 of 145 in *Cafeteria*."
   - Replacement: "against 11 of 145 in *Cafeteria* and 35-39 in each of three other bicosoecids (149-172 genes each)." Source: `review/terminal_codons_relatives.tsv`.
4. **L20: the branch inference.** "Arose on the branch leading to *Bilabrum*" is a parsimony inference from four sampled relatives, stated inside "Not known".
   - Current: "so the change arose on the branch leading to *Bilabrum*; no other *Bilabrum* genome exists."
   - Replacement: "so the change arose after *Bilabrum* split from them; unsequenced close relatives may share it, and no other *Bilabrum* genome exists."
5. **L3: "single-celled".**
   - Current: "mostly ciliates and a handful of other single-celled lineages."
   - Replacement: "mostly ciliates and a handful of other lineages, nearly all single-celled."
   - Why: Cladophorales, one of the known code-change green algae, are multicellular filaments.
6. **L15 (optional addition).**
   - Current: "Some Labyrinthulea change their code in mitochondria only."
   - Replacement: "Some Labyrinthulea, and the Arctic *Plocamiomonas* (reported in 2026), change their code in mitochondria only."
   - Source: *Sci Rep* 2026, doi 10.1038/s41598-026-58736-z, *Plocamiomonas* mitochondrion on code 4.
7. **Note L42.** The Logan count is stale: README 12:05 says 58 of 129 at 12:02. Update it with the final numbers. The parent will send those.

## Could not check

- Lukeš 2025 *Curr Biol* Table S1 (403; nobody in the project has seen it).
- The SRA/Logan follow-up and the *Pirsonia*/*Mutimo* gap runs (still running), and *Maristentor* RNA. Pending by instruction.
- The exact motif rule behind "2,714 of 2,776" (`evidence.md` L20 says 2,776 of 3,189 entries carry the motif). I get 97-98% serine and 0 valine by two methods, but not the identical count.
- Protein-level evidence (no proteome exists), and RNA from this strain itself. The post already lists both as not known.

## Delta re-check of the Note (2026-11-12 14:08 EST; post body and figure unchanged)

Verdict: **FIX** (Note only; four wording fixes). Diffed against `scratchpad/bilabrum_post_prev2.md`; only L23 and L42-44 changed.

Holds, re-derived:
- *Maristentor* RNA (`maristentor/rna/ERR14711482.sites.tsv`, clean == 1, reads >= 3, same/reads): TAA at Glu 1,265 sites, 57,875/58,384 = 99.1%; TAG at Gln 365, 24,078/24,397 = 98.7%; TGA at Trp 184, 2,836/2,900 = 97.8%; GAA at Glu 99.4% (2,259); CAA at Gln 99.5% (854).
- eRF1 (`erf1_motifs.py` re-run): NODE_6418 (cov 137) KASNIKS, YLCDNKF, M51Q; NODE_12363 (cov 20) TASCIKS, YLCDSKF; both about 60% identical to human (`erf1_tx.tsv`).
- *Pirsonia* SRR34555172: 487,449 inner columns, TAG n = 0 (TAA 3, TGA 3). Ends (send == slen, not de-duplicated): TAA 353, TAG 223, TGA 124, other 296 = 996. *Mutimo* DRR963038: 81,253 inner, unpowered, TAG 0; ends TAA 64, TAG 61, TGA 70, other 52 = 247.
- Known codes: 8 units labelled known (*Remanella*, *Cristigera*, *Licnophora*, *Cephaleuros*, *Trentepohlia*, *Scotinosphaera*, *Boodlea*, *Nassula*); the other 6 have 276-94,695 inner columns.
- *Reophax* (`stage2/SRR12012189_TA[AG].txt`): GC3 0.53 vs 0.23 (p 1.2e-6), use distance 0.41 vs CGG control 0.07. *Dysteria* (`SRR9841576_TGA.txt`): dispersion 4.7; W sites in COX1, CYB, NU4M, COX2, ATP6, NU5M.
- 12 queued runs (`data/gap_nologan_order2.tsv`): none a stramenopile.
- Other stramenopiles in the Logan lane, all TAG n = 0: *Wobblia* 170k inner, *Incisomonas* 225k, *Develocauda* 118k, *Develocanicus* 71k, *Haloplacidia* 67k, *Feodosia* 25k.

Fixes:
1. "85 with enough aligned codons" -> "83". 85 of 132 units includes the calibration copies SRR33614924_L and SRR6032744_L. Of the 129 runs, 83 have >= 100k inner columns.
2. "the only two stramenopile orders in the gap list with no genome" -> "the only two stramenopiles among the 14 gap runs that Logan has no contigs for". Pirsoniales also has *Feodosia* in the Logan lane, and Logan covered other stramenopile orders.
3. "Two new signals failed" names 2 of 5. *Ploeotia*, *Eucyrtidium* and *Pterocladiella* (triage_full.md L55-57) also failed.
4. eRF1 "within the range of other ciliates": Q51 and C61 occur in none of the ciliate eRF1s in `erf1_queries.fa`; those positions vary among ciliates (51: L/M/R/Y; 61: N/S/H).

Optional: "the other 6 have under 100,000 aligned codons" implies power explains the misses. Yet 7 of the 8 hits were also under 100k, and *Furgasonia* has 0 TAA/TAG in 94,695. Separately, *Hemimastix* SRR6032744 is in the 12-run "no Logan contigs" queue, but SRR6032744_L was scored.

## Delta check 2 (14:25)

Verdict: **FIX** (two wording fixes; every number re-derives from raw files).

Holds: 129 Logan runs (132 units minus 3 `_L`), 83 with >= 100k inner columns (= column 31 "yes"). 14 queued = 13 in `gap_nologan.tsv` + *Hemimastix* calibration (has `SRR6032744_L`); 12 in order2, none a stramenopile, none started. *Pirsonia* 487,449 inner, TAG 0, ends 223/996; *Mutimo* TAG 0, 61/247. 8 known recoveries, 7 under 100k. `N_in` TAA/TAG/TGA: *Furgasonia* 39/31/44, 20/19/27 (TAA ends 103/218, 113/268); *Pseudomicrothorax* 19/21/25, 25/17/19; *Nassula* 1,767/1,056/8, no TAA/TAG ends. Five failed signals match triage. The two eRF1s are 55% identical to each other (two genes). Pan 2019 Methods quote confirmed; its Results: "three species representing two orders, i.e., Nassulida and Microthoracida". "Simplest reading" is fair parsimony.

Fixes:
1. "in 8 of the 14 gap units from lineages that have one, 7 of them from under 100,000 aligned codons; the other 6 did not call it. One of those, *Furgasonia*," -> "in 8 of the 12 gap units from lineages with a reported stop change, 7 of them from under 100,000 aligned codons (*Cristigera* only as a watch). Three of the other 4 did not call it (276 to 53,000 aligned codons); the fourth, *Furgasonia*,". Why: the 14 includes *Nyctotherus* and *Anteclevelandella* only via NCBI code 6 / phylum Ciliophora; the lab's ledger lists both as standard (Rotterova 2024), so "did not call it" reads as six misses.
2. "they differ from human at positions that also vary among ciliates, and nothing in them points to this code." -> "at the stop-recognition motifs they differ from human only at positions that also vary among ciliates, and no change is known to go with this code." Why: they are ~60% identical to human overall; Q51 and C61 (in TASNIKS) occur in no ciliate eRF1 compared, so "nothing points" overstates.

Optional: Pan 2019 also charted nassophorean stop-codon use (Fig. 4B, not seen): "...with the standard code and charted their stop-codon use". Give the reason: "*Furgasonia* and *Nassula* are both Nassulida, so the simplest reading...".
