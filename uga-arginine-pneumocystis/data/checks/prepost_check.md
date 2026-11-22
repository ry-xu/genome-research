# Pre-post check: mitochondrial UGA read as arginine in *Pneumocystis jirovecii* (new-stops-euk) (independent check, 2026-11-11)

Claim checked, from the outbox draft `outbox/2026-11-11-mito-uga-arginine/post.md` (file time 20:44; read again at 21:04):
- In the mitochondria of *Pneumocystis jirovecii*, UGA most likely means arginine, not tryptophan. *P. jirovecii* is filed under NCBI table 4, where UGA = tryptophan.
- "No NCBI genetic code reads UGA as arginine."

**What I read:**
- the post and `image.png`;
- README, the lead section ("Lead (2026-11-11): mitochondrial UGA read as arginine") and lines 718–756;
- `results/tga_arg/`: `rna/summary.txt`, `pneumo_cols.tsv`, `isolates.tsv` and the figure;
- `results/newleads/`: `lit.md`, `blind/report.md`, `blind/key.tsv` and the blinded agent's `final/` tables;
- the scripts behind them: `src/codes_census/sitescan_org.py`, `src/codes.py` (`blinded()`), `src/tga_arg/*.py` and `src/codes_census/reads_codons.py`.

**How I checked it:**
- I re-derived every number with my own scripts in `prepost/uga-arg_work/`, working from the census `.npz` files, the per-run RNA TSVs, the GenBank records and `others.fa`.
- A separate agent did most of the literature search. I re-checked the key points myself (section 1).

**What I touched:**
- I edited nothing else in `lab/` and messaged no agent.
- One slip: a bad glob wrote two Europe PMC downloads (`cisse2021.xml`, `ma2013.xml`) into `work/new-stops-euk/` for a few minutes. I moved them to my scratchpad, and the folder listing is back to what it was.

## Verdict: FIX, THEN POST

- **The core claim survives.** Every census number reproduces from the raw files with my own code (section 2).
  - *P. jirovecii* ranks first of all 22,698 genomes for UGA at conserved arginine positions.
  - Leaving *P. jirovecii* out of the consensus, or all 9 *Pneumocystis*, gives the same 23 positions. The family weighting is therefore not circular.
  - The RNA rules out editing. Every gene is full length and every UGA is internal, so UGA is not a misread stop or a frameshift.
- **Novelty holds** (section 1).
  - No report has UGA = arginine, RNA editing, or a genetic-code question in *Pneumocystis*. Ma et al. 2013 only noted the missing tRNA and read UGA as tryptophan.
  - No NCBI table reads UGA as arginine (gc.prt v4.6, 23 Sep 2024).
- **One bullet has wrong numbers, and four statements need fixing** (section 5):
  - **Isolates bullet:** one of the "19" isolates is a duplicate of the reference genome. All 18 real ones come from one NIH group, not "several labs". The "one site not found" is in fact found.
  - **RNA bullet:** CGA reads number 3, not 4.
  - **Editing comparison:** hornwort and fern organelles edit *this* codon, UGA to CGA, not "other codons".
  - **Disease line:** "often deadly in people with HIV" overstates.
  - **Figure:** the title and legend are clipped at the right edge.

## 1. Novelty

A separate agent searched by organism and by mechanism (web, Europe PMC, bioRxiv, PubMed). I re-checked the points the post depends on.

**By organism (*Pneumocystis*):** nothing earlier.
- **Ma et al. 2013** (*FASEB J* 27:1962–72, doi 10.1096/fj.12-224444) is the only paper to discuss UGA in *Pneumocystis*.
  - It sequenced the *P. jirovecii*, *P. carinii* and *P. murina* mitochondrial genomes (I confirmed this from the abstract).
  - Quoted via the agent and `lit.md`; I could get only the abstract myself: "there is no tRNA gene containing the corresponding anticodon UCA. It is possible that the C in the wobble position of the tRNA-Trp with anticodon CCA … is modified".
  - So the post's "Known before" line about Ma is accurate.
- **No other paper** reports UGA = arginine, mitochondrial RNA editing or any code question in *Pneumocystis*. That includes Cissé 2018 *mBio*, Cissé 2021 *Commun Biol* and the later isolate papers.
- **Cissé et al. 2021** (*Commun Biol* 4:305, doi 10.1038/s42003-021-01799-7). I read the text: "*P. macacae*, the closest living sister species to *P. jirovecii* identified to date". The post's use of it is correct. (MT726212, the prairie-dog genome, is not from this paper. It is a direct submission.)
- **PRIDE:** the keyword "Pneumocystis" returns only PXD030103 (mouse B-cell exosomes). The post's "no *Pneumocystis* proteome" holds. I checked PRIDE only, not MassIVE or jPOST.

**By mechanism (UGA = Arg anywhere; UGA read through editing):**
- **Yılmaz et al. 2026** (*Sci Rep*, doi 10.1038/s41598-026-62093-2), *Craticula* diatom mitochondrion.
  - The paper calls it a "putative reassignment", so "proposed" fits better than "reported" (optional fix).
  - My census agrees. *Craticula* PV031684.1 is the next-highest genome after *P. jirovecii*, with UGA at 15 conserved arginine positions.
- **Maturová, Žárský & Eliáš 2026** (bioRxiv 10.64898/2026.11.04.764322). I read the abstract from the bioRxiv API.
  - *Licnophora* and "certain other spirotricheans" read UGA as arginine in-frame and as a stop at gene ends.
  - Also: "the mechanism of UGA decoding as arginine in Licnophora … remains elusive". No tRNA was found there either, which is a useful parallel for the post's "Which tRNA reads UGA" (optional fix).
  - The same preprint reports UGA = leucine in *Balantidium*. The post cites it correctly.
- **Fang et al. 2021** (*Front Plant Sci* 12:780054), Trentepohliales plastids, and **Pichon, Luscombe & Plessy** (*F1000Research* 8:2072, v2 2020), *Bathochordaeus*, 3 codons: both accurately described.
- **Hornwort and fern organelles** edit UGA to CGA (U-to-C editing). This is the closest look-alike and the reason the RNA check matters.
  - In my census, hornwort and fern genomes are the next class after *Craticula*, at 8–14 UGA per genome at conserved arginine positions.
  - The post's "as plant organelles do with other codons" is wrong on this point (required fix 3).
- **Alfonzo et al. 1999** (*EMBO J*): in *Leishmania* mitochondria, the tRNA-Trp anticodon is edited C34→U so that it reads UGA as tryptophan.
  - This is a real precedent for Ma's modified-tRNA-Trp idea. A missing UCA tRNA is therefore weak evidence either way.
  - The post doesn't need to cite it, but the Risks note should know about it.
- **Baum & Beier 1998** (*NAR* 26:1390): a tRNA-Arg with a U\*CG anticodon reads UGA. *P. jirovecii* has a tRNA-Arg(UCG) gene (my scan finds it), so a candidate reader exists (optional fix).
- **NCBI gc.prt** v4.6 (23 Sep 2024): the standard code plus 26 alternatives. None has UGA = R. Both "26 alternatives" and "no NCBI genetic code reads UGA as arginine" hold.

**Bottom line:** "the first case in fungi, not the first anywhere" is right and fairly framed.

## 2. Correctness: the headline numbers

**How I re-derived them:**
- I reimplemented the consensus myself: family weight = 1/(number of genomes sharing the taxonomy minus its last rank).
- Each genome votes with its declared code, with sense-to-stop erased.
- Per-read RNA counts come from the 10 per-run TSVs. Isolates were done two ways: the project's flank method, and my own whole-genome minimap2 (asm20, `--cs`).

| post claim | source | my result | check |
|---|---|---|---|
| About 22,700 mitochondrial genomes | the 3 census `.npz` files | 22,698 genomes, 84,424,517 placed sites, no accession counted twice | match |
| NCBI lists 26 alternatives; none has UGA = Arg | NCBI gc.prt v4.6 | 27 tables (standard + 26); no UGA→R | match |
| UGA at 23 positions with ≥ 90% arginine | census | 23 | match |
| At none with tryptophan; this genome spells Trp UGG at 37 | census | UGA at conserved Trp positions 0; UGG at conserved Trp positions 37 | match |
| 36 UGA in aligned regions, 32 where arginine is most common | census | 36 placed; the most common amino acid there is R 32, K 2, Y 1, S 1. Every UGA at a strictly conserved position is at an arginine one. | match |
| Relatives at the 23: Arg codon 179 of 184, UGA 4 (3 *P. carinii*, 1 *P. murina*), lysine codon 1, UGG never | census, 8 relatives | AGA 113, CGA 48, CGT 9, AGG 5, CGG 4 (= 179); UGA 4 (*carinii* at cox3:164, nad4:373 and cox1:191; *murina* at nad2:310); AAA 1 (PP484879); UGG 0 | match |
| Relatives' hosts (mouse, two "rat species", rabbit, dog, macaque, prairie dog, squirrel) | GenBank records | *murina* mouse; *carinii* and *wakefieldiae* both rat; *oryctolagi* rabbit; *canis* dog (NC_060320, Finland); MT726214 rhesus macaque; MT726212 black-tailed prairie dog (Kansas); PP484879 Finlayson's squirrel (Italy) | match, but "two rat species" reads as two species of rat (optional fix) |
| **19 more isolates, "from several labs", UGA 595 of 608; 12 AGA at one site; one site not found** | `others.fa` | See the isolate notes below the table. | **mismatch** |
| RNA: 7,109 of 7,298 reads (97.4%) read UGA at the 32 sites | 10 per-run TSVs (PRJEB3063) | 7,109 / 7,298 = 97.41% | match |
| **"and 4 read CGA"** | same | **3** at the 32 sites. The 4th CGA read is at nad5:537, a site at a tyrosine-consensus column and not one of the 32. One more read is lost to the top-3 truncation in `reads_codons.py`. | **mismatch** |
| The patient's strain has AGA at that site, 129 of 130 | same | nad2:282: AGA 129 of 130, with AGA in every one of the 10 runs | match |
| Without that site, 99.2% UGA | same | 7,109 / 7,168 = 99.18%. At the strictly conserved sites only: 5,028 / 5,063 = 99.31%. Lowest depth at any site: 15. | match |
| No tRNA with anticodon UCA | my anticodon-loop scan; the blinded agent's ARAGORN | Trp(CCA) present, also Arg(UCU) and Arg(UCG); no UCA | match |
| Blinded re-check: Arg at 26 of 29 conserved positions, p 2e-13; both controls table 4 | agent `final/tests.tsv`, `columns.tsv`, `report.md`, `key.tsv` | G04 = *P. jirovecii*. TGA at 29 conserved columns: R 26, K 1, W 1, Y 1. p = 2.15e-13. Codetta also calls R. Pass rule written down (14:36) before launch: PASS. | match (but see the p-value note below) |
| *P. carinii* RNA: 50 of 50 UGA at 4 sites | SRR2156857 tally | 17 + 22 + 8 + 3 = 50, all UGA | match |
| Macaque species is the closest relative (Cissé 2021) | the paper | quote verified (section 1) | match |
| Figure bars and footnote | census | per genome, UGA at its own conserved-Arg positions / regular Arg codons there / UGG at conserved-Trp positions: *jirovecii* 23/21/37, *carinii* 4/41/35, *murina* 1/45/35, *wakefieldiae* 0/47/37, *oryctolagi* 0/47/37, *canis* 0/45/36, macaque 0/46/37, prairie dog 0/48/37, squirrel 0/47/37. No genome has UGA at a conserved Trp position. | match |

**Isolate notes:**
- The project's flank method reproduces 595 / 12 / 1. My whole-genome alignment finds the missing site, JX855938 atp8:41, and it is UGA, so all 608 sites are found: 596 UGA + 12 AGA.
- **JX499143 is the reference genome itself.** The RefSeq comment on NC_020331 reads "identical to JX499143", and JX499143.1 is in `others.fa`. That leaves 18 independent isolates: UGA **564 of 576** (97.9%) and AGA 12, all at nad2:282 (6 UGA and 12 AGA there).
- **All 18 were sequenced or submitted by the Kovacs/Ma group at NIH:**
  - JX855936–8: Ma 2013;
  - MH010437–46: Cissé 2018 *mBio*, patients in the US, Germany and Switzerland;
  - MT319762: China, submitted from NIH;
  - MW530525–8: Azar et al. 2021, Connecticut.

**Not re-checked:** "all of RefSeq plus one per additional species from GenBank" (I checked only the genome count). Also the README's earlier "697 times", which the Note already flags as not reproducing.

**The p-value.** The blinded p tests "UGA is not used like its table-4 synonym" (i.e. not tryptophan). It does not test "is arginine". Codetta's call is what says arginine.

**The blinded set-up has three caveats:**
- The agent identified each coded genome to species by blastn against RefSeq (99.95–100%), so it knew G04 was *P. jirovecii*. Its ortholog references excluded *Pneumocystis*.
- Its calling rule was settled after it saw the data. G04 passes under both versions: relaxed 31/34, strict 22/23.
- The post's "standard fungal code" for the controls means table 4 (the fungal mitochondrial code), not the standard code (table 1).

**Is UGA a stop or a frameshift artefact?** No.
- *P. jirovecii*'s aligned span for every gene matches its relatives' spans. For example cox1 covers model columns 23–529 (relatives 22–518 to 22–530) and cob 38–417 (relatives 38–416 to 38–419).
- Every one of its 36 aligned UGA codons is internal and has aligned, in-frame columns downstream: cob:106 has 311 downstream and nad2:33 has 407. The shortest is atp8:42 in a 43-column model.
- If UGA were a stop, the 10 of its 15 protein genes that contain UGA would all be cut short. They are not.
- In the blinded agent's gene-end table for G04, the genes whose end is followed directly by a stop use TAA or TAG. None ends at UGA.
- The RNA shows the same UGA, and so do the 18 isolates.

**Is the family weighting circular?** No.
- *P. jirovecii* votes **tryptophan** at its own UGA codons: it declares table 4, and the census votes with the declared code. So it works against the claim, not for it.
- The 9 *Pneumocystis* genomes share one family weight. At the 23 positions they carry 0.08–0.24% of the total weight, out of 4,162–6,470 genomes in 424–1,228 families.
- At those positions, the arginine share runs from 0.938 to 0.999 and the tryptophan share is at most 0.0017.
- Leaving out *P. jirovecii*, or all 9 *Pneumocystis*, gives exactly the same 23 positions and the same counts.

**How unusual is it?**
- *P. jirovecii* (23) tops all 22,698 genomes. Next come *Craticula* (15, published) and hornwort and fern organelles (8–14, known editing).
- Among 4,630 fungi: 4,557 have 0, 71 have 1, *P. carinii* has 4 and *P. jirovecii* 23. No fungus outside *Pneumocystis* has more than 1.
- Using UGG only for tryptophan is common in itself: 1,059 of 4,486 fungi have no UGA at conserved Trp positions.
- So "no UGA at Trp" alone means little. The 23 positions do the work.

## 3. Weak parts

1. **No protein evidence.** The post says so. The one live alternative is Ma's: a modified tRNA-Trp(CCA) that reads UGA as tryptophan.
   - That is not far-fetched. *Leishmania* edits its tRNA-Trp anticodon to read UGA (Alfonzo 1999). A missing UCA tRNA is therefore weak evidence for either reading.
   - **What argues against tryptophan is the positions.** The 23 positions are 94–99.9% arginine across about 4,200–6,500 genomes.
   - One of them is cox1 model column 457: 99.5% arginine across 6,188 genomes. My pairwise alignment to cow COX1 puts it at **R439 of the conserved MPRR pair**, where R438/R439 bind the heme propionates. The codons there read `ATG CCC AGA TGA`, i.e. M P R [UGA].
   - "Most likely" is the right strength.
2. **The isolate evidence is less independent than it reads.**
   - All 18 assemblies come from one group, possibly with a shared pipeline, so they cannot fully rule out that the reference's codons were carried into them.
   - Against that, the isolates differ from the reference at nad2:282 and nad1:33, so they are not copies.
   - The RNA (PRJEB3063) is independent: a University of Lausanne lung-wash sample, total RNA, NuGEN-amplified, HiSeq 2000.
3. **Within-species swaps support UGA = Arg, and the post doesn't use them** (optional fix). UGA and a regular arginine codon swap within *P. jirovecii* in both directions, as synonymous codons do:
   - nad2:282: UGA in the reference, AGA in 12 of 18 isolates and in the RNA patient.
   - nad1:33: CGA in the reference (strictly conserved Arg column nad1:40), but UGA in 41 of 41 RNA reads and in 6 of 18 isolates.
4. **RNA mapping is loose but adequate.**
   - Mapping quality 0 is accepted, against a gene-only reference with ≤ 5% divergence.
   - My controls read at background level: UGG sites read UGG 98.0%, CGU sites CGU 99.1%. The 97.4% at the UGA sites is at that level, and 99.2% once the polymorphic site is dropped.
   - Reads are counted one by one, so the two overlapping halves of a read pair can be counted twice. That inflates n, not the percentages.
5. ***P. carinii* is thin.**
   - 4 sites, one with only 3 RNA reads.
   - Its 4 is more than any fungus outside *Pneumocystis* (max 1), so it may mean something.
   - 3 of its 4 fall on *P. jirovecii*'s 23 positions, about what chance gives (expected ≈ 2). That is no evidence of a shared origin.
   - "Too few to call a code" is right.
6. **"Changed in the *P. jirovecii* lineage"** rests on the macaque species being the closest relative (verified) and having 0. This is fine as worded ("most of the 23 sites most likely").
7. **"Never where they have tryptophan" depends on the column definition.** It is true for the census's strict columns, but with its own column definition the blinded agent found 1 UGA at a tryptophan-majority column (against UGG 41). "At none of the 37" is the safer wording.

## 4. Image

`outbox/.../image.png` is identical (md5) to `results/tga_arg/fig/pneumo_uga_arg.png`. It shows horizontal bars per *Pneumocystis* genome: UGA (blue) plus regular arginine codons (grey) at that genome's conserved-arginine positions.
- **Correct:** every bar and label matches my per-genome counts (section 2), and the footnote is correct.
- **Broken:** at 1720 px, the title is cut off at the right ("…in the human specie"), and so is the second legend entry ("regular arginine codons (CGN, A").
  - In `src/tga_arg/fig_pneumo.py` (title at line 34, `savefig` at line 40), shorten the title or save with `bbox_inches='tight'`.
  - A shorter title: "UGA replaces arginine codons in human *Pneumocystis*".
- **Does it read at a glance?** Mostly. The blue 23 against a column of zeros is clear. But the other half of the claim (never at tryptophan) lives only in the footnote.
  - Optional: a narrow second column per row, "UGA at conserved Trp: 0".
- **Possible lead image** (in the style of the gene-barcode figure Ryan used before):
  - A codon strip at cox1's MPRR arginine pair (cow R438/R439, the heme-binding arginines; 99.5% Arg across 6,188 mitochondrial genomes).
  - Rows: cow, a yeast, *P. murina* or the macaque species with a regular arginine codon, then *P. jirovecii* `ATG CCC AGA UGA`.
  - The *P. jirovecii* UGA cell in the accent colour, sitting in the same column as the others' arginine codons.
  - I did not pull the relatives' codons at that column or build the figure.

## 5. Wording fixes

Required:
1. **Isolates bullet (line 8).**
   - Current: "19 more *P. jirovecii* mitochondrial genomes, from several labs, carry UGA at the 32 arginine-position sites 595 of 608 times. The exceptions are 12 AGA (an arginine codon), all at one site, and one site not found."
   - Replace with: "18 more *P. jirovecii* mitochondrial genomes, from patients in the US, Europe and China (all sequenced by one NIH group), carry UGA at the 32 arginine-position sites 564 of 576 times; the 12 exceptions are AGA (an arginine codon), all at one site."
   - Also change "19 more" to "18 more" in "Done → Isolates and relatives" (line 20), and change "19 isolates" and "595 of 608" in the Note.
   - Reason: JX499143 is the reference genome, and whole-genome alignment finds the "not found" site (it is UGA).
2. **RNA bullet (line 9):** change "and 4 read CGA" to "and 3 read CGA".
3. **Editing comparison (line 9).**
   - Current: "as plant organelles do with other codons."
   - Replace with: "as hornwort and fern organelles do with this very codon (UGA to CGA)."
   - Reason: hornwort and fern organelles edit UGA itself, which is exactly why the RNA check matters.
4. **Disease (line 5).**
   - Current: "a pneumonia that is common and often deadly in people with HIV or on immune-suppressing drugs".
   - Replace with: "a pneumonia that can be deadly in people with weakened immune systems, from HIV or immune-suppressing drugs".
   - Reason: with treatment, hospital mortality is about 4% in AIDS against about 27% in other patients (Roux et al. 2014, *Emerg Infect Dis*). Recent reviews give 5–30% in HIV and 4–76% without HIV. It is deadlier *without* HIV.
5. **Blinded bullet (line 11):** change "the standard fungal code" to "the usual fungal mitochondrial code (table 4)". "Standard" reads as the standard code, table 1.
6. **tRNA bullet (line 10).** This is Ma et al.'s 2013 observation, not a new finding, and a missing UCA tRNA does not by itself argue for arginine (see the *Leishmania* precedent).
   - Move it into "Known before" (Ma 2013 is already there), or reword it to: "As Ma et al. (2013) noted, the genome has no tRNA with anticodon UCA, the usual reader of UGA as tryptophan in mitochondria."
7. **Figure:** fix the clipped title and legend (section 4).

Optional:
- **Add the within-species swaps** as a bullet. For example: "Inside *P. jirovecii*, UGA and regular arginine codons swap places at two sites: one where some isolates have AGA, and one where the reference genome has CGA but the RNA patient's strain reads UGA (41 of 41 reads)."
- **Add scale** to the per-position bullet: "(94–99.9% arginine across about 4,200–6,500 genomes)".
- **Per-position bullet:** change "never where they have tryptophan" (line 5) to "at none of the 37 positions where they have tryptophan" (weak part 7).
- **Relatives bullet:** change "from mouse, two rat species, rabbit, dog, …" to "two from rats, and one each from mouse, rabbit, dog, macaque, prairie dog and squirrel". Both rat *Pneumocystis* live in the same host species.
- **Blinded bullet:**
  - "p 2e-13" → "p 2e-13 against tryptophan".
  - Optionally add "(it identified the species itself by BLAST)".
- **"Known before":** change "Reported in the mitochondrion of the diatom *Craticula vanensis*" to "Proposed in…". The paper calls it a "putative reassignment".
- **"Not known → Which tRNA reads UGA":** add "The genome does have a tRNA-Arg with anticodon UCG, which can read UGA in some settings (Baum & Beier 1998); *Licnophora*'s UGA-arginine reader is also unknown (Maturová et al. 2026)."
- **Note → Risks:** add that a tRNA-Trp with an edited anticodon reads UGA in *Leishmania* mitochondria (Alfonzo 1999). Ma's idea has a precedent, and only peptides can exclude it.
- **Define terms** on first use for a general reader: tRNA and anticodon ("the three letters on a tRNA that pair with the codon"), and ENA and PRIDE (public sequence and protein archives).
- ***P. carinii*:** worth one clause: its 4 is more than any other fungus (all others have 0 or 1).

Not errors:
- The "Known before" list is accurate (Ma 2013, Maturová, Fang, Pichon). Only the Yılmaz verb is optional.
- The Note's description of the blinded test (PASS, controls called table 4) matches the agent's files.
