# Fact-check: Bilabrum TAG = Gln post

Fresh subagent a6d654dfeaa62a243, finished 2026-11-11 ~23:43. It could not write files, so Myrtle saved its text report here at 23:45 without changes to content. Line numbers (L) refer to scratchpad/bil_post_draft.md as it stood at 23:3x (pooled 1KP RNA numbers).

## Verdict

The claim appears new and is right in substance, but the draft needs fixes before posting. No earlier report of a nuclear code change in *Bilabrum*, any bicosoecid, Bigyra/Opalozoa or any stramenopile was found (2014-2026, preprints included). The only earlier sighting is Vierock et al. 2022, who saw this TAG in a *Bilabrum* channelrhodopsin and called it a premature stop. The draft credits that paper to the wrong authors.

Core numbers match the project files: 166 in-frame TAG vs 1 TAA and 4 TGA; Gln at 44 of 91 clean sites; 0 of 149 gene ends TAG; eRF1 S60V; pooled RNA 41/60 (recomputed, no ties, all strict majorities) and 12.3 / 14.4 / 0.8 per 1,000. The 0.8 comparison is conservative: all protein-changing C->T changes at codon position 1 give 0.4 per 1,000.

Must fix:
1. The *Sci Adv* 2022 paper is Vierock et al., not Govorunova et al. (post L16, "Would care" L30, evidence.md L61, run_1kp.sh comment).
2. "Tested on 67 genomes of known code": the 67 are genomes and transcriptomes (43 positives, 24 negatives), and the rules were set on them (README L55-64). *Spironucleus* was relabelled after its score was seen (README L61).
3. "Four other bicosoecid genera": three named genera (*Cafeteria*, *Caecitellus*, *Halocafeteria*) plus one uncultured Cafeteriaceae genome (GCA_964026335.1).
4. The tRNA sentence is overstated (item 14).

## Per-sentence check

1. L1 claim - OK (evidence.md L5-8, L42-48; review/bil.analysis.txt L4-8). Marine: Schoenle et al. 2020; Baričević et al. 2023.
2. L1 "No change ... in stramenopiles" - OK (notes/known-nuclear-stop-codes.md §1.2; Lukeš 2025 doi 10.1016/j.cub.2025.09.071; novelty search). Keep "nuclear": stramenopile mitochondrial changes exist.
3. L3 "mostly ciliates and a handful of other lineages" - OK (ledger rows 1-34).
4. L3 "fixed before scoring and tested on 67 genomes" - Wrong/overstated. Use "call rules set on 67 genomes and transcriptomes of known code and frozen before the screen ran".
5. L3 "bacteria-eating" - Unsupported. Schoenle 2020 (doi 10.1016/j.ejop.2019.125665): heterotrophic *Cafeteria*-like nanoflagellates. Use "heterotrophic".
6. L3 "released March 2025 (GCA_965178125.1, 43 chromosomes)" - OK: 2025-03-04, Wellcome Sanger, 43 chromosomes, 48.53 Mb, PRJEB86245 (factcheck/dl/web1/asm_bil.json). GCA_965178095.1 is the alternate haplotype.
7. L3 "clearest new case" - OK as judgement (README L171-178, Lead 1); fully checkable only once [N] and the screen are final.
8. L5 166 / 1 / 4 - OK (bil.analysis.txt L4; evidence.md L44-46).
9. L5 no TAG in *Cafeteria* - OK (TAG 0, TGA 13). NCBI: *Bilabrum* under "Bicosoecida incertae sedis".
10. L6 Gln 44 of 91, Lys 11 - OK (bil.analysis.txt L5).
11. L6 "pattern of CAG" - OK: best-fit codon profile Gln 1.005/site vs Lys 0.120 (bil.analysis.txt L8).
12. L7 "It no longer ends genes" - Overstated; 0 of 149 does not exclude rare TAG stops and conflicts with L19. Use "It ends none of 149 genes checked".
13. L7 0 of 149 (TGA 136, TAA 13) vs 11 of 145 *Cafeteria* - OK (evidence.md L8-9; terminal_codons.py). *Cafeteria* has the lowest TAG share of the relatives (others 20-24%), so this is the conservative comparison.
14. L8 tRNA "scores as glutamine-like as the genome's own glutamine tRNA" - Overstated. Gln 73.1 vs tRNA-Gln(TTG) 70.4 (review/iso.tbl), but it scores His 73.2 equally; the genome's own tRNA-Gln(CTG) scores His 66.5 > Gln 56.1, so the models do not separate Gln/His here. By sequence it is closest to tRNA-Tyr (68.5% vs 58.9% to tRNA-Gln; evidence.md L36-37). A second CUA copy is a pseudogene (trna.out L11). Use "a tRNA with anticodon CUA, which pairs with UAG; identity models score it as glutamine- or histidine-like, not tyrosine-like".
15. L8 eRF1 valine at S60 kept by 2,714 of 2,776 - OK (evidence.md L19-20; human numbering, TASNIKS 58-64). Recounts with other motif patterns: 2,725/2,766 and 2,761/2,821; Val 0 every time. Add "in UniProt reference proteomes". The eRF1 gene itself has an in-frame TAG at L69.
16. L9 no RNA-seq for this genome - OK: ENA has only 7 host runs (PacBio, Hi-C) for *Bilabrum* and *Dinothrix rugata* (factcheck/dl/web2/ena_runs.tsv).
17. L9 related *Bilabrum* in 1KP algal cultures - OK: 18S 99-100% to *B. latius* HFCC35; this genome 99.6% (evidence.md L12-13, L61-63). Hosts: one red, three green algae.
18. L9 "three cultures that carry the same strain" - Slightly overstated. Runs agree at 97.9-99.0% of codons vs ~89.5% with the genome (strain_compare.txt); no pair disagrees TAG vs CAG at sites covered in two runs (runagree.py). But B1ChR2 has two haplotypes differing at 32 of 1,171 nt, including TAG/CAG at codon 99. Use "what looks like one strain".
19. L9 "about 96% identical" - OK: mean read identity 95.6-95.8% (ERR*.summary.txt); read-level, includes sequencing error, mapped regions only.
20. L9 TAG 41 of 60, CAG 18 - OK (all3.table.txt; 60th site CAA; covered = >=2 pooled reads). No ties; 3 sites have both TAG and CAG reads; at >=5 reads 33 TAG, 12 CAG, 1 CAA. 58 of 60 sites come from LLXJ; ALZF alone gives 5 TAG, 5 CAG.
21. L9 "CAG becomes TAG at 12.3 ... CTG to TTG 14.4" - Numbers OK (23/1,875 vs 36/2,506, ALL3.pairs.txt; two-sided Fisher p = 0.60); "becomes" implies a direction. Use "where this genome has CAG, the 1KP strain has TAG at 12.3 per 1,000 sites". 21 of 23 events from LLXJ; the small runs show TAG at about half the Leu rate (5.7 vs 10.9; 8.1 vs 15.8). Reverse: TAG->CAG 270 per 1,000, TTG->CTG 215.
22. L9 His->Tyr, Arg->Trp, Pro->Ser 0.8 per 1,000 - OK (3/3,945; p = 3.8e-9). All C->T at codon position 1, like CAG->TAG. Hand-picked but conservative: all such missense give 3/7,320 (0.4 per 1,000).
23. L9 "A stop codon ... would be weeded out like those; TAG is not" - OK as inference. Direct test underpowered: CAA->TAA + CGA->TGA 1/375 (p = 0.16). Sense -> TAG 25 times, -> TAA 0, -> TGA 1 (ALL3.pairs.txt).
24. L11 image (MCM2 Q607, OZ237759.1:64848) - clean site, but the only 1KP read there reads CAG. Better sites with strong TAG support and no CAG reads: ENOA Q166 (OZ237754.1:82155, 32 TAG reads); PSA6 Q123 (OZ237746.1:74495, 35); PRS6B Q340 (OZ237753.1:3443, 22). Check *Cafeteria* has CAG at the chosen site.
25. L14 known lineages - OK (ledger; Lukeš 2025).
26. L14 TAG(+TAA)=Gln arisen many times - OK (Pánek 2017 BMC Biol doi 10.1186/s12915-017-0353-y).
27. L15 Swart 2016 MMETSP - OK (Cell 166:691, PMC4967479: "no new codes in the remaining 265 eukaryotes"); *Cafeteria* not named but MMETSP includes *Cafeteria* and Bicosoecida sp. CB-2014.
28. L15 Labyrinthulea mito only - OK; use "Some Labyrinthulea" (Žihala 2020 doi 10.1016/j.ympev.2020.106908).
29. L15 *Blastocystis* poly(A) - OK (Klimeš 2014 doi 10.1093/gbe/evu146, ~15% of nuclear genes).
30. L15 NCBI standard code for every *Bilabrum* taxon - OK (5 taxa, code 1, NCBI default).
31. L16 "Govorunova et al. 2022 (*Sci Adv*)" - Wrong. Vierock J, ..., Rozenberg A, ..., Béjà O, Hegemann P. Sci Adv 8:eadd7729 (2022), doi 10.1126/sciadv.add7729, PMC9733931.
32. L16 "One transcript ... second copy" - Overstated. The quote applies to the original transcript "and some of other B1ChR2 variants"; the protein came from a phased stop-free haplotype assembled from reads.
33. L16 "TAG where the second copy has CAG and the solved structures have glutamine" - Facts OK (8JXO/8JXP model residues 2-263), but the structures used the stop-free construct, so not independent evidence. Say "haplotype", drop the structure clause.
34. L17 Done - OK (rules set on the controls).
35. L19 - OK; conflicts with L7 heading.
36. L20 "four other bicosoecid genera" - Wrong count: *Cafeteria*, *Caecitellus* (GCA_965644315.1), *Halocafeteria* (GCA_057816405.1), uncultured Cafeteriaceae (GCA_964026335.1).
37. L26 Sanger ASG cobiont, unpublished - OK (PRJEB86245 under PRJEB84068 "*Galeidinium rugatum*"); no genome note found.
38. L27 "95.8% identity" - Stale; use "about 96% (95.6-95.8%)".
39. L30 Would care - Wrong/unsupported. Channelrhodopsin: Vierock and Hegemann (Berlin), Rozenberg and Béjà (Technion); structures Zhang and Pei. Arndt group (Cologne) described *Bilabrum* (Schoenle 2020); "holds cultures" unverified. Optional: Pfannkuchen group (Baričević 2023, doi 10.1016/j.ejop.2023.125999); Konstanz source of the *Dinothrix* culture.

## Novelty search

No prior report found (Europe PMC incl. 2024-2026 preprints, web; saved in factcheck/dl/web2/).

- Earlier sightings, not called a code change: Vierock 2022 (doi 10.1126/sciadv.add7729), treated as premature stop; Zhang...Pei bioRxiv 2023 (doi 10.1101/2023.07.22.550149) and PDB 8JXO/8JXP, stop-free construct.
- Later channelrhodopsin papers mention B1ChR2, nothing on codons: Govorunova 2022 mBio (PMC9765531); Govorunova 2024 (PMC10865875); Shan 2024 Nat Commun (PMC11344840); Tajima 2023 Cell; Mermet-Joret 2025 eLife.
- Code surveys/reviews, no stramenopile nuclear change: Swart 2016; Pánek 2017; McGowan 2023/2024; Záhonová 2016/2025; Lukeš 2025 (Table S1 not seen); Maturová 2026 bioRxiv (doi 10.64898/2026.11.04.764322, ciliates only); Codetta (Shulgina & Eddy 2021, bacteria/archaea); EukPhylo 2025; Bellas 2023 PNAS.
- Not nuclear code changes: Žihala 2020; Klimeš 2014.
- Bicosoecid/*Bilabrum* papers, standard code or no code analysis: Schoenle 2020; Baričević 2023; *Halocafeteria* genome preprint (doi 10.64898/2026.06.25.734631); bicosoecid virophage preprint (doi 10.64898/2026.07.28.741278); Yamada 2025 PNAS.
- Remaining risk: Sanger could publish the cobiont with its own code analysis; nothing yet (no genome note, no RNA-seq at ENA).

Check scripts: <agent scratchpad>/fc_scripts/ (ties3.py, ct_rates.py, runagree.py, topq.py). Downloads: factcheck/dl/web1/, factcheck/dl/web2/.
