# Audit of the AGA = Ala lead (fresh subagent, finished 2026-11-11; saved verbatim by Clover)

## Verdict
No code bug creates or inflates the signal. Frame, strand, codon indexing, the strict rule, the Fisher contrast, the class
definitions and the tRNA parse (positions 3:70 and N73) are all correct. The consensus reads these genomes' own AGA as Arg
(table 4), which works against the claim. Hand-drawn stems check out: T. inkin has G3·U70, A73 and anticodon UCU, and it
shares its acceptor with a separate tRNA-Ala. A reproduction that doesn't use the census or HMMs agrees in every genome
tested. Almost all numbers reproduce exactly. The problems are wording and inference.

## Wrong or loose
1. "Pooled AGA-at-Arg counts come from the Arg-reading genomes" is wrong. In Cystobasidiomycetes all 3 come from MT366950.1,
   and in Exobasidiomycetes 1 of 27 comes from BK074393.1. MT366950.1 is the same genome as BK072865.1, which is in the 32
   (rotated, 99.6% identical). Its 3 AGA-at-Arg codons are CGA in BK072865, so they are assembly differences. Corrected,
   Cystobasidiomycetes is about 9/0/0/247, and no genome in that class reads AGA as Arg.
2. The p-values (1e-154 etc.) are inflated by pseudoreplication. The 32 genomes are about 29 species, and several Tilletia
   are near-identical. The 124 AGA-at-Ala sites are only 24 distinct columns. Counting each column once per lineage gives
   p = 5e-16 (Trichosporonales), 3e-12, 3e-8, 3e-7 and 0.028 (Sporidiobolales).
3. Pitch "sits where related proteins have alanine": 78% of the time (124 of 159: G 15, S 7, P 5, L 4, V 2, T 1), while
   GCN is 97.5% Ala.
4. Pitch "never where they have arginine" is true by construction (class requires 0 AGA at Arg columns).
5. "0 of 4,387": two carriers have AGA at Arg columns (MT366950.1, BK074393.1), so the class definition must be attached.
   22 of the 4,387 have no ARAGORN output and 65 have no UCU tRNA. Fair count: 0/4,322 with a UCU tRNA.
6. "(67, 39, 9)" mixes units. Class-wide counts are 80, 39 and 12; 67 is Trichosporonales alone.
7. "AGG is 0 at strict columns in all three" is off by one: Exobasidiomycetes has 1, at a Met column.
8. "5 same-class genomes" are 4 Exobasidiomycetes plus Solicoccozyma. Within Tremellomycetes almost all relatives have no
   AGA at all, so the comparison there is mostly against AGA-free genomes (codon lost, then reused as Ala?).
9. "Probably three origins" is untested. 9 AGA-at-Ala columns are shared across the three groups, vs ~4 expected by chance
   (P ≈ 7e-4). Those columns are no less conserved than other Ala columns. Shared ancestry or codon hotspots both open.
10. Auricularia NC_058230.1 has a tRNA(UCU) identical over 74 nt to T. fuciformis's: likely contamination; don't list it
    as a carrier.

## Independent reproduction (tblastn of NCBI reference proteins, code 4; AGA at unanimous Ala / Arg columns)
| genome | AGA at Ala / Arg (of unanimous) | CGN at Arg |
|---|---|---|
| T. inkin | 18 / 0 (27) | 71/77 |
| A. gamsii | 26 / 0 (38) | 74/75 |
| T. mesenterica | 29 / 1* (35) | 76/77 |
| T. indica | 14 / 0 (21) | 65/65 |
| C. benthicum | 2 / 0 (4) | 79/79 |
| C. minutum | 3 / 0 (4) | 78/78 |
| Controls (Solicoccozyma, Exobasidium, Symmetrospora) | 0 / 17 (18) | — |
*method artefact (hit ran into intron; real exon has CGA). Pooled AGA vs CGN p = 4e-104. AGA at Ala 71% vs GCN 93%.

## Decoy
Every other codon is >= 83% its standard residue in the census; AGA is the only one with a non-standard majority. In the
reproduction the only other hit is ATA at internal Met columns in A. gamsii (14/26) and T. mesenterica (9/19), already
flagged by the blinded re-check. ATG->ATA is one change; GCN->AGA needs two.
Side finding: Sympodiomycopsis PP995850.1 puts AGA at Ser 20 of 25 times and never at Arg (possible AGA = Ser).
Not rechecked: Salmacisia and Pascua explanations, literature, protein run. Process slip: one 10-line read-only print with
system python3.

## Appendix
- Code checked: census_org.py, sitescan_org.py, codon_contrast.py, trna_vs_usage.py, org_trnas.py acceptor(), codes.py.
  Acceptor parse calls G:T on 4,585 of 4,657 tRNA-Ala. T. mesenterica's tRNA reading (offset 1, trim 1) supported by U8
  and T-loop TTC. All 36 G3:U70 tRNA(UCU) have A73. Hannaella's UCU tRNA only at threshold 0.8 (minus 18525–18594, G:T,
  A73), overlapping a weak plus-strand call.
- Exact matches: per-genome counts for all 4,638 fungal genomes, lineage table, all three contrast p-values, all-columns
  shares (AGA R 13 of 1,321), raw reads 54/54 (depth 2–9), CP116791.1 60/7 (6 A, 1 L), Ceraceosorus ATA = Met.
- Reproduction references: Cryptococcus NC_004336.1, Jaminaea NC_023248.1, Rhodotorula NC_036340.1. 0 mismatches between
  extracted codons and tblastn translation. An earlier tblastx attempt was invalid (aligned rRNA, off-frame).
- Shared columns: cox2:32, nad1:200/212, nad3:111, cox1:112/222/331/343/393.
- Files: star.summary.tsv, star.decoy.tsv, star.stats.txt, decoy.tsv, cache.npz (556 MB, rebuilt by cache.py).
