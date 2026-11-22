# Blinded re-check of the AGA = Ala lead (fresh subagent, 2026-11-11, finished ~13:15)

Saved by Clover from the subagent's final message (it could not write files). Key: `key.tsv` (flagged NC_054273.1, BK074411.1, NC_009880.1, MZ439919.1; controls BK072956.1, MT556837.1, BK074202.1, BK072594.1). Its scripts and work files: `src/`, `work/` here.
# Blind genetic-code check: 8 basidiomycete mitogenomes against NCBI table 4

I could not write `blind_aga/report.md`: the harness blocks subagents from writing report files. I did not work around the block, so the full report is below.

Five codons differ from table 4 with both methods agreeing:
- **ATA reads as Met** in *Ceraceosorus*.
- **AGA reads as Ala** in *Tremella*, *Apiotrichum*, *Tilletia* and *Trichosporon*.

ATA is probably Met in *Tremella* as well, and leans Met in *Apiotrichum*, but Codetta did not call either. No stop codon changes meaning in any genome.

**About the folder-name hint.** The name `blind_aga` points at AGA, so I treated all 64 codons the same way:
- **Method 1 (Codetta v2.0, run with -m)** infers every codon on its own.
- **Method 2 (ortholog-column tally)** ran one pre-planned Fisher test per codon per genome. The cut-off was Bonferroni over 512 tests: p below 0.05/512 = 9.8e-5.

ATA came up without any reference to the hint.

## Overall table

How to read it:
- The "strict" and "lenient" column sets are defined under Method 2 below.
- p = p_dev, a one-sided Fisher test: does this codon sit at columns of the odd residue more often than its table-4 synonyms do?
- Codetta logP is a natural log. Codetta only makes a call when the probability is above 0.9999.
- N = the number of aligned columns Codetta had for that codon.

| Genome | Inferred code | Method 1 (Codetta): call, N | Method 2: strict / lenient |
|---|---|---|---|
| BK072594.1 *Ceraceosorus guamensis* | **table 4 + ATA = Met** | ATA = M, N=62; logP(I) = −54.3 | M at 10/15 columns (p 8.8e-9); M at 18/28 (p 2.4e-16) |
| BK072956.1 *Exobasidium vaccinii* | table 4 | AGA = R, N=15; ATA ?, N=5 | AGA: R 7/7; R 11/12. ATA too rare |
| BK074202.1 *Solicoccozyma terrea* | table 4 (many codons too rare) | no deviations | nothing flagged; AGA leans R (4/4, too rare) |
| BK074411.1 *Tremella mesenterica* | **table 4 + AGA = Ala**; ATA = Met probable | AGA = A, N=40; logP(R) = −119.9. ATA ?, N=19: P(M) 0.84, P(I) 1e-6 | AGA: A 11/13 (p 8.7e-11); A 25/29 (p 1.3e-19). ATA: M 6/10 lenient (p 2.5e-7); strict 2/4, not significant |
| MT556837.1 *Meira* sp. | table 4 | AGA ?, N=6, P(R) 0.99 | AGA: R 4/4; R 5/5 |
| MZ439919.1 *Apiotrichum gamsii* | **table 4 + AGA = Ala**; ATA leans Met, not called | AGA = A, N=33; logP(R) = −86.2. ATA ?, N=31: P(M) 0.91, P(I) 0.09 | AGA: A 7/11 (p 5.3e-7); A 12/18 (p 1.0e-10). ATA lenient: L 7, M 6, I 1 of 19; planned test not significant |
| NC_009880.1 *Tilletia indica* | **table 4 + AGA = Ala** | AGA = A, N=24; logP(R) = −55.2. ATA = I, N=163 | AGA: A 5/6 (p 1.3e-6); A 11/18 (p 9.9e-10). ATA: I 21/28 strict |
| NC_054273.1 *Trichosporon inkin* | **table 4 + AGA = Ala** | AGA = A, N=27; logP(R) = −55.2 | AGA: A 5/6 (p 9.3e-7); A 7/12 (p 1.7e-7) |

**Stop codons and AGG:**
- TAA and TAG are too rare in every genome to test.
- TGA reads as Trp, as table 4 says, but this is confirmed only in *Ceraceosorus* (Codetta W, N=7; W at 9/9 strict columns). Elsewhere TGA is too rare.
- AGG can't be called anywhere: 0–4 copies in genes per genome, at most 1 on a conserved column. So I can't say whether AGG follows AGA to Ala.

## Methods

### Method 1: Codetta v2.0
- I used the named install and the exact command, run on copies of the genomes in `work/codetta/`.
- My comparison of its output against table 4 is in `work/codetta_vs_t4.tsv`.

### Method 2: ortholog-column tally

**Relatives.** 11 annotated basidiomycete mitogenomes, spread across orders:
- NC_003049.1 *Schizophyllum commune*
- NC_004336.1 *Cryptococcus neoformans* var. *grubii*
- NC_008368.1 *Ustilago maydis*
- NC_009905.1 *Pleurotus ostreatus*
- NC_020353.1 *Microbotryum lychnidis-dioicae*
- NC_021436.1 *Rhizoctonia solani* AG3
- NC_036340.1 *Rhodotorula mucilaginosa*
- NC_044103.1 *Puccinia triticina*
- NC_054271.1 *Trametes versicolor*
- NC_086480.1 *Malassezia dermatis*
- PX614074.1 *Wallemia sebi*

I re-translated their genes with table 4 and got exactly the GenBank proteins, with no internal stops. The genes used are cox1-3, cob, atp6/8/9, nad1-6, nad4L and rps3. Each gene had 8–11 relatives, aligned with MAFFT.

**Placing the target genes:**
- tblastn with table 4, e-value cut-off 1e-5.
- Hits were merged per gene, strand and frame; overlaps between genes were resolved by bitscore.
- Every codon was read through (stops become X), and the segments were added to the relatives' alignment with `mafft --addfragments`.
- All 15 genes were found in all 8 genomes: 15–34 segments and about 4,500–4,800 codons per genome.
- A coordinate check confirmed that all 37,000 tallied codons match the genome at their stated position and strand.

**Conserved columns:**
- **Strict:** at least 6 relatives have a residue and all of them agree.
- **Lenient:** at least 6 relatives, with the top residue at 80% or more.
- I excluded columns where 2 or more relatives start their protein, to keep start codons out.
- That leaves 1,275–1,426 strict and 2,238–2,489 lenient columns per genome.

**Planned tests**, for each codon X against its table-4 synonyms:
- **p_dev:** does X sit at columns of the top non-table-4 residue more often than its synonyms do?
- **p_lossA:** does X sit at columns of its own table-4 amino acid less often than its synonyms do?
- A codon is flagged as a DEVIATION when p_dev is below 9.8e-5 and the odd residue outnumbers the table-4 one.

**Added after I had seen the data (post hoc):**
- A pair test of ATA at Met columns against Ile columns, and of AGA at Ala against Arg, each compared with the synonyms.
- A fit of each codon's residue profile against every amino acid's profile, run on all codons.
- A minimum of 5 conserved columns before I will call a codon.

## Per-genome results

In the "too rare" lists, each codon is written as codon(N / n / in):
- N = Codetta aligned columns
- n = lenient conserved columns
- in = copies in the gene segments

**BK072594.1 *Ceraceosorus*: table 4 + ATA = Met**
- **Method 1:** Codetta calls ATA as M from 62 columns. logP is 0 for M, −41.7 for L and −54.3 for I.
- **Method 2:**
  - Strict: M at 10 of 15 columns (also L 3, V 1, I 1). The synonyms ATT+ATC sit at M 0/54 and I 46/54. p_dev 8.8e-9, p_lossA 2.7e-8.
  - Lenient: M at 18 of 28 (L 8, V 1, I 1), against synonyms at M 1/143. p 2.4e-16.
- **Where:** the M columns are internal, not start positions. They include six cox1 columns where all 11 relatives have M, plus columns in cob, nad4, nad5 and nad6.
- **How ATA is used:**
  - At Ile columns: strict ATA 1 against ATT 28 and ATC 18; lenient ATA 1 against ATT 64 and ATC 52. ATA has almost left Ile.
  - At Met columns: strict ATA 10 against ATG 3; lenient 18 against 20.
- **Pair test** (post hoc): p 2.9e-20.
- **Other codons:** AGA is too rare (3 in genes). TGA = W.
- **Too rare:** ACG(3/3/11), AGA(2/3/3), AGG(0/0/4), CCG(4/4/4), CGC(1/1/2), GTC(13/3/25), TAA(0/1/3), TAG(0/0/4), TGC(3/1/5), TTA(0/0/7), TTG(0/0/3).
- **Codetta uncalled, but method 2 agrees with table 4:** CTC, GCG, GGC, TCC, TCG.

**BK072956.1 *Exobasidium*: table 4**
- **AGA = Arg:**
  - Method 1: Codetta calls R from 15 columns; every other amino acid has logP of −25 or lower.
  - Method 2: R at 7/7 strict and 11/12 lenient columns. AGA sits at 11 Arg columns and no Ala columns.
- **ATA, unresolved:**
  - Codetta P(M) 0.52, P(I) 0.35, from only 5 columns.
  - Lenient: M 2 of 7, plus scattered residues. ATA occurs 25 times in genes.
  - Too rare to call.
- **TTG flag:** G at 2 of 3 columns, both at segment edges (likely misalignment at intron boundaries). Not called.
- **Too rare:** ACC(3/4/8), ACG(11/4/14), AGG(0/0/0), ATA(5/7/25), CAG(3/4/6), CGC(0/2/2), CTC(2/2/8), GTC(1/1/5), TAA(0/1/4), TAG(0/2/4), TCG(7/4/13), TGA(2/2/6), TGC(2/3/7), TTA(0/1/5), TTG(1/3/7).
- **Codetta uncalled, but method 2 agrees with table 4:** AAG, AGC, CCC, GAG, GCC, GGC, GTG, TCC.

**BK074202.1 *Solicoccozyma*: table 4, many codons too rare**
- No deviation from either method.
- AGA has 5 copies in genes and is R at 4/4 lenient columns (Codetta P(R) 0.95, N=2). It leans Arg but is too rare to call.
- **Too rare:** AAG(0/0/2), ACC(1/1/2), ACG(0/0/0), AGA(2/4/5), AGG(0/0/0), ATA(0/2/7), CAG(0/0/2), CCC(0/0/1), CGC(0/0/0), CTC(5/4/9), GAG(1/1/4), GGC(0/0/0), GTC(3/2/4), TAA(0/0/0), TAG(0/0/0), TCC(2/0/4), TCG(1/0/3), TGA(0/0/0), TGC(0/0/1), TTG(0/1/1).
- **Codetta uncalled, but method 2 agrees with table 4:** AGC, CGT, GCG, GTG.

**BK074411.1 *Tremella*: table 4 + AGA = Ala; ATA = Met probable**
- **AGA reads as Ala:**
  - Method 1: Codetta calls A from 40 columns. logP is 0 for A, −40.2 for S and −119.9 for R.
  - Method 2, strict: A at 11 of 13 (S 1, R 1). Synonyms CGN+AGG sit at A 0/52 and R 51/52. p_dev 8.7e-11, p_lossA 4.1e-11.
  - Method 2, lenient: A at 25 of 29, p 1.3e-19.
  - Where: spread across 11 genes and mostly internal; many are columns where all 11 relatives have Ala.
  - At Arg columns: AGA 1 against CGA 33, CGT 14 and CGG 4. At Ala columns: AGA 11 against GCN 56 (strict) and 25 of 141 (lenient).
  - Pair test: p 6.4e-22. AGA occurs 71 times in genes.
- **ATA reads as Met (probable):**
  - Method 1: Codetta P(M) 0.84 and P(I) 1e-6, from 19 columns. It doesn't make a call because the remaining 0.16 sits on its "?" (no-call) outcome, not on Ile.
  - Method 2, lenient: M at 6 of 10 (L 2, F 1, C 1, I 0), against synonyms at M 2/152. p_dev 2.5e-7, which passes the planned test.
  - Method 2, strict: M 2 of 4, not significant.
  - Where: internal columns in nad5 and cob (11/11 Met in both), plus nad1, nad4L and nad2.
  - At Ile columns: ATA 0 of 143 (ATT 93, ATC 34). At Met columns: ATA 6 against ATG 32.
  - Pair test (post hoc): p 3.7e-9.
- **TTG flag:** G at 2/2 columns, one at a segment edge. Not called.
- **Too rare:** AAG(3/2/3), AGC(4/1/14), AGG(0/1/1), CAG(2/2/2), CCC(3/4/9), CGC(0/0/1), GGC(0/1/2), TAA(2/0/6), TAG(0/0/0), TGA(0/1/1), TGC(1/0/0), TTG(1/2/3).
- **Codetta uncalled, but method 2 agrees with table 4:** ACG, GTC, TCC, TCG.

**MT556837.1 *Meira*: table 4**
- **AGA = Arg:**
  - Codetta P(R) 0.99 from 6 columns, just under its call threshold.
  - Method 2: R at 4/4 strict and 5/5 lenient. AGA sits at 5 Arg columns (CGA 40) and no Ala columns.
- **ATA:** too rare (5 in genes, 0 conserved columns).
- **TTA flag:** F at 3 of 4 internal columns, with Codetta N=2. Leu-for-Phe swaps are common at such columns, so I don't call it.
- **Too rare:** AAG(3/3/5), AGC(4/3/17), AGG(1/0/2), ATA(1/0/5), CAG(3/1/8), CGC(0/0/0), GGC(1/1/3), GTC(6/0/6), TAA(1/0/2), TAG(1/0/4), TGA(1/1/2), TGC(0/1/3), TTA(2/4/11), TTG(1/2/3).
- **Codetta uncalled, but method 2 agrees with table 4:** ACC, ACG, AGA, CTC, GAG, TCC.

**MZ439919.1 *Apiotrichum*: table 4 + AGA = Ala; ATA leans Met, not called**
- **AGA reads as Ala:**
  - Method 1: Codetta calls A from 33 columns. logP is −22.5 for "?", −34.8 for S and −86.2 for R.
  - Method 2, strict: A at 7 of 11 (S 2, G 1, R 1). Synonyms sit at A 0/53 and R 51/53. p_dev 5.3e-7, p_lossA 4.6e-9.
  - Method 2, lenient: A at 12 of 18, p 1.0e-10.
  - At Arg columns: AGA 1 against CGA 40, CGT 9 and CGG 2. At Ala columns: AGA 7 against GCN 61 (strict) and 12 of 141 (lenient).
  - Pair test: p 2.2e-13. AGA occurs 71 times in genes.
- **ATA:**
  - Method 1: Codetta P(M) 0.91, P(I) 0.09, from 31 columns.
  - Method 2, lenient: L 7, M 6, then I, F, T, S, D and V once each. The planned test picks the top non-Ile residue, which is L, and gives p 1.1e-3 (not significant). The 5 strict columns are scattered.
  - At Ile columns: ATA 1 of 140 (ATT 117). At Met columns: ATA 6 against ATG 33.
  - Post hoc: pair test p 7.8e-9. ATA's residue mix (M 32%, L 37%) looks like ATG's (M 65%, L 24%), not ATT's (I 82%); the profile fit prefers M over I by 30.7 nats.
  - Verdict: ATA is probably not Ile and is likely Met, but this is not called.
- **Too rare:** AAG(1/2/4), ACC(3/2/4), ACG(2/0/2), AGC(1/2/4), AGG(0/0/2), CAG(1/0/2), CCC(2/0/2), CGC(0/0/0), CGG(2/2/3), CTC(1/1/5), CTG(2/4/4), GAG(2/1/2), GCG(3/3/7), GGC(0/1/2), GTC(3/1/2), GTG(8/4/14), TAA(1/0/1), TAG(0/0/1), TCC(3/3/9), TCG(2/2/5), TGA(0/2/3), TGC(1/0/0), TTG(0/3/6).
- **Codetta uncalled, but method 2 agrees with table 4:** ATC, GGG.

**NC_009880.1 *Tilletia*: table 4 + AGA = Ala**
- **AGA reads as Ala:**
  - Method 1: Codetta calls A from 24 columns. logP is −10.4 for "?", −17.7 for G and −55.1 for R.
  - Method 2, strict: A at 5 of 6 (S 1). Synonyms sit at A 0/52 and R 51/52. p_dev 1.3e-6, p_lossA 1.7e-7.
  - Method 2, lenient: A at 11 of 18, p 9.9e-10.
  - At Arg columns: AGA 0 against CGA 30, CGT 15 and CGG 5 (strict); 1 of 71 (lenient). At Ala columns: AGA 5 against GCN 65 (strict) and 11 of 144 (lenient).
  - Pair test: p 1.3e-12. AGA occurs 76 times in genes.
- **ATA = Ile, firmly:**
  - Codetta calls I from 163 columns; logP(M) is −172.5.
  - Method 2: I at 21/28 strict columns.
  - ATA is the most common codon at Ile columns (62 of 146, against ATT 55).
- **Too rare:** ACC(1/2/6), AGG(0/0/0), CAG(2/2/4), CCG(2/3/5), CGC(1/2/3), GGC(1/2/3), GTC(9/4/13), TAA(0/0/2), TAG(0/1/1), TCG(5/4/11), TGA(1/1/3), TGC(0/1/1), TTA(0/1/6), TTG(0/1/1).
- **Codetta uncalled, but method 2 agrees with table 4:** AAG, ACG, AGC, CTC, GTG, TCC.

**NC_054273.1 *Trichosporon*: table 4 + AGA = Ala**
- **AGA reads as Ala:**
  - Method 1: Codetta calls A from 27 columns. logP is −13.8 for S, −14.3 for "?" and −55.2 for R.
  - Method 2, strict: A at 5 of 6 (S 1). Synonyms sit at A 0/56 and R 49/56. p_dev 9.3e-7, p_lossA 2.8e-5.
  - Method 2, lenient: A at 7 of 12, p 1.7e-7.
  - At Arg columns: AGA 0 against CGA 47 (strict); 1 of 73 (lenient). At Ala columns: AGA 5 against GCN 58 (strict) and 7 of 148 (lenient).
  - Pair test: p 5.4e-9. AGA occurs 55 times in genes.
- **ATA:** too rare (4 in genes).
- **Flags, not called:**
  - TTG reads as F at 3/3 internal columns (Codetta N=3). Probably Leu-for-Phe swaps.
  - TCG reads as W at 2/2 columns near segment edges (Codetta N=0).
- **Too rare:** AAG(1/2/3), ACC(0/0/1), ACG(0/1/2), AGC(0/2/5), AGG(0/0/0), ATA(0/2/4), CAG(3/3/4), CCC(1/1/4), CCG(0/1/1), CGC(0/0/0), CGG(1/1/1), CGT(4/4/6), CTC(0/1/2), GAG(0/2/3), GCG(2/3/6), GGC(1/2/2), GTC(5/2/8), GTG(9/2/15), TAA(0/1/7), TAG(0/0/0), TCC(2/3/6), TCG(0/2/3), TGA(1/3/3), TGC(0/1/1), TTG(3/3/3).
- **Codetta uncalled, but method 2 agrees with table 4:** CTG, GCC, GGG.

## Low-count flags

Five codons passed the planned Bonferroni test on only 2–4 columns:
- *Exobasidium* TTG → G, 2 of 3
- *Tremella* TTG → G, 2 of 2
- *Meira* TTA → F, 3 of 4
- *Trichosporon* TTG → F, 3 of 3
- *Trichosporon* TCG → W, 2 of 2

They pass only because the synonyms almost never carry the odd residue. Codetta had at most 3 columns for any of them. Several sit at segment edges, and the F cases look like ordinary Leu-for-Phe substitutions. I treat all five as too rare to call.

## Codetta code strings

64 codons in TCAG order; ? = not inferred.
```
BK072594.1  FF??S?S?YY??C?WWL?LLPPP?HHQQR?RRIIMMTTT?NNKKSS??V?VVAAA?DDEEG?GG
BK072956.1  FF??S?S?YY??C??WL?LLP?PPHHQ?R?RRII?MT?T?NNK?S?R?V?V?A?AADDE?G?GG
BK074202.1  FFL?S?S?YY??C??WL?LLP?PPHHQ???RRII?MT?T?NNK?S???V?V?AAA?DDE?G?GG
BK074411.1  FFL?S?S?YY??C??WLLLLP?PPHHQ?R?RRII?MTTT?NNK?S?A?V?VVAAAADDEEG?GG
MT556837.1  FF??S?SSYY??C??WL?LLPPPPHHQ?R?RRII?MT?T?NNK?S???V?VVAAAADDE?G?GG
MZ439919.1  FFL?S?S?YY??C??WL?L?P?PPHHQ?R?R?I??MT?T?NNK?S?A?V?V?AAA?DDE?G?G?
NC_009880.1 FF??S?S?YY??C??WL?LLPPP?HHQ?R?RRIIIMT?T?NNK?S?A?V?V?AAAADDEEG?GG
NC_054273.1 FFL?S?S?YY??C??WL?L?P?P?HHQ???R?II?MT?T?NNK?S?A?V?V?A?A?DDE?G?G?
```

## Limitations and what I could not do

- **report.md was not written**, because the harness blocks report files from subagents. The content above is the report.
- **System python slip.** One early Bash command accidentally called the system `python3 -` with empty input, next to a perl edit. It ran nothing and touched no data, but it breaks the uv-only rule. Every analysis ran through `uv run --no-project ... python -I`.
- **Post hoc tests.** The ATA Met-vs-Ile pair test, the profile fit and the 5-column minimum were chosen after I had seen the data.
  - The firm calls (*Ceraceosorus* ATA → Met and the four AGA → Ala) don't depend on them: each passes both the planned test and a Codetta call.
  - The *Tremella* and *Apiotrichum* ATA verdicts do depend on them, which is why I graded those "probable" and "leans".
- **Relatives are annotated with table 4.** If a relative really reads AGA as Ala or ATA as Met, its GenBank protein shows R or I there. That makes the method conservative, not prone to false calls.
- **No shuffled or decoy null** was built for the ortholog tally. The table-4 synonyms are the only control.
- **Alignment caveats.** Some flags sit at segment edges, where intron-split genes can misalign. The lenient consensus lets common substitute residues through (L/M/I/V, F/L).
- **Uncallable everywhere:** AGG, TAA and TAG. TGA is callable only in *Ceraceosorus*.
- **Not attempted:**
  - tRNA anticodon analysis.
  - Any comparison with the literature.
  - An independent subagent review of the scoring code. I did check every tallied codon against the genome coordinates.

Everything is in `scratch/blind_aga/`. Scripts are in `src/`. The files to open first are `work/codetta_vs_t4.tsv`, `work/ortho/all.codon_table.tsv` and `work/merged_calls.tsv`.
