# Blinded mitochondrial genetic-code calls (R01–R09)

## Calling rule (frozen 2026-11-11 15:06 EST (file birth time), before any per-genome result was looked at)

Codon meanings are compared with NCBI table 1 (standard code). "std" = the table-1 meaning of the codon.

### Method 1: Codetta (Pfam-A_enone, `-m`, defaults otherwise)
Codetta v2 run per genome with the Pfam database, `-m` (mitochondrial Pfams included), e-value 1e-10,
decoding probability threshold 0.9999, max fraction 0.01. Codetta "infers X" for a codon when its output letter
is X (not `?`). Score reported: N consensus columns used, and the log-probability gap
logP(X) − logP(std) (or logP(X) − logP(2nd best) for table-1 stop codons) from the `.genetic_code.out` file.

### Method 2: reference-protein alignment, per-codon likelihood
1. Six-frame translation of each genome with every table-1 stop codon (and any codon with a non-ACGT base)
   written as `X`, so reassigned stops do not break reading frames.
2. Reference mitochondrially encoded proteins (cox1/2/3, cob, atp6/8/9, nad1–6, nad4L, ribosomal proteins etc.)
   fetched from UniProt (curated Swiss-Prot, taxonomically diverse across eukaryotes, subsampled per gene and
   lineage). Searched against the six-frame translations with BLASTP (e-value ≤ 1e-5).
3. Each aligned genome codon is a "site". Sites in the 5 aligned columns at each end of an HSP are dropped. Where
   one reference has several HSPs over a site, its best HSP is kept. Each site's reference residues get total
   weight 1 (split evenly over the references aligned there; reference gaps ignored).
4. For each amino acid a, a profile P_a(b) (reference residue b given genome amino acid a) is estimated from this
   genome's own sites of codons whose table-1 meaning is a, leaving out the codon under test (pseudocount 1 per
   residue). Pass 1 uses table-1 meanings; codons that reach "M2-strong" in pass 1 are moved to their inferred amino
   acid and the profiles are re-estimated (pass 2). Pass 2 numbers are reported.
5. For codon c with sites S_c: LL(a) = Σ_s Σ_b w_s(b) · ln P_a(b), for all 20 amino acids. Reported: n = |S_c|,
   best amino acid X, LLR_std = LL(X) − LL(std) (nats), LLR_2nd = LL(X) − LL(second best).

Method-2 levels for a codon whose best amino acid X differs from std:
- **M2-strong(X)**: n ≥ 10, LLR_std ≥ 10, LLR_2nd ≥ 3.
- **M2-support(X)**: n ≥ 5, LLR_std ≥ 3.
For a table-1 stop codon (no std amino acid) the internal aligned sites are themselves evidence of sense use:
- **M2-strong(X)**: n ≥ 10 and LLR_2nd ≥ 3.  **M2-support(X)**: n ≥ 5 and LLR_2nd ≥ 1.

**Contradiction** of X: Codetta infers a different letter than X (including the std amino acid), or Method 2's best
amino acid is std / another amino acid with n ≥ 5 and LLR over X ≥ 3.

### Tiers (one per codon per genome)
- **call**: no contradiction, and either (a) M2-strong(X) with Codetta inferring X or `?`, or
  (b) Codetta infers X and M2-support(X).
- **weak, not a call**: evidence for a non-standard X that does not reach "call": Codetta infers X alone
  (Method 2 below support), M2-support(X) alone, or the two methods conflict. Also: a sense codon with ≤ 1
  internal aligned site while the other table-1 codons of its amino acid have ≥ 100 sites in that genome
  (possible stop/unassigned codon) — flagged weak, never a call.
- **no call**: data too thin for both methods (Codetta `?` and Method 2 n < 5). Listed for the three table-1
  stop codons and for any codon with weak/conflicting evidence; codons not listed are consistent with table 1
  or uninformative.

Control (run before final calls, thresholds not adjusted by it): Method 2 and Codetta on the codetta example
mitogenome `examples/NC_013935.1.fna` (README states UGA = Trp, AUA = Met); expect those two and nothing else.

### Rule change v1 → v2 (2026-11-11 ≈15:09 EST, after the control run, before any R01–R09 result was looked at)
The control exposed a bug in Method 2 step 4 (v1). Leave-one-codon-out empties the profile of an amino acid that
has a single table-1 codon (Met = ATG, Trp = TGG), so in pass 1 ATG scored as Ile (M2-strong, wrong), was moved to
Ile, and contaminated the Ile profile, which then gave a false pass-2 ATT → Val (n = 208, LLR_std 14.3).
v1 control result: TGA → W (strong), plus false ATG → I (pass 1) and ATT → V (pass 2); ATA → M was found in pass 1.
**v2 (used for all reported results):** the codon under test is left out of its own amino acid's profile only
if ≥ 30 sites remain for that amino acid; otherwise it stays in. Thresholds and tiers are unchanged.
v2 control result: TGA → W (n = 37, LLR_2nd 61.2, strong), ATA → M (n = 76, LLR_std 16.7, LLR_2nd 12.5, strong);
no other codon has a non-standard best amino acid; TTA and TTG have 0 aligned sites (absent flag → weak only).

### Post-hoc change v2 → v3 (2026-11-11 ≈15:12 EST, AFTER looking at R01–R09 results)
Looking at the v2 results showed one failure mode: in R05, TAA had 33 Method-2 "sites" that are BLAST extensions
through poorly aligned stretches (mostly 1–2 references per site, local identity often < 0.2; Codetta aligned
0 columns at TAA). Diffuse hydrophobic noise fits the Met profile best, so v2 gave TAA → Met as M2-strong, and
moving TAA into Met in pass 1 then pushed ATG → Leu in pass 2.
**v3 (post-hoc):** a Method-2 site is kept only if ≥ 3 references align to it and the median local identity
(± 10 aligned columns around the site, X columns ignored) is ≥ 0.30. Everything else unchanged.
v3 control result: TGA → W (n = 34, LLR_2nd 67.5), ATA → M (n = 66, LLR_std 17.8); nothing else non-standard.
**How the two versions are combined:** a codon is a final "call" only if it is a call under both the frozen v2 rule
and v3. A call under only one version is reported as "weak, not a call". The post-hoc change can therefore
only remove calls, never add them. (It removed exactly one: R05 TAA → Met.)

## Codon table

Scores: **M2** = Method 2 (n sites; LLR_std = ln-likelihood of the called amino acid minus table-1 amino acid;
LLR_2nd = minus the runner-up; v2 numbers, v3 in brackets). **CD** = Codetta letter, consensus columns used, and
log-probability gap to the table-1 amino acid / to the runner-up. Codons not listed agree with table 1 or are uninformative.
For table-1 stop codons, "No call" with n ≥ 5 means aligned sites exist but no amino acid is preferred
(LLR_2nd < 1), as expected for stray alignment through a real stop; nothing indicates sense use.

| Genome | Codons read differently from table 1 (tier: evidence) |
|---|---|
| R01 | **AGG → Met (call)**: M2 n=105, LLR_std 75.3, LLR_2nd 26.5 vs Leu [v3 n=76, 72.3, 21.3]; CD M, 48 cols, gap 158.8 vs Arg, 40.1 vs Leu. **CGG → Leu (call)**: M2 n=31, LLR_std 18.4, LLR_2nd 3.5 vs Met [v3 n=23, 24.2, 4.8]; CD `?` (17 cols, best L, gap 66.1 vs Arg, 9.8 vs Ile; below Codetta's 0.9999 threshold). *Weak:* TAA → Ala (v2 only: n=9, LLR_2nd 2.2; v3 n=2; CD 0 cols, so noise). *Weak, possible stop/unassigned (absent from gene cores):* TCA (v2 n=2, v3 n=1; first stop-like codon after 4/12 genes vs TAA 5/12; 0–3/13 in other genomes), TCG (v3 n=0), TTA (v3 n=1), AGA (v3 n=1), GGG (v2 n=1, v3 n=0). *No call:* TAG, TGA (M2 n ≤ 3, CD 0 cols; no sense use seen, so consistent with stop). |
| R02 | **TGA → Trp (call)**: M2 n=88, LLR_2nd 100.0 vs Phe [v3 n=63, 102.9]; CD W, 42 cols, gap 110.7. **AGG → Met (call)**: M2 n=50, LLR_std 72.3, LLR_2nd 15.7 vs Leu [v3 n=41, 67.9, 15.7]; CD M, 25 cols, gap 68.8 vs Arg, 12.6 vs Leu. *Weak, possible stop/unassigned:* CGG (v2 n=2, v3 n=0). *No call:* TAA (v2 n=5, diffuse, LLR_2nd 0.5; v3 0), TAG (0). |
| R03 | **TGA → Trp (call)**: M2 n=78, LLR_2nd 95.5 vs Phe [v3 n=54, 98.2]; CD W, 34 cols, gap 85.5. *No call:* TAA (v2 n=18 and TAG v2 n=12, both diffuse, LLR_2nd < 1; v3 n=1 each; CD 0 cols; consistent with stop). |
| R04 | **TAG → Leu (call)**: M2 n=97, LLR_2nd 23.8 vs Ile [v3 n=77, 25.2 vs Met]; CD L, 60 cols, gap 33.8 vs Ile. **AGG → Ala (call)**: M2 n=92, LLR_std 72.5, LLR_2nd 7.0 vs Cys [v3 n=62, 66.7, 12.9 vs Ser]; CD A, 40 cols, gap 100.4 vs Arg, 15.4 vs Ser. *Weak, possible stop:* TCA (v2 n=9 low-quality sites, best Thr, LLR_std 1.7; v3 n=0; CD 0 cols; 400 TCA trinucleotides in the genome; TCA is the first stop-like codon after 10/12 genes' aligned ends, vs 0–3/13 in the other genomes). *Weak, absent from gene cores (no terminator signal):* TTA, TCG, CTG, CCG, CGG, ACG, AGA (n ≤ 1 in v2 and v3), ATA (v2 n=2, v3 n=0). *No call:* TAA (n=2/0), TGA (n=3/1); CD 0 cols. |
| R05 | **TGA → Trp (call)**: M2 n=58, LLR_2nd 73.4 vs Met [v3 n=34, 58.2]; CD W, 18 cols, gap 42.8. *Weak, not a call:* TAA → Met (call under frozen v2: n=33, LLR_2nd 4.0 vs Tyr, CD `?` with 0 cols; under v3 n=5, best Leu, LLR_2nd 1.1). ATG → Leu (v2 only, n=81, LLR_std 11.8, contradicted by CD M with 31 cols; a side effect of TAA being moved into Met; v3 ATG = Met). *No call:* TAG (v2 n=19 diffuse, LLR_2nd 0.0; v3 n=4). Data thin: CGG (n=4, CD `?`). |
| R06 | **TGA → Trp (call)**: M2 n=86, LLR_2nd 119.1 vs Tyr [v3 n=64, 111.5]; CD W, 46 cols, gap 118.1. **AGG → Met (call)**: M2 n=52, LLR_std 73.5, LLR_2nd 17.0 vs Leu [v3 n=40, 58.5, 15.6]; CD M, 26 cols, gap 85.7 vs Arg, 19.2 vs Leu. *Weak, possible stop/unassigned:* CGG (n=0 in v2 and v3; CD 0 cols). *No call:* TAA (n=5 diffuse/0), TAG (2/0). |
| R07 | **TGA → Trp (call)**: M2 n=75, LLR_2nd 98.5 vs Tyr [v3 n=61, 103.7]; CD W, 40 cols, gap 100.1. *No call:* TAA (n=3), TAG (n=1). Data thin: CGG (n=2). |
| R08 | **TGA → Trp (call)**: M2 n=77, LLR_2nd 103.4 vs Phe [v3 n=64, 99.0]; CD W, 40 cols, gap 53.9. *No call:* TAA (v2 n=11 diffuse, LLR_2nd 0.2; v3 0), TAG (2/0). Data thin: CTC, CCG, CGC, CGG, GTC, GAC. |
| R09 | **TGA → Trp (call)**: M2 n=76, LLR_2nd 113.1 vs Phe [v3 n=64, 107.7]; CD W, 43 cols, gap 101.2. *No call:* TAA (v2 n=9 diffuse; v3 n=4), TAG (2/0). Data thin: CTC, CTG, CAC, CGG. |

## Frozen calls (final)

These calls are final. Each is a call under both the frozen v2 rule and v3. Codetta and Method 2 infer the same amino acid for each, except R01 CGG, where Codetta outputs `?` but its top amino acid is also Leu.

| Genome | Final calls |
|---|---|
| R01 | AGG = Met; CGG = Leu (the only call made by Method 2 alone; Codetta's top amino acid is also Leu but below its threshold) |
| R02 | TGA = Trp; AGG = Met |
| R03 | TGA = Trp |
| R04 | TAG = Leu; AGG = Ala |
| R05 | TGA = Trp |
| R06 | TGA = Trp; AGG = Met |
| R07 | TGA = Trp |
| R08 | TGA = Trp |
| R09 | TGA = Trp |

Not final: R05 TAA → Met was a call under the frozen v2 rule only. It is downgraded to "weak, not a call" by the
post-hoc v3 filter, which is disclosed above. All codons not listed are consistent with table 1 where there are data.
No sense-to-stop change is called anywhere: the frozen rule has no call path for it. The strongest such signal is
R04 TCA (see the table).

## Method notes

- **Codetta v2.0**, Pfam-A_enone, `-m`, other settings default (e-value 1e-10, P ≥ 0.9999, max fraction 0.01). It was
  run as `codetta_align` → `codetta_summary` → `codetta_infer`, with `HMMER_NCPU=1` and 4 genomes at a time. The
  installed copy crashes under numpy 2 (`dtype 'a1'`), so summary/infer ran with `uv --python 3.12 --with 'numpy<2'
  --with 'scipy<1.14'`. Codetta gave a letter for 36–50 of 64 codons per genome. It aligned 0 columns at TAA in
  every genome. It also aligned 0 columns at TAG, and at TGA, wherever those were not reassigned.
- **Method 2.** References: 603 mitochondrion-encoded Swiss-Prot proteins (UniProt REST, `organelle:mitochondrion
  AND reviewed:true`, eukaryotes, no fragments, ≥ 40 aa). Up to 3 Vertebrata, 3 Arthropoda, 6 other Metazoa,
  10 Fungi, 8 Viridiplantae and 12 other eukaryotes per gene family, one per genus. Families: cox1–3, cob, atp1/6/8/9,
  nad1–6, nad4L, mito ribosomal proteins. Each genome was translated in 6 frames with table-1 stops as X, and
  searched with BLAST+ blastp (refs as queries, e-value ≤ 1e-5, max 20 HSPs per target). Every genome had hits to the 13
  core families (cox1–3, cob, atp6, atp9, nad1–6, nad4L); R03 also hit rps3. That gave 3,900–4,300 sites per genome
  (2,600–3,400 after the v3 filter). The likelihood was computed per codon from profiles estimated on the
  genome's own sites, in 2 passes, as in the rule. Scripts are in
  `scratchpad/blindagg.4ceG/scripts/`
  (`m2.py`, `tiers.py`; `m2_v1.py` and `m2_v2.py` are the earlier versions).
- **Control** (Codetta's example NC_013935.1; the README states UGA = Trp and AUA = Met): Codetta and Method 2 v2/v3
  both returned exactly TGA = Trp and ATA = Met as calls. TTA and TTG had 0 aligned sites (absent flag, weak only).
- **Terminator scan (post-hoc diagnostic, not part of the rule).** For each gene family, take the best HSP that
  reaches within 20 aa of the reference C-terminus. Report the first table-1 stop or flagged codon within 60
  codons downstream, in frame. As a null, in R02, R03 and R05–R09 TCA came first after 0–3 of 12–13 genes (TCG after 0–2).
- **Blinding.** No genome or derived sequence was searched against any public database. Only reference proteins
  were downloaded, and they were searched against the genomes. No species was looked up or inferred.
- **Data too thin to call** (Method 2 n < 5 and Codetta `?`, table-1 sense codons):
  - R01: TCA TCG CCG CAG AAG GAG GGG
  - R02: CTG CGG
  - R03: none
  - R04: TTA TCG CTG CCG CGG ATA ACG AGA
  - R05: CGG
  - R06: CTC CTG CGC CGG GAC GAG
  - R07: CGG
  - R08: CTC CCG CGC CGG GTC GAC
  - R09: CTC CTG CAC CGG

  Also thin everywhere: whichever of TAA/TAG/TGA is not reassigned in a genome (0 Codetta columns). ATA, AAA and
  AGA are table-1 meanings wherever there are data. The exceptions are R04 (ATA, AGA thin) and R01 (ATA n=11, best
  Ile, CD `?`; AGA thin).
