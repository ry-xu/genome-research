# RULES — proline-codon yeasts (Saiph, frozen 20 Nov before any target analysis)

Target: the post `outbox/2026-11-11-mito-proline-codons/` (Clover; work in `work/new-stops-euk/`).
Claim: in the mitochondria of *Wickerhamiella siamensis* BK074538, *W. infanticola* BK074532, *W. verensis* BK074542 and
*Sugiyamaella pinicola* BK074262, the codons CUU/CUA are read as Pro.
Task (BRIEF 20 Nov): (T) find the tRNA that reads CUU/CUA as Pro, (R) search all public RNA to test U→C editing,
then either get the post to "needs image" or kill it.

## Seen before freezing (disclosure)
- Clover's README: all four genomes keep tRNA-Pro(UGG), tRNA-Leu(UAA) and a 72–73 nt tRNA(UAG). Most Dipodascomycetes
  also have a UAG tRNA, including the CUN = Leu genome BK074539. The *W. siamensis* UAG tRNA is 0.67 identical to its own
  tRNA-Pro and is no closer to it than to other tRNAs.
- `aga_ala/fungi.aragorn.tsv`: the *Ashbya gossypii* (NC_005789.1) UAG tRNA has G:U at 3:70 and A73.
- ENA read_run, tax_tree(*Wickerhamiella*) and tax_tree(*Sugiyamaella*) with library_source != GENOMIC: 29 runs.
  All are RNA-Seq of *W. versatilis*, *W. domercqiae* or *S. americana* (all CCN users in Clover's grid), plus one
  16S amplicon run. There are none for the four target species.
- Tree (`wickerhamiella_tree/wt.treefile`): the three target *Wickerhamiella* form a clade. Their next clade adds
  BK074528, BK074537 and BK074540. *S. pinicola*'s smallest clade adds BK074255, BK074256, BK074258, BK074260 and BK074265.

## Arm T — which tRNA reads CUU/CUA as Pro

**T0. Calls**
- pyaragorn, with the code of `new-stops-euk/src/codes_census/org_trnas.py`, run on all 49 Dipodascales mitogenomes
  in `wickerhamiella_tree/accs.txt` and on the control genomes. Genomes are read as circular.
- The output adds each tRNA's anticodon position.

**T0b. Missed-tRNA sweep**
- Query: all 49 genomes' tRNA-Pro and tRNA-Leu sequences. Search each target genome with
  `blastn -task blastn-short` (or word 7) at e ≤ 1e-3.
- Any hit that overlaps no ARAGORN call by ≥ 30 nt is a candidate. Report its anticodon (the triplet at the aligned
  34–36) and its ARAGORN call when ARAGORN is re-run on hit ±50 nt with relaxed settings.

**T1. Origin of the CUN tRNA (S2, the main statistic)**
- Query q is each target genome's UAG tRNA (and any T0b candidate).
- Sisters are the CCN-using genomes named in the disclosure above.
- D = mean id(q, sisters' UAG tRNAs) − mean id(q, sisters' tRNA-Pro(UGG)).
  - id = `trna_origin.py`'s global identity (match 1, mismatch −1, gap −2/−1), with the anticodon triplet set to N in
    every sequence.
- D > 0.05: q is the ortholog of the relatives' CUN tRNA, recharged.
- D < −0.05: q is Pro-derived.
- Otherwise: unresolved.
- Calibration (same-shape decoys): the same D for every other tRNA X in the target genome:
  D_X = mean id(X, sisters' same-anticodon X) − mean id(X, sisters' tRNA-Pro).
  - It shows how far real orthologs sit from Pro.
  - It is reported as a distribution; no p-value is claimed.

**S1. Own-genome nearest class**
- The class of q's best-matching other tRNA in the same genome, with the anticodon masked as above.
- Leu/UAG tRNAs are excluded from the candidates.

**S3. Identity bases**
- Positions read: acceptor-stem pairs 1:72–7:66, discriminator 73, and the 7-nt anticodon loop (32–38, located from
  the anticodon position).
- "Pro-specific" bases are computed from non-UAG tRNAs of the 49 genomes only. A base qualifies if:
  - ≥ 90% of tRNA-Pro(UGG) carry it, and
  - < 20% of all other tRNAs carry it.
- Count the Pro-specific bases in each q against its sisters' UAG tRNAs.
- The same count for every other amino acid's specific-base set is the decoy: the q should not gain those.
- **Convergence** = the same Pro-specific base gained in both origins (the *Wickerhamiella* clade and *S. pinicola*).

**Positive controls (same code, run before reading any target output)**
- **PC-T1, *Ashbya gossypii* NC_005789.1**
  - CUU/CUA = Ala, via a tRNA(UAG) derived from the yeast-mito tRNA-Thr(UAG) (Ling et al. 2014, NAR 42:499).
  - Sisters: the nearest CUN = Thr Saccharomycetaceae mitogenomes in `aga_ala/fungi.aragorn.tsv` that have a
    tRNA(UAG), up to 5. They are chosen by the census's table assignment and family, before any identity is computed.
  - Pass:
    - D (vs sisters' UAG tRNA-Thr, against sisters' tRNA-Ala(UGC)) > 0.05.
    - S3 shows G3:U70, absent from all the sisters' UAG tRNAs.
- **PC-T2, *Saccharomyces cerevisiae* mito tRNA-Thr(UAG)**
  - Derived from tRNA-His (Su et al. 2011, NAR 39:4866).
  - Pass: S1 nearest class = His.
  - If PC-T2 fails, S1 is uninformative and is dropped, and the post says so.
- If PC-T1 fails, T1/S3 cannot speak, and I post that before reading the targets.

## Arm R — RNA and U→C editing

**R0. Inventory (all logged in `data/inventory/`)**
- **ENA read_run**, tax_tree(410830, Trichomonascaceae), with library_source in {TRANSCRIPTOMIC, METATRANSCRIPTOMIC,
  TRANSCRIPTOMIC SINGLE CELL, OTHER}.
- **NCBI SRA E-utilities**: `txid410830[Organism:exp] AND (rna[Filter] OR transcriptomic)`. Sent with no email or
  personal parameter.
- **Y1000+ transcriptome BioProjects (PRJNA372xxx series)**: list every species, to catch target species filed under
  an older name.
- **PRIDE and MassIVE**: keyword search for *Wickerhamiella*, *Sugiyamaella*, *Starmerella*, *Candida infanticola*,
  *C. sorbophila* and *C. siamensis*.
- **Metatranscriptomes**:
  - Logan Search, at most 8 manual queries (1-kb cox1/cob windows of the target genomes), only if usable without an
    account. Otherwise skipped and stated.
  - The NCBI STAT parquet, only if a tax_id-filtered read of the 4 target taxids costs < 5 GB of transfer.
  - Hits must pass a read-level species check: ≥ 95% identity over ≥ 100 nt to the target mitogenome, and better than
    to any of the 45 other Dipodascales mitogenomes.

**R1. Editing test in a target species (only if R0 finds RNA)**
- Map reads with `minimap2 -ax sr` to the target mitogenome plus the 45 others as decoys. Keep MAPQ ≥ 20.
- Read the base at the 2nd position of every CUU/CUA codon at a Pro column (Clover's `cun_CTT.tsv`/`cun_CTA.tsv`
  sites) and at every CCN Pro and UUA/UUG Leu codon.
- Editing call: median C fraction ≥ 0.5 at CUU/CUA-Pro sites with depth ≥ 10.
- No editing: median ≤ 0.1, with the same depth gate and ≥ 20 sites.

**R2. Does the family edit at all? (relatives' RNA)**
- Runs: up to 2 per species of *W. domercqiae* (Y1000 NRRL Y-6692 and PRJNA1146581), *W. versatilis* and
  *S. americana* NRRL YB-2067, each against its own mitogenome (BK074530, BK074543, BK074251).
- Each run is streamed with `fastq-dump -X 4000000 --split-spot --stdout`. Before freezing a run, check that it holds
  the species: ≥ 1% of reads map to the mitogenome at ≥ 95% identity.
- Caller:
  - Positions with depth ≥ 20, base quality ≥ 20, MAPQ ≥ 20.
  - An RNA–DNA difference (RDD) is a position where ≥ 20% of reads disagree with the assembly.
  - A WGS run of the same strain, put through the same caller, removes assembly and variant sites. If no same-strain
    WGS run exists, that is stated.
  - All 12 substitution types are tabulated. Strand comes from the gene's orientation (genes from Clover's census
    coordinates, or blastx of Sc mito proteins when absent).
- Result: the count of U→C (and C→U) RDDs inside coding sequence, by codon position.

**PC-R (editing-caller control, same code)**
- An *Arabidopsis thaliana* Col-0 Illumina RNA-Seq run against NC_037304.1.
  - The run is chosen as the first ENA run matching: Col-0, total or rRNA-depleted RNA, ≥ 10 M spots.
  - Its gene strands come from the RefSeq annotation.
- Pass:
  - ≥ 100 C→U RDDs in coding sequence.
  - ≥ 80% of all coding RDDs are C→U.
  - The majority are at codon positions 1–2.
- If PC-R fails, R2 cannot speak, and I say so.

## Verdict rules (fixed now)
- **KILL** if either:
  - R1 calls editing in any target species; or
  - a prior-art pass finds the reassignment already published (search 2025–26 literature and GitHub, including Sarhan
    et al. 2025 and its supplements).
- **UNRESOLVED, stays "not ready"** if R2 finds ≥ 10 U→C coding RDDs in relatives' mito mRNA while PC-R passes. The
  family edits, so editing becomes a live explanation.
- **"Needs image"** if all of the following hold:
  - no target RNA calls editing (or no target RNA exists);
  - R2 finds < 10 U→C coding RDDs in all relatives, with PC-R passing;
  - Arm T is reported, whatever it shows;
  - a fresh subagent (raw data) and Sabik have checked the post edits.
- **Wording when no target RNA exists**: "editing not excluded directly; no U→C editing in the transcripts of three
  related yeasts". A request for Ryan goes in `outbox/` (RNA or proteome of CBS 13331 / NRRL Y-17858).

## Machine
- One heavy job at a time: 2 threads if the 1-min load is < 8, 1 if < 12, wait if > 12.
- ≤ 4 GB of downloads planned (cap 10). Streams are deleted after counting.

## Amendment A1 (10:45 20 Nov; before any R2/relative RNA read was mapped)
- PC-R run DRR917134 (first ENA match) gave 4,019 MAPQ>=20 mito reads of 8 M; 0 coding positions at depth >= 20.
  Its first 200k spots put 0.03% of reads on the chloroplast and 0.04% on the mitogenome: an organelle-poor
  (in effect polyA) library despite "Inverse rRNA" metadata. The control was empty, not failed.
- New rule: walk the same ENA list (data/inventory/ath_inverse_rrna.tsv, file order, Col-0 rows) and take the first run
  whose first 200k spots put >= 1% of reads (MAPQ >= 20) on mito + chloroplast (NC_037304.1 + NC_000932.1).
  Mapping for the caller stays against NC_037304.1 alone, as frozen. Caller code and pass bars unchanged.
- The same input pre-check applies to R2 runs already (species check, >= 1% at >= 95% identity).

## Amendment A2 (11:07 20 Nov; before the S. americana and W. versatilis PYCC 7102 calls were read)
- A relative's run counts toward "< 10 U→C coding RDDs" only if it passes the species check AND has >= 1,000 coding
  positions at depth >= 20, of which >= 100 are U (transcript strand) at codon position 2. Runs below this are
  reported with their numbers as "too little mito RNA to test". (PC-R had 29,613 positions.)
- If no relative's run passes, R2 cannot speak, and the post wording says so.
