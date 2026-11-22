<!-- Copied from the lab record (work/silent-operon/, 22 Nov 2026). Paths such as out/, src/, recheck/ or frozen-core/ refer to the original working folder; the small tables from out/ are in this folder (see ../README.md), the rest is not included. -->
# C_rna rules: is "silent" true beyond three RNA samples from one species? (frozen 22 Nov 2026, before any new read is mapped)

## What the search found (metadata only, no reads looked at)
- All Crystal Geyser RNA in SRA (`search/cg_all.csv`, 729 CG runs): 8 JGI metatranscriptomes + 1 dsRNA run.
  - Used before: SRR11558543 (CG04), SRR11558542 (CG08), SRR11558425 (CG16).
  - New, same 2015 series: **SRR11558439 (CG05, 25 May 2015)**.
  - New, older JGI series "Crystal Geyser CG12/CG13/CG22 metaT" (PRJNA416596-8, Gp0115573-5, date missing):
    **SRR6261398 (CG12), SRR6241672 (CG13, selection RT-PCR), SRR6241673 (CG13, selection "other"), SRR6261448 (CG22)**.
  - **SRR7589086**: dsRNA purified from Crystal Geyser Site 5 mineral matrix (Feb 2015, PowerSoil RNA kit). Secondary.
- Esser et al. 2025 (Environ Microbiol Rep 17:e70096) used CG05, CG08, CG16 as *cellular* metatranscriptomes: "the erupted
  water from CG was sequentially filtered on PTFE filters, which were frozen on site on dry ice". So the "viral communities"
  title may be a GOLD label, not a viral-fraction prep. Which filter pore size was used is not stated there.
- Other carriers: only the CADDYS01 source project (PRJNA666025, China acid-mine-drainage sediment) has RNA:
  **SRR33092280 (Dabaoshan, Guangdong), SRR33092281 (Dexing, Jiangxi), SRR33092282 (Anhui)**, 2023, METATRANSCRIPTOMIC,
  cDNA. Their DNA source runs are from Dashan/Wuyi (Guangxi). Branchwater finds CADDYS01 sp021794835 at cANI 0.98-0.99 in
  Guangdong AMD sediment (PRJNA910809), so presence in at least SRR33092280 is plausible but unknown.
- No RNA in the source projects of DASQQM01 (PRJNA317932), DASXVE01 (PRJNA620721), JASEIR01 (PRJNA983765),
  Kmv33 (PRJNA636628) or CADDYS01 sp902810715 (PRJEB36447).

## Data
- `mags fetch sra RUN --spots N` (streams the first N spots; whole-run download of 20-40 Gbp runs would break the 5 GB cap).
  These runs are unaligned, so the first N spots are sequencing order, not genome order. Both loci of a species are hit
  alike, so ratios are unaffected.
- Pilot: 1 M spots of each of the 10 new runs above. Then extend (up to the 5 GB total) only the runs whose pilot gives
  >= 3 normal-23S pairs for a target species, aiming at >= 30 normal-23S pairs per run. The extension rule uses only the
  normal counts (the divergent count is not looked at to decide depth).

## Reference and counting (the earlier method, unchanged)
- Reference: `orthogonal-ribosome/out/t5/ref.fa` (20,795 targets: the family 23S loci + 17,329 census 23S decoys; family
  16S loci with 150 nt flanks + 16S decoys). Both copies and decoys in one reference.
- `minimap2 -x sr -t 1 -c --secondary=no`; a read counts if primary, >= 60 aligned columns, identity (PAF col10/col11)
  >= 0.98. A read pair counts once, at the locus of its higher-AS mate.
- Targets (species groups):
  - **AS** = *Ca.* Aquicultor secundus, the 8 Crystal Geyser genomes: 23S classes from
    `frozen-core/.../expr/loci_class_all.tsv`, 16S classes from `orthogonal-ribosome/out/t5/class16.tsv` (as T5).
  - **CD** = CADDYS01 sp021794835 (GCA_027345605.1, GCA_027392615.1, GCA_021794835.1): 23S div = DAHYRC010000072
    (full length); 23S norm = DAHWGC010000001, DAHYRC010000115, JAJYQA010000029 (identical). 16S div = DAHYRC010000072
    16S (same operon as the div 23S); 16S norm = DAHWGC010000001 and JAJYQA010000029 16S; the DAHYRC010000125 16S fragment
    is classed by >= 97% identity to either (else decoy). Classes fixed in `ref/classes.tsv` before mapping.
  - Everything else = decoy.
- 16S: body pair (all counted mates within the 16S +-10 nt), spacer-crossing pair (a counted mate reaches >= 20 nt past
  the 16S 3' end = precursor), edge otherwise. 23S loci have no flank; all counted pairs are body.
- Code check before any new run is read: the new counter must reproduce T5 on `out/t5/SRR11558543.paf`
  (16S div 5 / norm 141 body, div 4 / norm 13 spacer; 23S div 0 / norm 125).

## Verdicts, per run and species group (fixed now)
- r23 = (div 23S pairs / norm 23S pairs) / D, D = expected div:norm pair ratio at equal copy number.
  AS: D = 545/577 (DNA, SRR6051694). CD: D = locus length ratio (div/norm 23S locus length), unless a DNA run is counted.
- Gate G1: norm 23S pairs >= 30, else UNINFORMATIVE (reported anyway).
- **SILENT**: upper end of the two-sided 95% Clopper-Pearson interval of div/(div+norm), converted to r23, < 0.10.
- **ACCUMULATES**: div >= 10 and r23 >= 0.25.
- otherwise **LOW / UNRESOLVED** (r23 and its interval reported).
- Pooled over the new runs of one group as well (sum of pairs), same rule.
- 16S: as T5 (same r16 conversion by length of GCA_009993245.1's copies, precursor guard with T5's k for AS). For CD the
  precursor guard is not calibrated, so CD 16S can only be SILENT or UNRESOLVED.
- Reported with every verdict: mismatch-free (NM=0) counts; spacer pairs per class (div spacer > 0 = the divergent operon
  is transcribed in that run).

## Controls (fixed now)
- C1 rRNA depletion: per run, the share of read pairs whose primary alignment to the rRNA reference is >= 85% identical
  over >= 60 columns (pairs counted on the first 1 M spots). Undepleted total RNA is typically > 80% rRNA. If the share is
  < 30%, the run is called rRNA-depleted. The CG13 pair (selection RT-PCR vs "other") is compared directly.
- C2 depletion probes: probe-based kits (Ribo-Zero type) target conserved rRNA. In silico, for each AS copy (div, norm),
  the fraction of its 23S and 16S covered by blastn hits (word 11) of >= 90% identity over >= 50 nt to E. coli K-12 and
  B. subtilis 168 rRNA. If the divergent copy is covered no more than the normal copy, probes cannot remove it
  preferentially, and depletion cannot explain div = 0 while norm > 0. Probe-free depletion (e.g. Zymo RiboFree / DSN
  normalisation) removes abundant sequences, which would also favour the rarer divergent copy.
- C3 same-species signal: norm 23S and 16S pairs show the species' rRNA survives the prep in that run; divergent-16S
  spacer pairs show the operon is present and transcribed.
- C4 recovery: T5's simulated reads (PC-T5a) recounted with the new counter (div 498/500 at div, norm 500/500).

## What would mean what
- New whole-cell-like runs (CG12/13/22) SILENT for AS: "silent" holds beyond the 2015 series, so it is not a viral-fraction
  artefact of three samples.
- Any run ACCUMULATES for AS: silence is condition-dependent; report the run, its date and r23.
- CD SILENT in AMD RNA: silence holds in a second genus. CD ACCUMULATES: the AS result is species-specific.
- All new runs UNINFORMATIVE: "silent" rests on the 2015 series only; say exactly what was searched.
- A lead (any ACCUMULATES, or a SILENT verdict in a second genus) gets a fresh-subagent recount.

## Amendments
- **A1 (19:00, before any read was mapped).** `mags fetch` warns that the CG runs are "stored aligned"
  (schema alignment_sorted). `vdb-dump --info` shows they were loaded from FASTQ (latf-load), with no REFERENCE table, and
  the first spots carry an N at base 2 (first flowcell tile): the spot order is sequencing order, not genome order, so the
  head is usable as planned.
- **A2 (19:16).** CG12's first count was void: editing src/run_one.sh while bash was executing it re-ran its mapping line and truncated the PAF to 0 lines. Re-mapped from the stored FASTQ with src/map_only.sh (same commands); no extra download.
