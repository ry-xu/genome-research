# Audit of the CUU/CUA = Pro finding (fresh subagent, 2026-11-11 ~12:55)

Saved by Clover from the subagent's final message (the harness blocked it from writing this file). Scripts in this
folder are the subagent's; `census_set_reads.out` and `iso_census.out` were rerun from here at 13:0x.

**Verdict.** No code bug creates or inflates the CUU/CUA-at-Pro signal. An independent check with no profile HMMs and
no census consensus (`indep.py`: table-4 six-frame ORFs, BLOSUM62 alignment to relatives' proteins) reproduces it in
all four genomes and not in the same-genus control. Errors were in bookkeeping and wording.

## Code
1. **Low.** The read step's consensus came from gbmitoorf only; the census, `pro_leu` and the contrast use mito +
   gbmito + gbmitoorf. 3,635/3,635 read-table rows match the gbmitoorf-only consensus. Site set 77/81/81/43 vs the
   census's 78/82/82/44 (cob:342 and cob:404 missing, Vir_cox1:462 extra). Also affects `reads_isolates/` (81 vs 82).
   Fix: refiltered to the census site set (`census_set_reads.py`, `iso_census.py`); docstring of `reads_codons.py` now
   says to pass all census files.
2. **Negligible.** Deletions count as non-matching; an insertion inside a codon is missed. Can only lower agreement.
3. **Notes.** MAPQ 0 comes from the duplicated ORF reference and cannot inflate agreement. Same-run reads test
   assembly fidelity only. OR703884 is the same run (identical after a 36,482 bp rotation, 0 mismatches). Possible
   double-counting of mates leaves percentages unchanged. Fisher synonyms include the other CTN codons (conservative).
4. **Checked and correct.** census_org frame/offsets, strand, HMM-walk indels, CIGAR walk incl. reverse strand,
   0x904 filter, strict (≥0.9) definition. The four queries have weight 1/105 each; dropping them changes 0 strict
   columns; lowest top share at CTN-at-Pro columns 0.903. `pro_leu`, `cun_CTA.tsv`, `cun_CTT.tsv` reproduce exactly.

## Independent reproduction (`indep.py`, relatives BK074529 and BK074524, 10 genes)

| query | at relatives' Pro (n) | CTN | CCN | at relatives' Leu UUA/UUG (n) | UUA/UUG | CTN |
|---|---|---|---|---|---|---|
| BK074538 | 114 | 108 | 0 | 462 | 395 | 0 |
| BK074532 | 114 | 106 | 2 | 462 | 395 | 1 |
| BK074542 | 114 | 109 | 0 | 462 | 405 | 0 |
| BK074262 | 114 | 49 | 48 | 464 | 345 | 0 |
| control BK074539 | 113 | 0 | 103 | 459 | 124 | 237 |

With *Sugiyamaella* relatives (BK074251, BK074259): *S. pinicola* 57 CTN and 55 CCN at Pro, 0 CTN at Leu.

## Reads on the census site set

| genome | sites | reads | matching | lowest site (≥3 reads) | MAPQ≥20 |
|---|---|---|---|---|---|
| BK074532 | 78/78 | 3,782 | 98.86% | 94.4% | 35 sites, 98.60%, lowest 87.5% |
| BK074538 | 82/82 | 46,781 | 99.60% | 96.2% | 40 sites, 99.60% |
| BK074542 | 82/82 | 13,970 | 99.57% | 96.9% | 2 sites, 98.95% |
| BK074262 | 44/44 | 16,746 | 98.92% | 94.3% | 20 sites, 99.32%, lowest 70% (nad5:245, 7/10) |
| ERR11532760 → BK074538 | 82/82 | 1,450 | 98.76% | 88.9% (8/9) | – |

*S. pinicola* CCN sites: lowest 90.3% (nad5:12, 28/31). Isolate: CTA 1,234/1,248, CTT 198/202.

## Wrong or loose sentences (all fixed in README and post)
- Read coverage numbers used the wrong site set (above).
- "Every site with ≥3 reads matches in ≥94%" holds only for CTN-at-Pro sites at MAPQ 0.
- Blinded controls were three same-genus genomes plus *Starmerella aceti* BK074222, not four same-genus.
- The blinded subagent was told table 3 was a candidate code (not "told nothing about the codon").
- "TGA = Trp and ATA = Ile rules out table 3": TGA is Trp in table 3 too; ATA = Ile and no CUN = Thr rule it out.
- Wolters 2023 used Codetta v2.0 on 353 mitogenomes and mention only "the CUN reassignment"; they listed OR703884
  (this *W. infanticola*) as table 4.
- Table 3 is all four CUN = Thr and is used beyond Saccharomycetaceae; *Ashbya* is itself in Saccharomycetaceae.
- tRNA-Pro similarity 0.67 is *W. siamensis* only.
- PYCC 7069 and the type strain CBS 13331 were both isolated from sugar-cane leaves in Thailand: possibly the same
  isolate, separate culture and sequencing.

Unverifiable: "no report of any other codon read as Pro" (the auditor's searches found none either); "no such editing
known in fungal mitochondria" (author's knowledge).

Sources checked: Ling et al. 2014 NAR 42:499; Wolters et al. 2023 Front Microbiol 14:1268944 (Table S1); Christinaki
et al. 2022 GBE evac073; Sarhan et al. Sci Data 13:28, PRJNA1367877; OR703884.1 (PRJNA998421); ENA ERR11532760,
SAMEA113572364, PRJEB62807; Khunnamwong et al. 2014 IJSEM (type strain DMKU-SE106T = CBS 13331).
