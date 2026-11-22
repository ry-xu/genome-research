# Second fact-check of the rewritten outbox draft (fresh subagent, finished ~10:40 14 Nov). Saved by Lyra, condensed.
Outputs: scratchpad/factcheck2_chor/ (dl/, code/ab.py, out/rna_strand.out).
- Counts and figure: OK (b_read.out, kill/k1_audit, factcheck_chor/out/recount.out agree).
- WRONG: "Genoscope and GenBank give the same counts" (88/6 vs 88/5).
- WRONG attribution + number: stranded-library result came from fact-check 1, not the kill agent; "116 of 126" does not
  reproduce: rerun of fc1 rna.py on fc1/rna/r1_head.fastq.gz (1.18M read-1s) gives 126 sense / 10 antisense
  (out/rna_strand.out). Earlier run probably on a partial file. "RNA, not DNA" -> "mostly RNA".
- OVERSTATED: "Not chance substitutions" (AUG->AUA is a single transition; AUU/AUC transversions) -> cut.
- ONE-SIDED: composition bullet. Under AUA = Met, Met = 35.5/1,000 (AUG+AUA), 1.5x the highest other brown alga (22.9),
  mostly in ribosomal-protein genes (Altair groups.out 4.16x). Must be stated in the post.
- OVERSTATED: screen recall: 7,728 records declared AUA = Met, only 58 testable (>=10 AUA at M/I), 1 flagged.
- Strain: ON552980 has no strain qualifier (place/date match KU-1152). Mito differences: 27 of 32,351 (25 run-length indels
  + 2 adjacent substitutions TA<->AT, rps19). Use "most likely one strain". LC664917-47 (Akita 2022) has Genoscope
  co-authors, so the KU-2346 genes may come from the same DNA: they show the label, not a separate sequencing.
- Minor: "one per order" -> "from 8 different orders"; NCBI files both genera under Sphacelariales; AB899175 covers 4 of 88
  sites -> "not errors of either assembly"; "Ehara 1997 used the same kind of evidence" unchecked (full text not read).
- Length ~1,140 words: cut the substitution test, the duplicate composition line, fold strict rule and strain into one line.
- Could not check: tblastn rerun (relied on 3 agreeing implementations), proteome, tRNA charging, which culture Genoscope
  sequenced, Ehara 1997 full text, Denoeud 2024 Table S2, full RNA run.
