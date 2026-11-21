# Arm T audit (fresh subagent, raw genomes; 20 Nov). Written by Saiph from the agent's reply (it could not write files).
Scripts a1-a8 and .out files in this folder hold the numbers.
- Re-extraction: own pyaragorn run = same 1,366 calls; anticodon_offset 1-based in all; T-loop-anchored acceptor parse
  agrees on all 1,177 parsable tRNAs.
1. D 0.262/0.262/0.214 CONFIRMED; 0.16-0.30 under other masking/gap choices (local alignment 0.09). Longer sister
   tRNAs lower D. Pro is the least similar class for 2 targets, 18th of 21 for the third. The UAG tRNA sits at the same
   locus (Thr-Ser-[UAG]-Arg-Trp) in all targets and sisters: synteny supports orthology. OVERSTATED as three results:
   two target sequences are identical, the third differs at 2 nt -> one event in one clade.
2. "Every tRNA-Pro has C1:G72" WRONG: 48 of 49 genomes; the exception is W. siamensis BK074538's own tRNA-Pro (A1:U72).
   C1:G72 occurs in no other class (0 of 1,083 non-Pro, non-UAG tRNAs). W. versatilis gained it without reassignment;
   S. pinicola lacks it. Not evidence of Pro charging on its own. N73 = A everywhere.
3. S. pinicola D = 0.105 CONFIRMED but weak: 2nd-lowest of 23 calibration values. One sister input is a bad 75-nt call
   at another locus; without it D = 0.123. G1:U72 in no sister at the shared locus.
4. No missed xAG tRNA CONFIRMED. A probable 2nd Pro(UGG) locus at 273-345 in S. pinicola (reads CCN).
5. Controls CONFIRMED. Sc control fragile: nearest class His under the frozen metric, Trp if gap columns excluded.
Bugs (none changes a verdict): t0b.py:25 overlap test misses <30-nt hits; t0b.py:30 drops loci with <5 hits (RULES: any
hit); t0b.py:28 rerun uses default not relaxed thresholds; origin.py:14 includes the bad S. pinicola sister call;
origin.py:99 calibration averages paralogs; origin.py:124 Ashbya calibration printed vs Pro while D uses Ala.
Saiph follow-up: Ashbya calibration recomputed vs Ala (TGC): 20 values, 0.122-0.419; Ashbya D 0.094 still below all.
No script produced t0b/q.fa and the BLAST hits (done inline); the auditor's a3_t0b.py reproduces the sweep.
