# Note on cat_members.tsv (2026-11-11 16:55 EST, Basil)

GCA_004144535.1 and GCA_004144565.1 are separate rows here. They are probably two bins of one genome (same four
DNA-SIP samples, PRJNA272922; only 7 shared clade-like markers where ~48 would be expected; see
`scratch/factcheck3e/REPORT.md` and `notes/candidates.md`), so they are counted once. Frozen Arm B on the merged bins
is in `armB_cat4m.tsv`: TTA rank 1, 6 families, 5.9x, log10 P -3.19, below the bar.

GCA_060838125.1 is the same strain as GCA_060838225.1 (ANI 99.95%) and is also not counted separately.

So the table's 17 rows with TTA at rank 1 make 15 genomes: the 7 Arm B passes in `cat_candidates.tsv` plus 8 below the bar.
