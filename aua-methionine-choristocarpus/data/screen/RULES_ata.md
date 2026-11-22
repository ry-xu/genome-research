# ATA read as Met where the declared table says Ile: test rules
Written 08:43 14 Nov by Lyra, after seeing: the GenBank K1 remainder (10 rows), the two unfollowed RefSeq remainder rows
(Pycnococcus NC_013935.1, Choristocarpus NC_085326.1), checks/k1_gb/b_codons.py (5 genomes) and the first survey
printout (c_survey.py; genomes with ATA M-share >= 0.5 listed by eye). So the lineages are already known to me; these
rules fix how they are counted and controlled, not which ones come out.

Observations: ledger cols, locus pid < 95 (K1-1), (acc, sid, spos, codon, res) counted once (K1-2).
Extra filters (new here):
- start: drop observations at anchor position spos <= 30 (alternative start codons align to anchor Met near the N end).
- echo: keep only anchors whose own genome's declared table reads ATA as Ile, or prokaryote anchors. Anchors from
  tables reading ATA as Met (2, 3, 5, 13, 21, ...; computed from Bio.Data.CodonTable) are dropped.
Population: genomes whose declared table (read_k1.declared, per data dir) reads ATA as Ile.
Flag (genome): ATA_M >= 10 and ATA_M >= 2 * ATA_I, and the control codon ATT behaves: ATT_M <= 0.05 * (ATT_M + ATT_I).
Controls, same filters and thresholds:
- false-flag rate: the same flag computed for ATT (ATT_M >= 10, >= 2*ATT_I) and for ATC, over the same population.
- recall: genomes declared with an ATA=Met table (expected flag) and the flag rate among them.
Counting: flagged genomes grouped by order and class; identical-count records (same organism, same counts) noted as
duplicates. Organelle kind reported apart (mito / plastid).
Then: per lineage, a fresh kill subagent (raw records; conserved-column test as in checks/k1_remainder/kill; tRNA
genes, e.g. trnI(CAU) / tilS presence) and a prior-art subagent, before any thread post calls anything new.
