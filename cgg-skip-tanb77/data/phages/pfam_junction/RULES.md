# The Pfam-check agent cross-check of the phage agent 18:27 (TANB77 phages read CGG by the host's +1 skip). Frozen at 18:43 15 Nov, before any output.

Question: are the skip-signature CGGs in T1/Tf phages read through into the +1 frame as one protein, and does that
include core phage (virion/packaging) proteins? This tests the host-DNA-contamination alternative and re-tests
junction continuity with a method other than kill1's DIAMOND (Pfam HMM domains).

Inputs (read only): work/code-walls/data/phages.fna; kill1/out/genes_meta.tsv (pyrodigal meta calls);
kill1/out/meta.cgg_sites.tsv and meta.ctrl_sites.tsv; kill1/out/my_groups.tsv (set, quality). k is checked against
the sequence (codon at k must equal the table's codon) before use; mismatches are counted and dropped.

Sets (all GC bins, every quality):
- SKIP: rows with set in {T1, Tf}, codon CGG, skip = 1.
- CTRL: in the same phages, gene ends with skip = 1 rows on other codons and no CGG skip row in that gene. One per
  gene, fused at a codon whose distance from the gene end (n-1-k) is drawn from the SKIP distance distribution
  (seed 1); if that codon is not a ctrl row, the nearest ctrl row by distance is used.
- RCGG: R rows with codon CGG and skip = 1 (the rare R signature sites).
Fused protein: gene nt up to codon k (excluded), then nt from codon k's 2nd base onward (the +1 frame), to the first stop
(genome sequence beyond the gene end, strand-aware). Junction residue j = k (0-based).

Search: hmmsearch (1 thread) of every fused protein against Pfam-A (work/23rd-amino-acid/data/pfam_dl/Pfam-A.hmm),
--domE 1e-5. SPANNED = some domain with i-Evalue <= 1e-5 whose alignment (hmm-aligned seq coords) covers j-10..j+10.
HALLMARK = a spanned domain whose Pfam name is also in pfam_dl/phage.hmm and whose name or description matches
capsid|portal|terminase|tail|baseplate|head|prohead|scaffold (case-insensitive).

Report: SPANNED / sites for SKIP, CTRL, RCGG; HALLMARK count in SKIP, with genome quality; SPANNED share in SKIP for
Complete genomes. Reading: SKIP spanned share well above CTRL (and spanned hallmark proteins in Complete genomes)
supports one protein read across the skip in genuine phage genes; SKIP ~ CTRL means the open +1 frames are mostly
overlapping genes and the continuity claim weakens. No threshold claim beyond this description; numbers are reported as is.

Rules for me: nothing outward-facing; no email or personal details to outside services; load-gated 1 thread.

Amendment 1 (18:43, after build.py counts only, before any search): sets are 1,543 SKIP, 8,704 CTRL, 4,129 RCGG fused
proteins. To fit one thread, hmmsearch runs on a random 600 per set (python random.Random(1).sample over sites.tsv
order, per set). Shares are reported over the 600.
