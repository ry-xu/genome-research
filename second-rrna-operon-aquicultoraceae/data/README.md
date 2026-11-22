# Data

Small tables only. The genomes (GTDB release 232), the read runs (SRA) and the large alignment and mapping files are
public or too large and are not included. Reports in `checks/` and `silent_operon/` were copied from the lab record;
paths in them that are not in this folder refer to the original working folder.

## `ptc/` (catalytic-centre bases, the first genome)
- `columns.tsv`: census 23S column table (E. coli position, consensus base, conservation `f_cons`, pairing partner).
- `fcons_ptc.out`: conservation of the changed PTC positions.
- `ptc_coords.tsv`, `ptc_pairwise.out`: the PTC bases of each copy, by structure-free pairwise alignment.
- `aq_genomes.tsv`, `aq_loci_class.tsv`, `aq_pairing.tsv`: the Aquicultorales genomes, their 23S loci (normal or
  divergent) and the base-pairing scores of each locus.
- `pile.out`, `dxG.pileup.txt`: read coverage of the DASXVE01 operon and contigs (reads at ≥ 97% identity).

## `kill_23s/` (the kill agent's own checks; code in `../src/checks/kill_23s/`)
- `prevalence.out`: 23 Aquicultoraceae genomes, 18 with a 23S locus, 15 with a divergent 23S.
- `pid_div.out`, `pid_norm.out`: identity between genera, divergent and normal copies.
- `q_vs_lsuall.tsv`: census hits of each divergent and normal 23S. `loci_blast.tsv`: identity of each family
  locus to the divergent and the normal copy.
- `pairing_free*.out`: model-free base pairing of the divergent copies and of functional controls.
- `synteny.out`, `contig_affinity.out`, `asd.out`: flanking genes, contig taxonomy and the anti-Shine-Dalgarno check.

## `scan/` (multi-copy 23S scan of 50,594 genomes)
- `class_full.tsv`, `check_full.tsv`, `chim_full.tsv`, `post_full.tsv`, `flag_tech.tsv`: flagged copies, their
  checks, chimera test and assembly technology. `disposition_outfam.tsv`: what happened to each out-of-family flag.

## `expression_23s/`, `expression_16s/` (RNA and DNA read counts, first three RNA runs)
- `*.loci.paf`, `loci.paf`: reads of the DNA run SRR6051694 and the RNA runs SRR11558543, SRR11558542 and SRR11558425
  (`loci.paf`) mapped to the loci; `loci_class_all.tsv` labels each locus. `src/frozen_core/fig_div23.py` counts them.
- `SRR11558543.counts.txt`, `class16.tsv`, `sim.paf`: the 16S body and spacer counts, the 16S locus classes and the
  simulated-read check.

## `contacts_16s/`, `tail_16s/` (16S tests)
- `t6*.txt`, `pc_t6a.txt`, `contacts_16s23s.tsv`: 16S–23S contact and covariation tests (PDB 6OG7 contacts).
- `carriers.jsonl`, `t1_*.log`, `pc2_*.log`: the tail-pairing (private message channel) test and its controls.

## `uba6777/`
- `ptc_ux.out`, `coab.out`: PTC bases and co-abundance read test of UBA6777 sp002731275.

## `checks/`
Kill, audit, fact-check, rules and prior-art reports of 14–19 Nov (see `../notes.md`, "Checks").

## `silent_operon/` (the checks of 22 Nov)
Each check's rules (`RULES_*.md`, written before the results) and result with its fresh-agent recheck (`RESULT_*.md`).
- `processing/` (`RESULT_processing.md`): `operons.tsv` (mature ends of each operon, full-CM cmsearch),
  `t1_stems.tsv`, `t1b_16S.tsv`, `t1b_23S.tsv`, `t1c_helix.tsv` (RNase III leader-trailer stems), `t2_upstream.tsv`,
  `t2_genes.tsv`, `t2_cons.tsv`, `t2_spacer_boxa.tsv` (promoter, boxA, upstream gene and conservation),
  `t3_content.tsv`, `5s.tbl` (5S and spacer tRNA content).
- `origin/` (`RESULT_origin.md`, `RESULT_neighbourhood.md`, `RESULT_selection.md`): `loci.tsv`, `units.tsv` (loci
  used), `t1_23S.contree`, `t1_16S.contree` (gene trees, UFBoot), `t1_au_*.txt` (AU tests), `t2_orthogroups.tsv`,
  `t2_unit_genes.tsv` (genes within 10 kb of each operon), `smallorf.faa`, `smallorf.aln` (the upstream peptide),
  `t3_*_pairs.tsv`, `t3_*_sites.tsv`, `t3_stems*.tsv` (pairing, site-class rates, processing stems),
  `contacts_7k00.tsv` (tRNA/mRNA/protein contacts from PDB 7K00), `aln23.sto` (the 23S cmalign alignment; also the
  input of `../src/fig_div_ptc_rows.py`).
- `rna/` (`RESULT_rna.md`): `<run>_head<N>.counts.txt`, the counts for the first N read pairs of each new RNA run
  (SRR11558439 = CG05, SRR6241672 = CG13 undepleted, SRR6241673 = CG13 depleted, SRR6261448, SRR6261398, SRR7589086,
  SRR33092280-82); `recheck_SRR6241672_counts.tsv`, the fresh agent's recount of CG13; `norm.blast`, `div.blast`, the
  probe-similarity check.
