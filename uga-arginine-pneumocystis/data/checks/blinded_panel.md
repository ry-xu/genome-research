# Blinded panel G01–G10: unblinded 2026-11-11 15:10 EST

Agent: fresh subagent (a0cf3c5f7a276fcce), told only to infer each genome's mitochondrial code from sequence (no
GenBank CDS or protein annotation). Methods: Codetta v2.0, and its own ortholog-column tally against 79 annotated
reference mitogenomes, with a Bonferroni cut of 0.05/640. Its full report is `agent_report.md`. Its scripts and key
outputs are in `agent/` (copied from the scratchpad).

Pass rule (fixed 14:36, before launch; in `key.tsv`):
- TGA called Arg in both TGA leads, by either method, with the other not contradicting;
- AGA called Ser in *Sympodiomycopsis*;
- no TGA or AGA deviation called in any of the 5 controls.

Unscored rows are reported but not scored.

## Verdict: PASS

| ID | Accession | Organism | Role | Prediction | Agent's call | Score |
|---|---|---|---|---|---|---|
| G04 | NC_020331.1 | *Pneumocystis jirovecii* | lead | TGA = Arg | TGA = Arg. Codetta R (N 29); columns R 26/29, p1 2e-13 | hit |
| G07 | PV031684.1 | *Craticula vanensis* | lead | TGA = Arg | TGA = Arg. Codetta R (N 19); columns R 36/38, p1 6e-29 | hit |
| G05 | PP995850.1 | *Sympodiomycopsis kandeliae* | lead | AGA = Ser (also ATA = Met) | AGA = Ser. Codetta S (N 24); columns S 16/21, R 0, p1 2e-17. Also ATA = Met | hit |
| G01 | MH800316.1 | *Proschkinia* sp. | control | table 4 | table 4 | ok |
| G06 | NC_020332.1 | *Pneumocystis murina* | control | table 4 | table 4 | ok |
| G08 | PX580640.1 | Stauroneidaceae sp. | control | table 4 | table 4 | ok |
| G09 | NC_060318.1 | *Pneumocystis wakefieldiae* | control | table 4 | table 4 | ok |
| G10 | PP995849.1 | *Parajaminaea phylloscopi* | control | table 4 (AGA = Arg) | no TGA or AGA change (AGA R 4/4). Extras: ATA = Met, TTA = stop | ok |
| G02 | BK073068.1 | *Golubevia* (TPA) | unscored | AGA too rare; leans Ser | AGA unresolved: 1 column (S); relaxed run S 3 of 6; 0 of 73 AGA at Arg columns | consistent |
| G03 | NC_013660.2 | *Pneumocystis carinii* | unscored | TGA too rare; leans Arg | TGA R at 4 of 4 columns, p1 1.5e-4 (misses the cut) | consistent |

## Extras (not scored)

**ATA = Met in G05 and G10.**
- G05: Codetta M; columns M 22/45, I 4.
- G10: columns M 8/11, I 0; Codetta leans M but is below its threshold.
- This was predicted for G05 (12:14 list) but not for G10.
- AUA = Met is already in my census for 36 genomes, including Pucciniales and *Melanotaenium*, and is a known
  reassignment type (yeast table 3). It is not new.

**TTA = stop (post hoc) in all three Ustilaginomycotina (G02, G05, G10).**
- In-frame use: TTA in frame 2–7 times against 239–270 at the +1/+2 shifts; TTG not used in frame.
- No tRNA with anticodon UAA or CAA.
- G10: all 7 in-frame TTA lie within 5 codons of an alignment end, p 1e-9.
- G05: TTA is the first codon after 4 of 15 aligned gene ends.
- In long-overhang ends, TTA sits at +1 or +2 while the next TAA/TAG is 18–119 codons on:
  - G05: nad6, cox2, orf248;
  - G10: nad4, nad1.
- This agrees with the open README lead "Basidiomycete mito UUR". The agent found it blind to the organism and the
  hypothesis.
- It does not separate "stop" from "unassigned". In these genomes TTA is absent from coding interiors and common in
  non-coding flanks, so any in-frame TTA is expected near a gene end either way.
- My earlier overrun test (20 vs 15.9, P 0.077) and pos0 test (21 vs 10.1, p 0.0018; rank 4 of 64) stand.

**Pneumocystis tRNA.** No genome in the panel has a default-threshold tRNA with anticodon UCA, which matches the
README.

## Caveats

- The agent finalized its calling rule after seeing the data; only the Bonferroni cut was planned. The scored calls
  pass under its relaxed and strict sensitivity settings too:
  - G04 TGA: 31/34 relaxed, 22/23 strict;
  - G07 TGA: 32/32 strict;
  - G05 AGA: S 16/21 under the primary rule.
- *Craticula* TGA = Arg is published (Yılmaz et al. 2026), so G07 is a recovery of known biology, not a new finding.
