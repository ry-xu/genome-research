# Aminoacyl-tRNA synthetases (aaRS) for TAG-reassigned lineages L2, L3, L5, L6

Extracted from the raw genome assemblies and verified by Prodigal (complete `partial=00` calls with canonical `ATG` start codons and `AGGAG`/`GGAG` Shine-Dalgarno motifs) and full Swiss-Prot DIAMOND searches.

## Files
- `aars_proteins.faa` — Full-length protein sequences (8 records across L2, L3, L5, L6)
- `aars_cds.fna` — Full-length nucleotide CDS (`ATG` to `TAA`/`TGA`)
- `aars_summary.tsv` — Coordinates, lengths, internal TAG counts, and Swiss-Prot matches
- Corresponding CUA tRNAs: [`../trnas/cua_trnas.fa`](../trnas/cua_trnas.fa)

## Key notes for expression / biochemical assays
- **Zero internal `TAG` codons:** None of these eight aaRS genes contains an internal `TAG` codon (all terminate cleanly at `TAA` or `TGA`), so the sequences can be synthesized/expressed directly under the standard genetic code without codon recoding inside the aaRS itself.
- **L2 (GluRS, `gltX`):** Present in `GCA_045599485.1` (486 aa); the other L2 assembly (`GCA_024404275.1`, 74.7% complete) is missing the `glnS`–`gltX` operon.
- **L3 (GluRS, `gltX`):** Present in both independent assemblies (`GCA_937907055.1` and `GCA_045628485.1`, 480 aa, differing at 4 amino acids).
- **L5 (GlnRS, `glnS`):** Present in both L5 species (`GCA_056232285.1`, 553 aa; `GCA_056261325.1`, 556 aa). Neither genome encodes GatCAB, confirming the direct GlnRS charging pathway.
- **L6 (TrpRS, `trpS`):** Present in `GCA_056209055.1` (341 aa) and independently confirmed in two CNGB MAGs of the same species (`QTP_TA20190011_bin.237`, 100% identical protein; `QTP_TA20190013_bin.296`, 1 aa difference).
- **3'-CCA note on the paired CUA tRNAs (`../trnas/cua_trnas.fa`):**
  - `L2_SFMI01_CUA` (76 nt) and `L3_ZC3RG03_CUA` (76 nt) **already include** a genomically encoded 3'-`CCA`.
  - `L5_Lachno_CUA` (73 / 74 nt) and `L6_HGM11604_CUA` (72 nt) **end at discriminator base N73** (no genomic 3'-`CCA`; append `CCA` for in vitro transcription if needed).
