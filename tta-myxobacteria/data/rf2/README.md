# Release factor 2 (RF2 / `prfB`) sequences from all 17 TTA-stop Polyangiaceae branch genomes

Full-length (`ATG` start to `TGA` stop) RF2 sequences reconstructed directly from the 17 TTA-stop branch genomes (10 GTDB + 7 outside-GTDB catalogue genomes, including the cultured predator–prey strain **`GCA_982529355.1` / 1-FT3.2**, listed first).

## Files
- `rf2_proteins.faa` — Full-length, frameshift-resolved protein sequences (380–384 aa)
- `rf2_cds_frameshift_fixed.fna` — Full-length nucleotide CDS with the single skipped `T` of the `CTTTGAC` +1 frameshift removed (ready for direct expression / synthesis without relying on programmed frameshifting)
- `rf2_cds_genomic_raw.fna` — Raw genomic DNA (`ATG` to `TGA`, retaining the native `CTTTGAC` +1 frameshift site)
- `rf2_summary.tsv` — Coordinates, key active-site/pocket residues (*E. coli* RF2 numbering), frameshift positions, and metadata

## Key notes for expression / biochemical assays
- **Programmed +1 frameshift (`CTTTGAC`):** All 17 genes carry the canonical bacterial `prfB` +1 frameshift motif `CTT-TGA-C` at residue 23 (or residue 26 in `bin496`, `spire_mag_00879468`, and `GCA_060838245.1`), where the ribosome slips +1 from `CTT` (Leu) into `GAC` (Asp). Use `rf2_cds_frameshift_fixed.fna` (or `rf2_proteins.faa`) for standard recombinant expression.
- **Zero internal `TTA` codons:** None of the 17 RF2 genes contains any internal `TTA` codon; all terminate cleanly at `TGA`.
- **Motif & key residues (*E. coli* numbering):**
  - All 17 carry **D143G**, **L201F**, **A329M**, and intact **GGQ**.
  - **15 of 17** carry **SPF → SPY** (including cultured strain `GCA_982529355.1` / 1-FT3.2).
  - **2 of 17** (`GCA_046351315.1` / `JAXLYQ01`, and `bin474`) retain **SPF** alongside **D143G + L201F + A329M**.
