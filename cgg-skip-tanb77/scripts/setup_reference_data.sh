#!/usr/bin/env bash
# Download the public reference data the analysis needs into data/ (all gitignored). Run from anywhere.
#   - GTDB R232 metadata                         -> data/gtdb_meta/{bac120,ar53}_metadata.tsv.gz
#   - UniProtKB/Swiss-Prot + NCBI taxdump         -> data/uniprot_sprot.fasta.gz, data/taxdump/rankedlineage.dmp
#   - Swiss-Prot clustered at 50% identity, prokaryote and virus clusters only, as a DIAMOND database
#                                                 -> data/sprot_c50.{faa,tsv}, data/sprot_p50.{faa,dmnd}
#     The original run used the Swiss-Prot release current on 2026-11-11; a newer release can shift individual hits.
#     The clustering command was not logged; `diamond cluster --approx-id 50` writes the same two-column
#     (representative, member) table that src/prok_db.py reads.
#   - Genomes of the scanned sets (NCBI)          -> data/atlas/genomes/, with ARAGORN tRNA calls in data/atlas/products/
#
# Not fetched here (large, or needed only to redo one part):
#   - The 24 TANB77 genomes of MGnify pig-gut v1.0 (RNA test): for each MGYG id in data/rna/pig_genomes.tsv,
#       https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/pig-gut/v1.0/species_catalogue/<id[:11]>/<id>/genome/<id>.{fna,gff,faa}
#     into data/proteomics/mgnify_pig/ (the default of src/rna_cgg.py; RNA_MDIR moves it).
#   - RNA and DNA reads: streamed from NCBI SRA by src/rna_stream.py (run lists in data/rna/*_runs.txt); the Illumina
#     run of the discovery genome is ERR11860625.
#   - The chicken-gut TANB77 genomes of the CUU replication (data/cuu_context/chicken_genomes.tsv): NCBI assemblies,
#     with pyrodigal gene calls, in data/rna_chicken/genomes/<acc>.{fna,gff}.
#   - UHGV v1.0 phages (Zenodo 10.5281/zenodo.17402089): votus_full.fna.zst and votus_metadata_extended.tsv.zst
#     into data/uhgv/. Pfam-A.hmm (InterPro FTP) into data/pfam/ for the Pfam junction check.
#   - The GTDB-wide ARAGORN tRNA screen read by src/analysis/cgg_reader_census.py (TRNA_SCREEN in src/paths.py) came
#     from an earlier lab project and is not redistributed; src/trna_products.py makes the same calls per genome.
# Needs: curl, diamond (2.2.8 was used). For the other parts: minimap2, sra-tools (fastq-dump), miniprot, BLAST+,
# tRNAscan-SE 2, HMMER, zstd.
set -euo pipefail
cd "$(dirname "$0")/.."
D=data
GTDB=https://data.gtdb.ecogenomic.org/releases/release232/232.0
mkdir -p $D/gtdb_meta $D/taxdump
for f in bac120_metadata ar53_metadata; do
  [ -s $D/gtdb_meta/$f.tsv.gz ] || curl -fL -o $D/gtdb_meta/$f.tsv.gz $GTDB/${f}_r232.tsv.gz
done
[ -s $D/uniprot_sprot.fasta.gz ] || curl -fL -o $D/uniprot_sprot.fasta.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
if [ ! -s $D/taxdump/rankedlineage.dmp ]; then
  curl -fL -o $D/taxdump/new_taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
  tar -xzf $D/taxdump/new_taxdump.tar.gz -C $D/taxdump rankedlineage.dmp
fi
if [ ! -s $D/sprot_p50.dmnd ]; then
  [ -s $D/sprot_c50.tsv ] || diamond cluster -d $D/uniprot_sprot.fasta.gz -o $D/sprot_c50.tsv --approx-id 50 --quiet
  uv run python -I src/prok_db.py                       # -> data/sprot_p50.faa
  diamond makedb --in $D/sprot_p50.faa -d $D/sprot_p50 --quiet
fi
for s in tanb77_set clos_lowgc_set cag698_set moll_set tanb77_reader_set; do
  uv run python -I src/trna_products.py $D/genome_sets/$s.tsv > /dev/null   # fetches each genome, then tRNA calls
done
echo "reference data ready in $D/"
