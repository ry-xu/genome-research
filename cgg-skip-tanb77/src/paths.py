"""Where the large inputs live. None of them ship with this folder (see scripts/setup_reference_data.sh); each can be
moved with an environment variable. The original run used the lab's shared genome cache (CODON_ATLAS) laid out as
  genomes/<acc[:8]>/<acc>.fna.gz            one gzipped FASTA per NCBI assembly (codon_atlas/cache.py fetches them)
  products/<acc[:8]>/<acc>/trna.parquet     ARAGORN tRNA calls (codon_atlas/trna.py; columns listed there)
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
ATLAS = Path(os.environ.get("CODON_ATLAS", DATA / "atlas"))
GENOMES = ATLAS / "genomes"
PRODUCTS = ATLAS / "products"
# GTDB R232 metadata (bac120/ar53_metadata_r232.tsv.gz)
GTDB_META = Path(os.environ.get("GTDB_META", DATA / "gtdb_meta"))
# UniProtKB/Swiss-Prot FASTA and NCBI taxdump rankedlineage.dmp, for src/prok_db.py
SPROT_FASTA = Path(os.environ.get("SPROT_FASTA", DATA / "uniprot_sprot.fasta.gz"))
RANKEDLINEAGE = Path(os.environ.get("RANKEDLINEAGE", DATA / "taxdump" / "rankedlineage.dmp"))
# ARAGORN screen of the GTDB R232 species representatives (per-genome anticodon counts), from an earlier
# lab project; not redistributed. Rebuild with codon_atlas/trna.py over the representatives.
TRNA_SCREEN = Path(os.environ.get("TRNA_SCREEN", DATA / "trna_screen" / "genomes.tsv"))
