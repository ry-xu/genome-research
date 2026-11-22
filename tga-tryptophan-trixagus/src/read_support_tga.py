"""Raw-read support for in-frame TGA sites (TGA read as Trp): stop-codon-census read_support.py unchanged (streamed
SRA reads, minimap2 -x sr, reads with > 3% mismatches ignored, per-site codon tally, TGG sites as the control),
except that genes are called with table 4 and the sites are in-frame TGA instead of TAG.
A real TGA = Trp genome shows TGA fixed in reads at its in-frame TGA sites; assembly errors or a mixed population don't.
usage: uv run --no-project --with pyrodigal,requests python -I src/read_support_tga.py <SRR> <acc> <total_spots> [streams=4]
-> results/read_support/<acc>_<SRR>.txt and .sites.tsv"""
import sys
from pathlib import Path

import pyrodigal

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE.parent / "stop-codon-census" / "src"))
import read_support as rs  # noqa: E402
import readthrough as rt  # noqa: E402

rt.GEN = HERE / "data" / "sym_genomes"
rs.OUT = HERE / "results" / "read_support"


def sites_t4(contigs, codon="TGA"):
    finder = pyrodigal.GeneFinder()
    finder.train(*[s.encode() for s in contigs.values()], translation_table=4)
    sites = {}
    for c, s in contigs.items():
        for g in finder.find_genes(s.encode()):
            nt = g.sequence()
            for i in range(0, len(nt) - 3, 3):  # the last codon (the stop) is excluded
                if nt[i:i + 3] == codon:
                    sites[(c, g.begin - 1 + i if g.strand == 1 else g.end - 3 - i)] = g.strand
    return sites


_summarise = rs.summarise
rs.tag_sites = sites_t4
rs.summarise = lambda sites, obs, codon: _summarise(sites, obs, "TGA" if codon == "TAG" else codon)

if __name__ == "__main__":
    rs.main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 4)
