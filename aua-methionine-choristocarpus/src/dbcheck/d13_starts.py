"""Start and stop codons of annotated CDS in NC_085326 and the Genoscope contig, and first residue of stored translations."""
from collections import Counter
from Bio import SeqIO
D = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/data/"
for acc in ["NC_085326", "CAXNAD020000071.1"]:
    rec = SeqIO.read(D + acc + ".gb", "genbank")
    st = Counter(); sp = Counter(); first = Counter(); ata_start = []
    for f in rec.features:
        if f.type != "CDS": continue
        s = str(f.extract(rec.seq)).upper(); st[s[:3]] += 1; sp[s[-3:]] += 1
        first[f.qualifiers.get("translation", ["?"])[0][0]] += 1
        if s[:3] == "ATA": ata_start.append(f.qualifiers.get("gene", f.qualifiers.get("locus_tag", ["?"]))[0])
    print(f"{acc}: start codons {dict(st)}; stop codons {dict(sp)}; first aa of /translation {dict(first)}; ATA-start CDS: {ata_start}")
