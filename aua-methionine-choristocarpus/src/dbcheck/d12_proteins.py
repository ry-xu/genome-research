"""NCBI Protein records (35 WBP from ON552980, 35 YP from NC_085326) and UniParc sequences: which code their sequences
follow, any transl_table / notes / comments; UniProtKB / proteome counts from saved query outputs."""
import re, csv
from Bio import SeqIO
from Bio.Seq import Seq
D = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/data/"
nc = SeqIO.read(D + "NC_085326.gb", "genbank"); on = SeqIO.read(D + "ON552980.gb", "genbank")
def cds_map(rec):
    out = {}
    for f in rec.features:
        if f.type == "CDS":
            s = str(f.extract(rec.seq)).upper()
            out[f.qualifiers["protein_id"][0]] = (f.qualifiers.get("gene", ["?"])[0], s)
    return out
cmap = {**cds_map(nc), **cds_map(on)}
def tr(s, ata_met):
    cod = [s[i:i+3] for i in range(0, len(s) - len(s) % 3, 3)]
    aa = str(Seq(s[:len(cod)*3]).translate(table=1))
    if ata_met:
        aa = "".join("M" if c == "ATA" else a for c, a in zip(cod, aa))
    return aa.rstrip("*")
for fn in ["on_proteins.gp", "nc_proteins.gp"]:
    raw = open(D + fn).read()
    recs = list(SeqIO.parse(D + fn, "genbank"))
    n1 = nM = nother = 0; nAta = 0
    for r in recs:
        pid = r.id; gene, s = cmap[pid]; p = str(r.seq)
        t1, tM = tr(s, False), tr(s, True)
        diff = sum(c == "ATA" for c in [s[i:i+3] for i in range(3, len(s) - 3, 3)])
        nAta += diff
        if p[1:] == t1[1:]: n1 += 1
        elif p[1:] == tM[1:]: nM += 1
        else: nother += 1
    print(f"== {fn}: {len(recs)} protein records; sequence == table-1 translation (ignoring aa 1): {n1}; == AUA->Met translation: {nM}; neither: {nother}; internal ATA codons: {nAta}")
    for term in ["transl_table", "genetic code", "AUA", "isoleucine", "methionine", "/note", "COMMENT"]:
        print(f"   raw-text hits '{term}': {len(re.findall(re.escape(term), raw))}")
    cm = sorted(set(r.annotations.get("comment", "").replace("\n", " ") for r in recs))
    print("   distinct COMMENT texts:", cm)
    coded = sorted(set(re.findall(r'/coded_by="([^"]+)"', raw)))[:3]
    print("   example coded_by:", coded, "| /transl_table in CDS qualifier lines:", len(re.findall(r"/transl_table=", raw)))
# UniParc
rows = list(csv.DictReader(open(D + "uniparc_all.tsv"), delimiter="\t"))
seqs = {}
for pid, (g, s) in cmap.items():
    seqs.setdefault(tr(s, False)[1:], set()).add(pid)
sM = {tr(s, True)[1:] for g, s in cmap.values()}
m1 = sum(r["Sequence"][1:] in seqs for r in rows); mM = sum(r["Sequence"][1:] in sM and r["Sequence"][1:] not in seqs for r in rows)
print(f"\n== UniParc: {len(rows)} entries for the 70 protein ids; UniProtKB accession present in {sum(bool(r['UniProtKB'].strip()) for r in rows)}; "
      f"sequence == table-1 translation: {m1}; == AUA->Met only: {mM}")
print("   organisms field values:", sorted(set(r["Organisms"] for r in rows))[:6])
for fn in ["uniprot_116065.tsv", "uniprot_376270.tsv", "uniprot_proteomes_116065.tsv", "uniprot_proteomes_376270.tsv", "uniprotkb_q2.txt"]:
    n = sum(1 for _ in open(D + fn)) - 1
    print(f"   {fn}: {n} data rows")
for fn in ["on_proteins.gp", "nc_proteins.gp"]:
    notes = {}
    for r in SeqIO.parse(D + fn, "genbank"):
        for f in r.features:
            for n in f.qualifiers.get("note", []):
                notes[n] = notes.get(n, 0) + 1
    print(f"   {fn} distinct /note values: {notes}")
