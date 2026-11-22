"""GenBank/RefSeq records: header info, CDS transl_table/notes, tRNA features with anticodons, any AUA/code text."""
import re
from Bio import SeqIO
from Bio.Seq import Seq
D = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/data/"
for acc in ["ON552980", "NC_085326", "CAXNAD020000071.1", "CAXNAE020000396.1"]:
    raw = open(D + acc + ".gb").read()
    rec = SeqIO.read(D + acc + ".gb", "genbank")
    print(f"\n===== {acc}: {rec.id} len={len(rec)} {rec.description}")
    print("  organism:", rec.annotations.get("organism"), "| date:", rec.annotations.get("date"))
    print("  taxonomy:", "; ".join(rec.annotations.get("taxonomy", [])))
    for ref in rec.annotations.get("references", []):
        print("  REF:", ref.title, "|", ref.journal)
    cm = rec.annotations.get("comment", "")
    print("  COMMENT:", cm.replace("\n", " ")[:800] if cm else "(none)")
    src = [f for f in rec.features if f.type == "source"][0]
    print("  source quals:", {k: v for k, v in src.qualifiers.items() if k not in ("db_xref",)}, src.qualifiers.get("db_xref"))
    cds = [f for f in rec.features if f.type == "CDS"]
    tt = {}
    for f in cds:
        k = f.qualifiers.get("transl_table", ["(absent)"])[0]; tt[k] = tt.get(k, 0) + 1
    print(f"  CDS count={len(cds)} transl_table values: {tt}")
    notes = [(f.qualifiers.get('gene', ['?'])[0], n) for f in rec.features for n in f.qualifiers.get("note", [])]
    print(f"  /note qualifiers ({len(notes)}):")
    for g, n in notes:
        print(f"     {g}: {n}")
    for term in ["AUA", "ATA", "genetic code", "transl_table", "isoleucine", "methionine", "code"]:
        print(f"  raw-text hits '{term}': {len(re.findall(term, raw, flags=re.I if term not in ('AUA','ATA') else 0))}")
    # translation consistency: does stored /translation match table-1 translation? count ATA codons inside CDS
    nata = 0; nX = 0; mism = 0
    for f in cds:
        s = f.extract(rec.seq)
        cod = [str(s[i:i+3]) for i in range(0, len(s) - len(s) % 3, 3)]
        nata += sum(c == "ATA" for c in cod)
        if "translation" in f.qualifiers:
            t1 = str(Seq(str(s)).translate(table=1)).rstrip("*")
            st = f.qualifiers["translation"][0]
            if t1[1:] != st[1:]:
                mism += 1
    print(f"  ATA codons in annotated CDS: {nata}; CDS whose stored /translation != table-1 translation (ignoring first aa): {mism}")
    trna = [f for f in rec.features if f.type == "tRNA"]
    print(f"  tRNA features: {len(trna)}")
    for f in trna:
        q = f.qualifiers
        print(f"     {str(f.location):28s} product={q.get('product',['?'])[0]:22s} gene={q.get('gene',['-'])[0]:10s} "
              f"anticodon={q.get('anticodon',['-'])[0]} note={q.get('note',['-'])[0]}")
