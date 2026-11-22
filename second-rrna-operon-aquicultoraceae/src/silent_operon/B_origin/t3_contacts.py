"""Classify every 23S and 16S nucleotide of E. coli 70S (PDB 7K00: mRNA, A- and P-site tRNA) by what it touches
(heavy atoms <= 4 A): tRNA/mRNA ('ligand'), r-protein ('protein'), the other subunit's rRNA ('bridge'), or none.
Writes out/contacts_7k00.tsv: kind, ecoli_pos, ligand, protein, bridge, proteins (names)
"""
import gzip, gemmi

st = gemmi.read_structure("lab/work/frozen-core/results/interface/pdb/7K00.cif.gz")
st.setup_entities()
model = st[0]
ent_of = {}
for ent in st.entities:
    for sub in ent.subchains:
        ent_of[sub] = ent.name if hasattr(ent, "name") else ""
# map chains to roles by entity description
desc = {}
doc = gemmi.cif.read("lab/work/frozen-core/results/interface/pdb/7K00.cif.gz")
b = doc.sole_block()
eid = list(b.find_values("_entity.id")); edesc = list(b.find_values("_entity.pdbx_description"))
emap = dict(zip(eid, [d.strip("'") for d in edesc]))
role = {}
for ch in model:
    sub = ch.get_polymer().subchain_id() if len(ch.get_polymer()) else None
    if sub is None:
        continue
    ent = st.get_entity_of(ch.get_polymer())
    d = emap.get(ent.name, "") if ent else ""
    role[ch.name] = d
print({k: v for k, v in role.items() if "rRNA" in v or "tRNA" in v or "mRNA" in v})

ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
out = open("out/contacts_7k00.tsv", "w")
out.write("kind\tecoli_pos\tligand\tprotein\tbridge\tproteins\n")
for ch in model:
    r = role.get(ch.name, "")
    if r not in ("16S rRNA", "23S rRNA"):
        continue
    kind = r.split()[0]
    other = "23S rRNA" if kind == "16S" else "16S rRNA"
    for res in ch.get_polymer():
        lig = prot = bri = 0; pn = set()
        for atom in res:
            if atom.element.name == "H":
                continue
            for m in ns.find_atoms(atom.pos, "\0", radius=4.0):
                cra = m.to_cra(model)
                if cra.atom.element.name == "H":
                    continue
                rr = role.get(cra.chain.name, "")
                if "tRNA" in rr or rr == "mRNA":
                    lig = 1
                elif "ribosomal protein" in rr:
                    prot = 1; pn.add(rr.split()[-1])
                elif rr == other:
                    bri = 1
        out.write(f"{kind}\t{res.seqid.num}\t{lig}\t{prot}\t{bri}\t{','.join(sorted(pn))}\n")
out.close()
