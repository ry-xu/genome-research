"""T6 contact sets: 16S-23S intersubunit contacts in PDB 6OG7 (E. coli 70S, non-rotated; file from work/rf2-interface).
usage: uv run --no-project --with gemmi python -I src/t6_contacts.py OUT.tsv
Heavy atoms, 4.0 A, author (E. coli) numbering. Rows: 16S residue, 23S residue, min distance."""
import sys
import gemmi

PDB = 'data/external/pdb/6og7.cif.gz'
CUT = 4.0
st = gemmi.read_structure(PDB); st.setup_entities(); st.remove_hydrogens()
doc = gemmi.cif.read(PDB); blk = doc.sole_block()
desc = dict(zip(blk.find_loop('_entity.id'), blk.find_loop('_entity.pdbx_description')))
model = st[0]; ent = {s: e for e in st.entities for s in e.subchains}
names = {}
for ch in model:
    p = ch.get_polymer()
    if len(p): names[ch.name] = gemmi.cif.as_string(desc.get(ent[p.subchain_id()].name, '')) if p.subchain_id() in ent else ''
s16 = [c for c, d in names.items() if '16S' in d]; s23 = [c for c, d in names.items() if '23S' in d]
print('16S chains', s16, '23S chains', s23, file=sys.stderr)
assert len(s16) == 1 and len(s23) == 1
ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
best = {}
for r in model[s16[0]].get_polymer():
    for a in r:
        for m in ns.find_atoms(a.pos, '\0', radius=CUT):
            cra = m.to_cra(model)
            if cra.chain.name != s23[0]: continue
            d = a.pos.dist(cra.atom.pos); k = (r.seqid.num, r.name, cra.residue.seqid.num, cra.residue.name)
            best[k] = min(best.get(k, 9), d)
with open(sys.argv[1], 'w') as o:
    o.write('n16\tb16\tn23\tb23\tdist\n')
    for k in sorted(best): o.write('\t'.join(map(str, k)) + f'\t{best[k]:.2f}\n')
print(len(best), 'pairs;', len({k[0] for k in best}), '16S residues;', len({k[2] for k in best}), '23S residues', file=sys.stderr)
# numbered E. coli sequences straight from the model (author numbering, insertion codes kept)
for c, tag in ((s16[0], '16S'), (s23[0], '23S')):
    with open(sys.argv[1].replace('contacts_16s23s.tsv', f'ecoli{tag}_num.tsv'), 'w') as o:
        for r in model[c].get_polymer():
            t = gemmi.find_tabulated_residue(r.name)
            o.write(f'{r.seqid.num}{r.seqid.icode.strip()}\t{t.one_letter_code.upper() if t else "N"}\t{r.name}\n')
