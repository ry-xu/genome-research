"""T2: neighbourhood of the divergent vs normal operons.
Orthogroups (single linkage, e<=1e-5, >=50% cover of the shorter protein), annotation (Swiss-Prot best hit, Pfam,
phage proteins), conservation per orthogroup across divergent and normal loci of Aquicultoraceae, suggestive classes.
Usage: python -I src/t2_flanks.py  -> out/t2_orthogroups.tsv, out/t2_unit_genes.tsv; prints the summary
"""
import csv, re
from collections import defaultdict

units = {r["unit"]: r for r in csv.DictReader(open("out/units.tsv"), delimiter="\t")}
fg = list(csv.DictReader(open("out/flank_genes.tsv"), delimiter="\t"))
alen = {r["gene_id"]: int(r["aa_len"]) for r in fg}

# orthogroups
par = {g: g for g in alen}
def find(x):
    while par[x] != x:
        par[x] = par[par[x]]; x = par[x]
    return x
for line in open("out/flank_all.m8"):
    q, s, pid, L, ql, sl, ev, bs = line.split("\t")
    if q == s:
        continue
    if int(L) >= 0.5 * min(int(ql), int(sl)):
        a, b = find(q), find(s)
        if a != b:
            par[a] = b
og = {g: find(g) for g in alen}

sp = {}
for line in open("out/flank_sprot.m8"):
    p = line.rstrip("\n").split("\t")
    if p[0] in sp:
        continue
    t = p[8]
    name = re.sub(r" OS=.*", "", t.split(" ", 1)[1]) if " " in t else t
    sp[p[0]] = f"{name} ({float(p[2]):.0f}%)"
ph = {}
for line in open("out/flank_phage.m8"):
    p = line.rstrip("\n").split("\t")
    if p[0] not in ph and float(p[2]) >= 30:
        ph[p[0]] = f"{p[8][:60]} ({float(p[2]):.0f}%)"
pf = defaultdict(list)
for line in open("out/flank_pfam.domtbl"):
    if line.startswith("#"):
        continue
    p = line.split()
    if float(p[12]) <= 1e-5 and p[3] not in pf[p[0]]:
        pf[p[0]].append(p[3])

CLASSES = {
    "ribosomal protein": r"ribosomal protein|Ribosomal_",
    "rRNA/tRNA modification": r"methyltransferase|Methyltransf|MTase|rRNA .*(pseudouridine|synthase)|PseudoU|RsmB|RsmC|RsmD|RsmE|RsmG|RlmN|SpoU|TrmH|Rlm|Rsm|KsgA|Erm",
    "RNase / nuclease": r"[Rr]ibonuclease|RNase|nuclease|Ribonuc|RNB|Nuclease|endonuc|exonuc|YbeY",
    "toxin-antitoxin": r"[Tt]oxin|[Aa]ntitoxin|PIN|RelE|ParE|MazF|PemK|HicA|HicB|HigB|HigA|VapB|VapC|Doc|Fic|YoeB|YafQ|AbiE|GNAT|ParD|Phd|MqsR|HipA|DUF1863|RelB|DinJ|CcdB|Txe|YefM",
    "defence": r"CRISPR|Cas[0-9]|Cas_|restriction|Restriction|Methylase_|Abi|retron|Gabija|DndB|Dnd|BREX|DISARM|Zorya|Septu|Thoeris|Lamassu|Wadjet|Druantia",
    "mobile": r"[Tt]ransposase|[Ii]ntegrase|[Rr]ecombinase|Phage_integrase|rve|DDE_Tnp|Transposase|IS[0-9]|resolvase|Resolvase|phage|Phage|prophage|Prophage|terminase|Terminase|capsid|portal|tail",
    "translation factor / ribosome-associated": r"[Ee]longation factor|[Ii]nitiation factor|GTP-binding|GTPase|Era|EngA|ObgE|HflX|RsfS|RaiA|hibernation|RMF|YfiA|HPF|Ribosome-binding|ribosome-associated|LepA|Der|RbfA|RimM|RimP",
}


def classify(g):
    txt = " ".join([sp.get(g, ""), " ".join(pf.get(g, [])), ph.get(g, "")])
    return [c for c, rx in CLASSES.items() if re.search(rx, txt)]


AQ = {u: r for u, r in units.items() if r["family"] == "Aquicultoraceae"}
gen_of = lambda u: AQ[u]["genus"]
pres = defaultdict(lambda: {"div": defaultdict(set), "norm": defaultdict(set), "ctl": set(), "order": set()})
for r in fg:
    u = r["unit"]; U = units[u]; o = og[r["gene_id"]]
    if r["side"] == "in":
        continue
    if U["family"] == "Aquicultoraceae":
        pres[o][U["label"]][U["genus"]].add(u)
    elif U["family"] == "control":
        pres[o]["ctl"].add(U["genus"])
    else:
        pres[o]["order"].add(U["genus"])

div_units = [u for u in AQ if AQ[u]["label"] == "div"]
norm_units = [u for u in AQ if AQ[u]["label"] == "norm"]
div_genera = sorted({AQ[u]["genus"] for u in div_units})
print("divergent units:", len(div_units), "in genera", div_genera, "| normal units:", len(norm_units))

members = defaultdict(list)
for g in alen:
    members[og[g]].append(g)


def og_name(o):
    names = [sp[g].split(" (")[0] for g in members[o] if g in sp]
    pfs = sorted({d for g in members[o] for d in pf.get(g, [])})
    best = max(set(names), key=names.count) if names else "-"
    return best, ",".join(pfs[:4])


rows = []
for o, d in pres.items():
    ndg = len(d["div"]); nng = len(d["norm"])
    ndu = sum(len(v) for v in d["div"].values()); nnu = sum(len(v) for v in d["norm"].values())
    cls = sorted({c for g in members[o] for c in classify(g)})
    nm, pfs = og_name(o)
    rows.append((ndg, ndu, nng, nnu, o, nm, pfs, cls, sorted(d["div"]), sorted(d["norm"]), sorted(d["order"]), sorted(d["ctl"])))
rows.sort(key=lambda r: (-r[0], -r[1]))
with open("out/t2_orthogroups.tsv", "w") as fh:
    fh.write("div_genera\tdiv_units\tnorm_genera\tnorm_units\tog\tname\tpfam\tclasses\tdiv_in\tnorm_in\torder_norm_in\tcontrols\n")
    for r in rows:
        fh.write("\t".join(str(x) if not isinstance(x, list) else ",".join(x) for x in r) + "\n")
print("\nOrthogroups next to the divergent operon in >= 3 genera (div genera / div units / norm genera / norm units):")
for r in rows:
    if r[0] >= 3:
        print(f"  {r[0]} / {r[1]:2d} / {r[2]} / {r[3]:2d}  {r[5][:55]:55s} {r[6][:40]:40s} {','.join(r[7])}  order:{len(r[10])} ctl:{','.join(r[11])}")
print("\nOrthogroups next to normal operons in >= 3 genera:")
for r in sorted(rows, key=lambda r: (-r[2], -r[3])):
    if r[2] >= 3:
        print(f"  {r[0]} / {r[1]:2d} / {r[2]} / {r[3]:2d}  {r[5][:55]:55s} {r[6][:40]:40s} {','.join(r[7])}  order:{len(r[10])}")
print("\nSuggestive-class genes within 10 kb, by unit type (count of genes; units with any):")
for lab, us in [("div", div_units), ("norm", norm_units), ("ctl", [u for u in units if units[u]["family"] == "control"])]:
    cnt = defaultdict(int); uu = defaultdict(set)
    for r in fg:
        if r["unit"] in us and r["side"] != "in":
            for c in classify(r["gene_id"]):
                cnt[c] += 1; uu[c].add(r["unit"])
    print(f"  {lab:5s} ({len(us)} units): " + "; ".join(f"{c} {cnt[c]} genes/{len(uu[c])} units" for c in CLASSES if cnt[c]))
with open("out/t2_unit_genes.tsv", "w") as fh:
    fh.write("unit\tgenus\tlabel\tgene\tside\tdist\trel_strand\taa\tog\tsprot\tpfam\tphage\tclasses\n")
    for r in fg:
        U = units[r["unit"]]
        g = r["gene_id"]
        fh.write(f"{r['unit']}\t{U['genus']}\t{U['label']}\t{g}\t{r['side']}\t{r['dist']}\t{r['rel_strand']}\t{r['aa_len']}\t{og[g]}\t"
                 f"{sp.get(g,'-')}\t{','.join(pf.get(g,[])) or '-'}\t{ph.get(g,'-')}\t{','.join(classify(g)) or '-'}\n")
