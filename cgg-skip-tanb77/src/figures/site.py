"""Extract the rpoC skip site in GCA_963931445.1 (OZ007016.1, minus strand) from the genome FASTA,
translate it with the C skipped, and align to B. subtilis and C. difficile RpoC translated from their genomes.
Prints the window used for the image. Reads only; writes site.json next to this script."""
import gzip, json, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

ATLAS = os.environ.get("CODON_ATLAS", "data/atlas") + "/genomes"   # config: the genome cache (src/paths.py)
GEN = f"{ATLAS}/GCA_9639/GCA_963931445.1.fna.gz"   # the blinded re-check's G4 (recheck/01_refs.py FOCAL)
REFS = {"Bacillus subtilis": ("GCF_000009045.1", "NC_000964.3", 125562, 129161, 1),   # from G4.paf
        "Clostridioides difficile": ("GCF_000009205.2", "NC_009089.1", 93918, 97403, 1)}
CTX = "TCAGGAGCAAGACGGAAATGTTAA"   # site_ctx in results/families/GCA_963931445.1.fs.tsv.gz

g = {r.id: str(r.seq).upper() for r in SeqIO.parse(gzip.open(GEN, "rt"), "fasta")}
chrom = g["OZ007016.1"]
lo, hi = 1548100, 1551800
rc = str(Seq(chrom[lo:hi]).reverse_complement())
i = rc.find(CTX); assert i >= 0 and rc.find(CTX, i + 1) < 0
c = i + 12                       # index of the C of CGG in rc
assert rc[c:c + 3] == "CGG"
plus0 = hi - 1 - c               # plus-strand 0-based coordinate of that C (as G on plus)
# frame: gene region on minus strand; find in-frame start upstream by walking back codons from c
up = rc[c % 3:c]                 # same frame as CGG
dn = rc[c + 1:]                  # shifted frame after skipping C
pu = str(Seq(up[len(up) % 3 and 0:len(up) - len(up) % 3]).translate())
pdn = str(Seq(dn[:len(dn) - len(dn) % 3]).translate())
# trim upstream to after last stop, downstream to first stop
pu = pu[pu.rfind("*") + 1:]
pdn = pdn[:pdn.find("*")] if "*" in pdn else pdn
skip_prot = pu + pdn
gpos = len(pu)                   # index of the glycine (first residue after skip)
noskip = str(Seq(rc[c - 3 * 0:c + 3 * 12]).translate())   # old frame from CGG onward

al = PairwiseAligner(); al.mode = "local"
al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -11, -1
out = {"genome": "GCA_963931445.1", "contig": "OZ007016.1", "strand": "-", "c_plus0": int(plus0),
       "dna_up": rc[c - 30:c], "c": rc[c], "dna_dn": rc[c + 1:c + 1 + 30],
       "prot_up": pu[-10:], "prot_dn": pdn[:10], "old_frame_from_cgg": noskip, "refs": {}}
print("skip protein len", len(skip_prot), "G index", gpos, skip_prot[gpos - 10:gpos + 10])
print("old frame from CGG:", noskip)
for name, (acc, cid, b, e, s) in REFS.items():
    sub = ATLAS + f"/{acc[:8]}/{acc}.fna.gz"
    with gzip.open(sub, "rt") as f:
        seq = {r.id: str(r.seq).upper() for r in SeqIO.parse(f, "fasta")}[cid]
    nt = seq[b - 1:e] if s == 1 else str(Seq(seq[b - 1:e]).reverse_complement())
    rp = str(Seq(nt).translate(table=11)).rstrip("*")
    a = al.align(skip_prot, rp)[0]
    # map query positions to subject
    qmap = {}
    for (qs, qe), (ss, se) in zip(*a.aligned):
        for k in range(qe - qs): qmap[qs + k] = ss + k
    win = range(gpos - 10, gpos + 10)
    row = "".join(rp[qmap[k]] if k in qmap else "-" for k in win)
    ident = sum(1 for k in win if k in qmap and rp[qmap[k]] == skip_prot[k])
    print(name, rp[:5], "len", len(rp), "score", a.score, "G maps to residue", qmap.get(gpos, -1) + 1,
          "codon", nt[3 * qmap[gpos]:3 * qmap[gpos] + 3] if gpos in qmap else None)
    print("  TANB77 ", skip_prot[gpos - 10:gpos + 10]); print("  ref    ", row, f"{ident}/20 identical")
    out["refs"][name] = {"acc": acc, "window": row, "res": int(qmap.get(gpos, -1)) + 1,
                         "codon": nt[3 * qmap[gpos]:3 * qmap[gpos] + 3], "ident20": int(ident)}
out["prot_window"] = skip_prot[gpos - 10:gpos + 10]
json.dump(out, open(os.path.join(os.path.dirname(__file__), "site.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "refs"}, indent=1))
