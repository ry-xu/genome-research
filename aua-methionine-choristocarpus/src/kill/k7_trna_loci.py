"""Kill check 4b: tRNA complement by blastn. Queries = all annotated tRNAs of 6 standard-code brown algae + Choristocarpus NC.
Targets = Choristocarpus (NC, WGS), Discosporangium contig, and 2 Sphacelariales controls (as a recall check: their own
tRNAs are removed from the query set when they are the target). Per locus: best query, its anticodon mapped through the
alignment, and the best hit from the Ile-CAU class (lysidine tRNA) at that locus.
usage: uv run python -I checks/chor_wgs/kill/k7_trna_loci.py GBDIR DMES_FA"""
import os, re, subprocess, sys
from Bio import SeqIO

GB, DMES = sys.argv[1:3]
tmp = "checks/chor_wgs/kill/tmp"
SP = {"NC_007683.1": "Fucus", "NC_007684.1": "Desmarestia", "NC_085327.1": "Cutleria", "NC_013473.1": "Saccharina",
      "ON552993.1": "Phaeostrophion", "MZ156064.1": "Protohalopteris", "NC_085326.1": "ChorNC"}
ILE_CAU = {("Fucus", 8394), ("Desmarestia", 9300), ("Cutleria", 8055), ("Saccharina", 8353), ("Phaeostrophion", 22257), ("Protohalopteris", 7793)}
MET_F = {("Fucus", 34901), ("Desmarestia", 37536), ("Cutleria", 35861), ("Saccharina", 36010), ("Phaeostrophion", 46360),
         ("Phaeostrophion", 49570), ("Protohalopteris", 40915), ("ChorNC", 2625)}
MET_E = {("Fucus", 3628), ("Desmarestia", 4412), ("Cutleria", 3287), ("Saccharina", 3464), ("Phaeostrophion", 17294), ("Protohalopteris", 2762)}

def acpos(s, want=None):
    c = [m.start() + 2 for m in re.finditer(r"(?=[ACGT]T[ACGT]{3}[AG])", s) if 30 <= m.start() + 2 <= 38]
    if want: c = [x for x in c if s[x:x + 3] == want]
    return c[0] if c else None
def want_of(lab, cls):
    if cls in ("IleCAU", "MetF", "MetE"): return "CAT"
    m = re.search(r"[(-]([ACGUacgu]{3})\)?$", lab)
    return m.group(1).upper().replace("U", "T") if m else None
Q = {}
for acc, sp in SP.items():
    r = SeqIO.read(f"{GB}/{acc}.gb", "genbank")
    for f in r.features:
        if f.type != "tRNA": continue
        st = int(f.location.start) + 1
        s = str(f.extract(r.seq)).upper()
        lab = (f.qualifiers.get("gene") or f.qualifiers.get("product") or ["?"])[0]
        cls = "IleCAU" if (sp, st) in ILE_CAU else "MetF" if (sp, st) in MET_F else "MetE" if (sp, st) in MET_E else lab.replace(" ", "")
        # anticodon offset: annotated if given, else loop pattern
        a = f.qualifiers.get("anticodon", [None])[0]
        off = int(re.search(r"pos:(\d+)", a).group(1)) - st if a and f.location.strand == 1 else acpos(s, want_of(lab, cls))
        Q[f"{sp}|{st}|{cls}"] = (s, off)
qf = os.path.join(tmp, "all_trna.fa")
targets = [("ChorNC", f"{GB}/NC_085326.1.gb"), ("ChorWGS", f"{GB}/CAXNAD020000071.1.gb"), ("Dmes", DMES),
           ("Phaeostrophion", f"{GB}/ON552993.1.gb"), ("Protohalopteris", f"{GB}/MZ156064.1.gb")]
for name, path in targets:
    if path.endswith(".gb"):
        g = str(SeqIO.read(path, "genbank").seq).upper()
    else:
        g = "".join(l.strip().upper() for l in open(path) if not l.startswith(">"))
    tf = os.path.join(tmp, f"tg_{name}.fa"); open(tf, "w").write(f">{name}\n{g}\n")
    with open(qf, "w") as fh:
        for k, (s, off) in Q.items():
            if k.split("|")[0] != name.replace("WGS", "NC") and not (name == "ChorWGS" and k.startswith("ChorNC")):
                fh.write(f">{k}\n{s}\n")
    out = subprocess.run(["blastn", "-task", "blastn", "-word_size", "7", "-evalue", "1e-3", "-dust", "no", "-query", qf, "-subject", tf,
                          "-outfmt", "6 qseqid qstart qend sstart send pident length evalue bitscore qseq sseq"], check=True, capture_output=True, text=True).stdout
    loci = []
    for l in out.splitlines():
        q, qs, qe, ss, se, pid, L, ev, bs, qa, sa = l.split("\t")
        lo, hi = min(int(ss), int(se)), max(int(ss), int(se))
        if int(L) < 50: continue
        hit = (float(bs), q, float(pid), int(L), ev, int(qs), int(ss), int(se), qa, sa)
        for loc in loci:
            if lo < loc[1] - 20 and hi > loc[0] + 20: loc[2].append(hit); loc[0] = min(loc[0], lo); loc[1] = max(loc[1], hi); break
        else: loci.append([lo, hi, [hit]])
    print(f"\n== target {name} ({len(g)} nt): {len(loci)} tRNA-like loci (hits >= 50 nt; queries exclude the target's own species)")
    for lo, hi, hits in sorted(loci):
        hits.sort(reverse=True)
        bs, q, pid, L, ev, qs, ss, se, qa, sa = hits[0]
        off = Q[q][1]; ac = "?"
        if off is not None:   # walk the alignment to the query anticodon (query offset is 0-based; qs 1-based)
            qp, sp_, step = qs - 1, ss, (1 if se > ss else -1); got = []
            for a, b in zip(qa, sa):
                if a != "-" and off <= qp < off + 3: got.append(b)
                if a != "-": qp += 1
            ac = "".join(got)
        ile = [h for h in hits if h[1].endswith("IleCAU")]
        ilestr = f"best IleCAU query: {ile[0][1]} {ile[0][2]}%/{ile[0][3]} e={ile[0][4]}" if ile else "no IleCAU query hit"
        classes = sorted({h[1].split("|")[2] for h in hits})
        print(f"  {lo}-{hi} {'+' if se > ss else '-'} best {q} {pid}%/{L} e={ev} anticodon(mapped)={ac} | classes hit: {','.join(classes)} | {ilestr}")
