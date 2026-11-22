"""Two follow-ups on the Throscidae symbionts after src/gene_content.py. Exploratory, no pass rule.
1. tRNA-Trp: tRNAscan-SE -B on each genome; aligns the tRNA-Trp genes (MAFFT) and lists alignment columns where the 3
   readers share one base that no non-reader has. Also prints the anticodon loop and the D-stem 11:24 pair.
2. In-frame TGA inside table-4 genes (pyrodigal): for each site, is it flanked by >= 10 aa of the gene's best Swiss-Prot
   alignment on both sides (DIAMOND), or near the gene's end? A sense TGA sits inside aligned sequence; a stop TGA sits
   at or past the alignment end, with the table-4 call running on to the next TAA/TAG.
usage: uv run --no-project --with pyrodigal,requests python -I src/throscid_trp.py   (after src/gene_content.py)
writes results/throscidae_tree/trna_trp.txt, tga_sites.tsv"""
import subprocess
from collections import Counter
from pathlib import Path

import pyrodigal

HERE = Path(__file__).resolve().parent.parent
TREE = HERE / "results" / "throscidae_tree"
GC = TREE / "gene_content"  # .fna files written by gene_content.py
SPROT = HERE.parent / "stop-codon-census" / "data" / "sprot.dmnd"
TRNASCAN = Path.home() / "dev" / "tools" / "bin" / "tRNAscan-SE"
READERS = ["Tser_sericeus", "Tder_dermestoides", "Ttur_turgidus"]
NONREADERS = ["Tcar_carinifrons", "Tduv_duvalii", "Tmey", "Aulonothroscus_AUS03"]
OUTGROUP = ["S_bostrichidophilus"]
ALL = READERS + NONREADERS + OUTGROUP


def read_fasta(p):
    d, n = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = []
        elif n:
            d[n].append(line)
    return {k: "".join(v) for k, v in d.items()}


def trna_trp():
    work = TREE / "trna"
    work.mkdir(exist_ok=True)
    trp, lines = {}, []
    for lab in ALL:
        ss = work / f"{lab}.ss"
        if not ss.exists():
            subprocess.run([str(TRNASCAN), "-B", "-q", "--brief", "-f", str(ss), "-o", str(work / f"{lab}.out"),
                            str(GC / f"{lab}.fna")], check=True, capture_output=True)
        seq = st = None
        typ = False
        for line in open(ss):
            if line.startswith("Type:"):
                typ = line.split()[1] == "Trp"
            elif typ and line.startswith("Seq:"):
                seq = line.split()[1]
            elif typ and line.startswith("Str:"):
                st = line.split()[1]
                typ = False
        trp[lab] = seq.upper()
        runs, i = [], 0  # stems in order: acceptor 5', D 5', D 3', anticodon 5', anticodon 3', ...
        while i < len(st):
            if st[i] in "<>":
                j = i
                while j < len(st) and st[j] == st[i]:
                    j += 1
                runs.append((i, j))
                i = j
            else:
                i += 1
        d5, d3 = seq[runs[1][0]:runs[1][1]].upper(), seq[runs[2][0]:runs[2][1]].upper()
        ac_loop = seq[runs[3][1]:runs[4][0]].upper()
        pair = f"{d5[1]}:{d3[-2]}" if len(d5) == 4 and len(d3) == 4 else "D-stem not 4 bp"
        lines.append(f"{lab}\tanticodon loop {ac_loop}\tD-stem 11:24 {pair}\t{seq}")
    fa = work / "trp.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in trp.items()))
    aln = _mafft(fa)
    lines.append("")
    for i in range(len(aln[READERS[0]])):
        r = {aln[k][i].upper() for k in READERS}
        n = {aln[k][i].upper() for k in NONREADERS}
        if len(r) == 1 and not r & n:
            lines.append(f"reader-only column {i + 1}: readers {r}, non-readers {n}, outgroup {aln[OUTGROUP[0]][i].upper()}")
    lines += [""] + [f"{k.ljust(22)}{v.upper()}" for k, v in aln.items()]
    (TREE / "trna_trp.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def _mafft(fa):
    out = fa.with_suffix(".aln")
    with open(out, "w") as fh:
        subprocess.run(["mafft", "--quiet", "--globalpair", str(fa)], check=True, stdout=fh)
    return read_fasta(out)


def tga_sites():
    work = TREE / "tga_sites"
    work.mkdir(exist_ok=True)
    recs, faa = [], []
    for lab in ALL:
        g = read_fasta(GC / f"{lab}.fna")
        finder = pyrodigal.GeneFinder(closed=True)
        finder.train(*[s.encode() for s in g.values()], translation_table=4)
        k = 0
        for s in g.values():
            for gene in finder.find_genes(s.encode()):
                k += 1
                nt = gene.sequence().upper()
                prot = gene.translate(translation_table=4).rstrip("*")
                gid = f"{lab}|{k}"
                faa.append(f">{gid}\n{prot}\n")
                recs += [(lab, gid, i // 3 + 1, len(prot)) for i in range(0, len(nt) - 3, 3) if nt[i:i + 3] == "TGA"]
    (work / "genes.faa").write_text("".join(faa))
    subprocess.run(["diamond", "blastp", "-q", str(work / "genes.faa"), "-d", str(SPROT), "-o", str(work / "genes.m8"),
                    "--quiet", "-k", "1", "-e", "1e-5", "--sensitive", "--outfmt", "6", "qseqid", "qstart", "qend",
                    "stitle"], check=True)
    hit = {}
    for line in open(work / "genes.m8"):
        q, qs, qe, t = line.rstrip("\n").split("\t")
        hit.setdefault(q, (int(qs), int(qe), t))
    summ = {lab: Counter() for lab in ALL}
    with open(TREE / "tga_sites.tsv", "w") as out:
        out.write("genome\tgene\ttga_codon\tprotein_len\tcodons_to_stop\tsprot_aln\tclass\tsprot_best\n")
        for lab, gid, pos, plen in recs:
            h = hit.get(gid)
            if h and h[0] + 10 <= pos <= h[1] - 10:
                cls = "inside alignment"
            elif h and pos > h[1]:
                cls = "past alignment end"
            elif h:
                cls = "alignment edge or before start"
            else:
                cls = "gene without Swiss-Prot hit"
            summ[lab][cls] += 1
            summ[lab]["all"] += 1
            out.write(f"{lab}\t{gid}\t{pos}\t{plen}\t{plen + 1 - pos}\t{f'{h[0]}-{h[1]}' if h else '-'}\t{cls}\t"
                      f"{h[2] if h else '-'}\n")
    cls = ["all", "inside alignment", "past alignment end", "alignment edge or before start", "gene without Swiss-Prot hit"]
    print("genome\t" + "\t".join(cls))
    for lab in ALL:
        print(lab + "\t" + "\t".join(str(summ[lab][c]) for c in cls))


if __name__ == "__main__":
    trna_trp()
    tga_sites()
