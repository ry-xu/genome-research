"""Which genes separate the TGA-reading Throscidae symbionts from the non-readers? Exploratory: a mechanism hint, not a test.
Proteins come from the tree run (results/throscidae_tree/<label>.faa, pyrodigal table 4). A protein of genome X counts as
present in genome Y if tblastn (genetic code 4, e <= 1e-5) finds it in Y's genome sequence, so missed gene calls and split
pseudogenes still count as present. Reported: genes found in every genome of one group and in none of the other.
Each is annotated by its best Swiss-Prot hit (DIAMOND).
usage: uv run --no-project --with pyrodigal,requests python -I src/gene_content.py   (after src/throscidae_tree.py)
writes results/throscidae_tree/gene_content.tsv"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE.parent / "stop-codon-census" / "src"))
import readthrough as rt  # noqa: E402

rt.GEN = HERE / "data" / "sym_genomes"
TREE = HERE / "results" / "throscidae_tree"
SPROT = HERE.parent / "stop-codon-census" / "data" / "sprot.dmnd"
GENOMES = {"Tser_sericeus": "GCA_040210545.2", "Tder_dermestoides": "GCA_040211365.2", "Ttur_turgidus": "GCA_040409785.2",
           "Tcar_carinifrons": "GCA_040211755.2", "Tduv_duvalii": "GCA_040210955.2", "Tmey": "GCA_037094695.1",
           "Aulonothroscus_AUS03": "GCA_040212085.2", "S_bostrichidophilus": "GCA_029851105.1"}
READERS = ["Tser_sericeus", "Tder_dermestoides", "Ttur_turgidus"]
NONREADERS = ["Tcar_carinifrons", "Tduv_duvalii", "Tmey", "Aulonothroscus_AUS03"]


def main():
    work = TREE / "gene_content"
    work.mkdir(exist_ok=True)
    for lab, acc in GENOMES.items():
        fna = work / f"{lab}.fna"
        if not fna.exists():
            fna.write_text("".join(f">{n}\n{s}\n" for n, s in rt.read_fasta(rt.fetch(acc))))
            subprocess.run(["makeblastdb", "-in", str(fna), "-dbtype", "nucl", "-out", str(work / lab)], check=True,
                           capture_output=True)
    present = {}  # (query genome, protein) -> set of genomes with a tblastn hit
    for q in READERS + NONREADERS:
        for lab in GENOMES:
            out = work / f"{q}__{lab}.tbn"
            if not out.exists():
                subprocess.run(["tblastn", "-query", str(TREE / f"{q}.faa"), "-db", str(work / lab), "-db_gencode", "4",
                                "-evalue", "1e-5", "-outfmt", "6 qseqid", "-max_target_seqs", "1", "-out", str(out)],
                               check=True)
            for line in open(out):
                present.setdefault((q, line.strip()), set()).add(lab)
    rows = []
    for (q, p), labs in present.items():
        if q in NONREADERS and set(NONREADERS) <= labs and not labs & set(READERS):
            rows.append((q, p, "non-readers only", labs))
        if q in READERS and set(READERS) <= labs and not labs & set(NONREADERS):
            rows.append((q, p, "readers only", labs))
    seqs = {}
    for q in READERS + NONREADERS:
        name = None
        for line in open(TREE / f"{q}.faa"):
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                seqs[(q, name)] = line.strip()
    qf = work / "candidates.faa"
    qf.write_text("".join(f">{q}|{p}\n{seqs[(q, p)]}\n" for q, p, _, _ in rows))
    ann = {}
    if rows:
        hits = work / "candidates_sprot.m8"
        subprocess.run(["diamond", "blastp", "-q", str(qf), "-d", str(SPROT), "-o", str(hits), "--quiet", "-k", "1",
                        "-e", "1e-5", "--sensitive", "--outfmt", "6", "qseqid", "stitle", "pident", "evalue"], check=True)
        for line in open(hits):
            qid, title, pid, ev = line.rstrip("\n").split("\t")
            ann[qid] = f"{title} ({pid} %, {ev})"
    with open(TREE / "gene_content.tsv", "w") as out:
        out.write("genome\tprotein\tpattern\tlen_aa\tin_S_bostrichidophilus\tswissprot_best\n")
        for q, p, pat, labs in rows:
            out.write(f"{q}\t{p}\t{pat}\t{len(seqs[(q, p)])}\t{'S_bostrichidophilus' in labs}\t{ann.get(f'{q}|{p}', '-')}\n")
    print(open(TREE / "gene_content.tsv").read())


if __name__ == "__main__":
    main()
