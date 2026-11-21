"""Align panel gene models to reference MSAs (MAFFT) and record, for every panel codon, the
reference residues in its alignment column.

usage: align_tally.py refs_dir msa_dir genes_dir out_dir ACC1,ACC2,... [threads]
Outputs out_dir/records.tsv (one row per panel codon) and out_dir/aln/*.
Reference residues in the tally come from the MASKED reference translations (contested codons
CUN/AUA/CGA/CGC and stops -> X), so the tally never relies on how the references' own contested
codons are read. Panel proteins given to MAFFT are also masked the same way ("mask" mode) or
plain table-4 ("t4" mode, sensitivity check).
"""
import os
import sys
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_fasta, write_fasta, translate_codons, codons_of, T4, CONTESTED  # noqa: E402

GENES = ["cox1", "cox2", "cox3", "cob", "atp6", "atp8", "atp9",
         "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"]
MASKSET = CONTESTED | {"TAA", "TAG"}


def mafft(args, out):
    with open(out, "w") as fh, open(out + ".log", "w") as lg:
        subprocess.run(["mafft", "--thread", "1"] + args, check=True, stdout=fh, stderr=lg)


def ref_msa(refs_dir, g, adir):
    src = os.path.join(refs_dir, f"ref_{g}.t4.faa")
    out = os.path.join(adir, f"ref_{g}.msa")
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        # MAFFT does not like '*' -> none present (checked in prep), but guard anyway
        recs = [(n, s.replace("*", "X")) for n, s in read_fasta(src)]
        tmp = os.path.join(adir, f"ref_{g}.in.faa")
        write_fasta(tmp, recs)
        mafft(["--localpair", "--maxiterate", "1000", "--amino", tmp], out)
    return out


def one(job):
    refs_dir, msa_dir, adir, acc, g, mode, models = job
    m = models[g]
    cods = [c[1] for c in m["codons"]]
    prot = translate_codons(cods, T4, mask=MASKSET if mode == "mask" else {"TAA", "TAG"})
    pf = os.path.join(adir, f"{acc}.{g}.{mode}.faa")
    write_fasta(pf, [(f"PANEL|{acc}", prot)])
    msa = os.path.join(msa_dir, f"ref_{g}.msa")
    out = os.path.join(adir, f"{acc}.{g}.{mode}.aln")
    if not (os.path.exists(out) and os.path.getsize(out) > 0):  # reuse an existing alignment
        mafft(["--add", pf, "--localpair", "--maxiterate", "1000", "--amino", msa], out)
    aln = dict(read_fasta(out))
    # masked/t4 reference sequences
    rmask = {n: s for n, s in read_fasta(os.path.join(refs_dir, f"ref_{g}.mask.faa"))}
    rt4 = {n: s for n, s in read_fasta(os.path.join(refs_dir, f"ref_{g}.t4.faa"))}
    refs = [n for n in aln if not n.startswith("PANEL")]
    L = len(aln[f"PANEL|{acc}"])
    # project masked residues onto aligned rows; record per-ref column of each residue
    col_res = [dict() for _ in range(L)]
    ref_col_of = {}
    for n in refs:
        row = aln[n]
        k = 0
        cols = []
        for j, ch in enumerate(row):
            if ch != "-":
                assert ch.upper() == rt4[n][k].replace("*", "X").upper(), (n, j, ch, rt4[n][k])
                col_res[j][n] = rmask[n][k]
                cols.append(j)
                k += 1
        ref_col_of[n.split("|")[0]] = cols  # gene models store the bare reference accession
    prow = aln[f"PANEL|{acc}"]
    rows = []
    k = 0
    for j, ch in enumerate(prow):
        if ch == "-":
            continue
        gp, cod, href, hq = m["codons"][k]
        # HSP-direct column (cross-check): column of the anchoring ref residue hq
        hcol = ref_col_of[href][hq - 1] if hq is not None and href in ref_col_of and hq - 1 < len(ref_col_of[href]) else -1
        res = "".join(col_res[j].get(n, "-") for n in refs)
        hres = "".join(col_res[hcol].get(n, "-") for n in refs) if hcol >= 0 else ""
        rows.append((acc, g, mode, k + 1, gp, cod, T4[cod], j, res, hcol, hres))
        k += 1
    assert k == len(m["codons"])
    return rows, refs


def main():
    refs_dir, msa_dir, genes_dir, out_dir, accs = sys.argv[1:6]
    threads = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    accs = accs.split(",")
    adir = os.path.join(out_dir, "aln")
    os.makedirs(adir, exist_ok=True)
    present = [g for g in GENES if read_fasta(os.path.join(refs_dir, f"ref_{g}.t4.faa"))]
    jobs = []
    for acc in accs:
        models = json.load(open(os.path.join(genes_dir, f"{acc}.genes.json")))
        for g in present:
            if g in models:
                for mode in ("mask", "t4"):
                    jobs.append((refs_dir, msa_dir, adir, acc, g, mode, models))
    with ThreadPoolExecutor(threads) as ex:
        results = list(ex.map(one, jobs))
    with open(os.path.join(out_dir, "records.tsv"), "w") as fh:
        fh.write("acc\tgene\tmode\tidx\tgpos0\tcodon\tt4aa\tcol\tref_res\thsp_col\thsp_ref_res\n")
        for rows, refs in results:
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
    # reference order per gene (same for all jobs of a gene)
    order = {}
    for (rows, refs), job in zip(results, jobs):
        order[job[4]] = refs
    json.dump(order, open(os.path.join(out_dir, "ref_order.json"), "w"), indent=1)
    print("jobs", len(jobs), "records", sum(len(r) for r, _ in results))


if __name__ == "__main__":
    main()
