"""One row per stage 2 report: the codon's tests next to its control's.

usage: stage2_table.py LIST > table.tsv     (LIST: lines "UNIT CODON MEANING"; reads results/stage2/UNIT_CODON.txt)
flag (triage only, the verdict is read per unit): split = GC3 p < 0.001 with a GC3 gap > 0.05, or codon-usage distance
> 2x the control's, or dispersion > 2x the control's; weak = share of the meaning at conserved columns < half the
control's own share; ok otherwise. org = share of the codon's contigs whose best hit is a mitochondrion- or
plastid-encoded protein (by Swiss-Prot mnemonic), where TGA=W (code 4) or other organelle codes are expected.
"""
import re
import sys

HOME = "lab/work/stop-sense-euk"
ORG = re.compile(r"\|(COX[123]|CYB|ATP[689]|NU[1-6L]+M|NU4LM|AI\d*M|MATK|PSB.|PSA.|RBL|ATP[ABEFHI]|RPOC?\d?|YCF\d+|"
                 r"RR\d+|RPS\d+|RPL\d+|CCMF?.?|TATC|NDH.|PET.)_")
COLS = ["n", "bodies", "disp", "gc3_with", "gc3_without", "gc3_p", "use_dist", "pid_with", "pid_without", "cons_n",
        "cons_share"]


def organelle(path):
    k = n = 0
    on = False
    for line in open(path):
        if line.startswith("contigs with"):
            on = True
        elif on and line.startswith("  "):
            n += 1
            k += bool(ORG.search(line.split("\t")[1]))
    return f"{k}/{n}"


def parse(path):
    rows, own = {}, None
    head = None
    for line in open(path):
        x = line.rstrip("\n").split("\t")
        if line.startswith("meaning's own codons"):
            own = float(line.split("ident=")[1])
        elif x[0] == "codon":
            head = x
        elif head and len(x) == len(head) and x[0] != "codon":
            rows[x[0]] = dict(zip(head, x))
        elif line.startswith("contigs with"):
            break
    return rows, own


def num(v):
    try:
        return float(v)
    except ValueError:
        return float("nan")


def main(lst):
    print("unit\tcodon\tmeaning\tcontrol\t" + "\t".join(f"{c}\tctl_{c}" for c in COLS) + "\tflag\torg")
    for line in open(lst):
        if not line.strip():
            continue
        unit, codon, meaning = line.split()[:3]
        path = f"{HOME}/results/stage2/{unit}_{codon}.txt"
        try:
            rows, _ = parse(path)
        except FileNotFoundError:
            print(f"{unit}\t{codon}\t{meaning}\tmissing")
            continue
        if codon not in rows or len(rows) < 2:
            print(f"{unit}\t{codon}\t{meaning}\tno_rows")
            continue
        ctl = [k for k in rows if k != codon][0]
        a, b = rows[codon], rows[ctl]
        split = ((num(a["gc3_p"]) < 1e-3 and abs(num(a["gc3_with"]) - num(a["gc3_without"])) > 0.05)
                 or num(a["use_dist"]) > 2 * num(b["use_dist"]) or num(a["disp"]) > 2 * num(b["disp"]))
        weak = num(a["cons_share"]) < 0.5 * num(b["cons_share"])
        flag = "split" if split else ("weak" if weak else "ok")
        print(f"{unit}\t{codon}\t{meaning}\t{ctl}\t" + "\t".join(f"{a.get(c, '')}\t{b.get(c, '')}" for c in COLS)
              + f"\t{flag}\t{organelle(path)}")


if __name__ == "__main__":
    main(sys.argv[1])
