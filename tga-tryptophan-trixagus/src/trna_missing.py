"""Side test for NCBI 'missing tRNA genes' genomes: which amino acids / codons lack any tRNA reader (ARAGORN via
pyaragorn), and are those codons still used in pyrodigal genes? A good assembly that uses a codon with no reader is
an unassigned-codon or unusual-reader lead; a fragmented one is just missing sequence.
Wobble model: a codon is 'read' if some tRNA's anticodon pairs with it under standard wobble (G:U at the third
codon position, anticodon I/U at position 34 treated as reading U/C/A or A/G); Ile2 CAU->AUA allowed for tRNA-Ile/Met.
usage: uv run --no-project --with pyrodigal,pyaragorn,requests python src/trna_missing.py <acc_list> -> results/trna/<acc>.json
"""
import json
import sys
from collections import Counter
from pathlib import Path

import pyaragorn
import pyrodigal

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "stop-codon-census" / "src"))
import readthrough as rt  # noqa: E402

rt.GEN = HERE.parent / "data" / "genomes"
OUT = HERE.parent / "results" / "trna"
COMP = {"A": "U", "C": "G", "G": "C", "U": "A", "T": "A"}


def reads(anticodon, codon):
    """anticodon 5'->3' (positions 34,35,36) pairs codon 5'->3' (1,2,3): 36:1, 35:2, 34:3 with wobble at 34:3."""
    a = anticodon.upper().replace("T", "U")
    c = codon.upper().replace("T", "U")
    if COMP[a[2]] != c[0] or COMP[a[1]] != c[1]:
        return False
    w, t = a[0], c[2]
    if a == "CAU" and c == "AUA":
        return True  # tRNA-Ile2 (lysidine / agmatidine C34); ARAGORN cannot tell it from tRNA-Met
    return (w == "G" and t in "UC") or (w == "C" and t == "G") or (w == "A" and t in "UCAG") or \
        (w == "U" and t in "AGU") or (w == "I" and t in "UCA")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for acc in [l.split()[0] for l in open(sys.argv[1]) if l.strip()]:
        if (OUT / f"{acc}.json").exists():
            continue
        try:
            cached = (rt.GEN / f"{acc}.fna.gz").exists()
            fn = rt.fetch(acc)
        except Exception as e:  # noqa: BLE001
            print(acc, "ERROR", e, flush=True)
            continue
        contigs = list(rt.read_fasta(fn))
        rna = pyaragorn.RNAFinder(translation_table=11, tmrna=False)
        trnas = []
        for _, s in contigs:
            for g in rna.find_rna(s.encode()):
                trnas.append((g.amino_acid, g.anticodon))
        total = sum(len(s) for _, s in contigs)
        finder = pyrodigal.GeneFinder(meta=total < 100000)
        if total >= 100000:
            finder.train(*[s.encode() for _, s in contigs], translation_table=11)
        use = Counter()
        for _, s in contigs:
            for g in finder.find_genes(s.encode()):
                nt = g.sequence()
                for i in range(3, len(nt) - 3, 3):
                    use[nt[i:i + 3]] += 1
        aa_present = sorted({a for a, _ in trnas})
        unread = {}
        for cod, n in use.items():
            aa = rt.CODE.get(cod, "X")
            if aa in "*X":
                continue
            if not any(reads(ac, cod) for a, ac in trnas if ac and len(ac) == 3):
                unread[cod] = (aa, n)
        tot = sum(use.values())
        res = dict(acc=acc, bp=total, contigs=len(contigs), n_trna=len(trnas), aa_with_trna=aa_present,
                   unread={c: dict(aa=a, n=n, per_1000=round(1000 * n / max(1, tot), 2)) for c, (a, n) in unread.items()})
        (OUT / f"{acc}.json").write_text(json.dumps(res))
        print(acc, total, len(contigs), "tRNAs", len(trnas), "unread", {c: v[0] + str(v[1]) for c, v in unread.items()},
              flush=True)
        if not cached:
            fn.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
