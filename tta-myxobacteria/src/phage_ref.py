"""Phage reference for pooled Arm B (PREREG.md, Changes 2026-11-10 21:31): MGV proteins (MGV's Prodigal calls,
partial=00 only) of the longest genome of every vOTU none of whose genomes MGV calls as reading TAG or TGA through,
minus every vOTU with a genome in any V-B pool (data/phage/pools/pools.tsv) -> data/phage/phageref.faa, .dmnd.

usage: curl -s https://portal.nersc.gov/MGV/MGV_v1.0_2021_07_08/mgv_proteins.faa | python -I src/phage_ref.py
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PH = ROOT / "data" / "phage"


def main():
    votu, length, recoded = {}, {}, set()
    with open(PH / "mgv" / "mgv_contig_info.tsv") as f:
        next(f)
        for line in f:
            r = line.rstrip("\n").split("\t")
            votu[r[0]], length[r[0]] = r[1], int(r[2])
            if r[9] != "NULL":
                recoded.add(r[1])
    pooled = set()
    with open(PH / "pools" / "pools.tsv") as f:
        next(f)
        for line in f:
            pooled |= {votu[g] for g in line.rstrip("\n").split("\t")[5].split(",") if g in votu}
    best = {}
    for g, v in votu.items():
        if v in recoded or v in pooled:
            continue
        if v not in best or (length[g], g) > (length[best[v]], best[v]):
            best[v] = g
    keep = set(best.values())
    print(f"{len(votu)} genomes, {len(set(votu.values()))} vOTUs; {len(recoded)} with a readthrough call, "
          f"{len(pooled)} in pools; reference: {len(keep)} genomes", flush=True)
    n = kept = 0
    on, seen = False, set()
    with open(PH / "phageref.faa", "w") as out:
        for line in sys.stdin:
            if line.startswith(">"):
                n += 1
                pid = line[1:].split()[0]
                on = pid.rsplit("_", 1)[0] in keep and "partial=00" in line
                if on:
                    kept += 1
                    seen.add(pid.rsplit("_", 1)[0])
                    out.write(f">{pid}\n")
            elif on:
                out.write(line.rstrip("\n").rstrip("*") + "\n")
    print(f"{n} proteins read, {kept} kept from {len(seen)} genomes", flush=True)
    if len(seen) < 0.99 * len(keep):
        sys.exit("protein stream looks truncated")
    subprocess.run(["diamond", "makedb", "--in", str(PH / "phageref.faa"), "-d", str(PH / "phageref"), "--threads",
                    "2", "--quiet"], check=True)


if __name__ == "__main__":
    main()
