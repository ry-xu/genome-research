"""Re-check, step 3a: list the count units for the main set (Amendment 4), one per line, for xargs.
Follow-up S copy (recheck_S/): runs from ../followS_runs.txt, done-check against recheck_S/counts/.
  .venv/bin/python -I units.py [CHUNK]      CHUNK = spots per unit (default 2,500,000)
Each line of ../followS_runs.txt is "dir 1 SPOTS 2 RUN"; spots 1..SPOTS of RUN are streamed (the same as
fastq-dump -N 1 -X SPOTS), cut into ranges of at most CHUNK spots. Output lines: "RUN --from A --to B".
Units already done (counts/<run>.<A>-<B>.tsv present and its summary has exit 0 and B - A + 1 spots) are left out,
so the same command resumes after an interruption.
"""
import os
import sys

R = "results/rna_cgg/recheck_S/"


def plan(chunk):
    units = []
    for line in open(R + "../followS_runs.txt"):
        x = line.split()
        if not x:
            continue
        assert len(x) == 5 and x[1] == "1" and x[3] == "2", line
        run, spots = x[4], int(x[2])
        for a in range(1, spots + 1, chunk):
            units.append((run, a, min(a + chunk - 1, spots)))
    return units


def done(run, a, b):
    stem = R + f"counts/{run}.{a}-{b}"
    if not (os.path.exists(stem + ".tsv") and os.path.exists(stem + ".summary")):
        return False
    meta = dict(x.rstrip("\n").split("\t") for x in open(stem + ".summary"))
    return meta["exit"] == "0" and int(meta["spots"]) == b - a + 1


def main():
    chunk = int(sys.argv[1]) if len(sys.argv) > 1 else 2_500_000
    for run, a, b in plan(chunk):
        if not done(run, a, b):
            print(f"{run} --from {a} --to {b}")


if __name__ == "__main__":
    main()
