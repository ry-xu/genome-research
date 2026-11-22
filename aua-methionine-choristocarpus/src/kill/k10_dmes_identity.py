"""Kill check 5b: is the Genoscope 'Discosporangium' assembly really Discosporangium (the sister-genus comparison depends on it)?
GenBank Discosporangium markers (cox1 AB899176, LM995414; rbcL AB252654, AB252655) blastn against the dmes and cten mito and
plastid contigs; best identity per marker per target.
usage: uv run python -I checks/chor_wgs/kill/k10_dmes_identity.py DISCO_DIR"""
import subprocess, sys
D = sys.argv[1]
W = "scratchpad/chor_work"
T = {"dmes_mito": f"{W}/dmes_mito.fa", "cten_mito_WGS": f"{W}/cten_mito.fa", "Chor_NC_085326": f"{W}/NC_085326.fa",
     "dmes_plastid": f"{D}/dmes_plastid.fa", "cten_plastid_WGS": f"{D}/cten_plastid.fa", "Chor_plastid_NC_085297": f"{W}/NC_085297.fa"}
for name, fa in T.items():
    out = subprocess.run(["blastn", "-task", "dc-megablast", "-query", f"{D}/disco_markers.fa", "-subject", fa, "-evalue", "1e-20",
                          "-outfmt", "6 qseqid qlen pident length bitscore", "-max_hsps", "1"], check=True, capture_output=True, text=True).stdout
    best = {}
    for l in out.splitlines():
        q, ql, pid, L, bs = l.split("\t")
        if q not in best or float(bs) > best[q][3]: best[q] = (int(ql), float(pid), int(L), float(bs))
    print(name, " | ".join(f"{q}: {v[1]}% over {v[2]}/{v[0]}" for q, v in sorted(best.items())) or "no hits")
