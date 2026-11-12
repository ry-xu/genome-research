"""Arm A over every GTDB R232 species representative, streamed (gtdb_proteins_nt_reps_r232.tar.gz, 187 GB, never
written to disk; GTDB's own Prodigal calls). Per genome: codonpos.profile() (codon counts by distance from the gene
end) and release-factor candidates (280-480 aa proteins with GGQ) searched against tag-bacteria's markers.dmnd
(RF1, RF2, aRF1 + decoys); RF1/RF2/aRF1 best hits are kept with their sequences for the motif screen (Arm C).

The download runs as K parallel HTTP range requests reassembled in order (one stream gets ~8 MB/s); members are
profiled by a process pool; DIAMOND runs on batches in a background thread.

usage: nice -n 10 uv run --no-project --with numpy,requests python -I src/gtdb_scan.py [K=4] [workers=4] [limit]
  -> results/gtdb/prof_NNNN.npz (CHUNK genomes each), results/gtdb/rf_hits.tsv, results/gtdb/scan.log
A restart re-reads the stream from the start and skips genomes already saved.
"""
import gzip
import os
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import codonpos as cp

URL = "https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_nt_reps_r232.tar.gz"
OUT = Path(os.environ.get("SCAN_OUT", ROOT / "results" / "gtdb"))
MARKERS = ROOT.parent / "tag-bacteria" / "data" / "markers.dmnd"
CHUNK = 2000
RFKEEP = ("RF1", "RF2", "aRF1")


class RangeReader:
    """File-like sequential reader over K concurrent HTTP range requests of `size` bytes, delivered in order."""

    def __init__(self, url, k=4, size=32 << 20):
        self.url, self.size = url, size
        self.total = int(requests.head(url, timeout=60).headers["content-length"])
        self.nchunks = (self.total + size - 1) // size
        self.slots = {}
        self.cv = threading.Condition()
        self.next_fetch = 0
        self.cur, self.cur_i, self.pos = b"", -1, 0
        self.ahead = 2 * k
        self.threads = [threading.Thread(target=self._worker, daemon=True) for _ in range(k)]
        for t in self.threads:
            t.start()

    def _fetch(self, i):
        lo = i * self.size
        hi = min(self.total, lo + self.size) - 1
        for attempt in range(200):
            try:
                r = requests.get(self.url, headers={"Range": f"bytes={lo}-{hi}"}, timeout=120)
                if r.status_code == 206 and len(r.content) == hi - lo + 1:
                    return r.content
            except requests.RequestException:
                pass
            time.sleep(min(60, 3 * (attempt + 1)))
        raise RuntimeError(f"range {i} failed")

    def _worker(self):
        while True:
            with self.cv:
                while self.next_fetch < self.nchunks and self.next_fetch > self.cur_i + self.ahead:
                    self.cv.wait()
                if self.next_fetch >= self.nchunks:
                    return
                i = self.next_fetch
                self.next_fetch += 1
            data = self._fetch(i)
            with self.cv:
                self.slots[i] = data
                self.cv.notify_all()

    def read(self, n=-1):
        out = []
        while n != 0:
            if self.pos >= len(self.cur):
                if self.cur_i + 1 >= self.nchunks:
                    break
                with self.cv:
                    while self.cur_i + 1 not in self.slots:
                        self.cv.wait()
                    self.cur_i += 1
                    self.cur = self.slots.pop(self.cur_i)
                    self.pos = 0
                    self.cv.notify_all()
            take = len(self.cur) - self.pos if n < 0 else min(n, len(self.cur) - self.pos)
            out.append(self.cur[self.pos:self.pos + take])
            self.pos += take
            if n > 0:
                n -= take
        return b"".join(out)


def work(acc, dom, raw):
    text = gzip.decompress(raw).decode("ascii", "replace")
    p = cp.profile(cp.prodigal_genes(text))
    return acc, dom, p


def diamond(rows):
    """rows: [(acc, gid, prot)] -> [(acc, gid, kind, pident, len, prot)] for RF1/RF2/aRF1 best hits."""
    if not rows:
        return []
    seqs = {f"{a}|{g}": p for a, g, p in rows}
    with tempfile.TemporaryDirectory() as td:
        q, h = Path(td) / "q.faa", Path(td) / "h.tsv"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(MARKERS), "-o", str(h), "--threads", "2",
                        "--quiet", "--evalue", "1e-20", "--max-target-seqs", "1", "--outfmt", "6", "qseqid",
                        "sseqid", "pident", "length"], check=True)
        out, seen = [], set()
        for line in open(h):
            qid, sid, pid, ln = line.rstrip("\n").split("\t")
            kind = sid.split("|")[0]
            if qid in seen:
                continue
            seen.add(qid)
            if kind in RFKEEP:
                a, g = qid.split("|", 1)
                out.append((a, g, kind, pid, ln, seqs[qid]))
    return out


def save_chunk(idx, items):
    accs = [a for a, _, _ in items]
    np.savez_compressed(OUT / f"prof_{idx:04d}.npz", acc=np.array(accs), dom=np.array([d for _, d, _ in items]),
                        pos=np.stack([p["pos"] for _, _, p in items]).astype(np.int32),
                        g30=np.stack([p["g30"] for _, _, p in items]).astype(np.int32),
                        g100=np.stack([p["g100"] for _, _, p in items]).astype(np.int32),
                        gany=np.stack([p["gany"] for _, _, p in items]).astype(np.int32),
                        stops=np.stack([p["stops"] for _, _, p in items]).astype(np.int32),
                        meta=np.array([[p["ngenes"], p["ncomplete"], p["gc"]] for _, _, p in items]))


def main(k=4, workers=4, limit=None):
    OUT.mkdir(parents=True, exist_ok=True)
    done, nchunk = set(), 0
    for f in sorted(OUT.glob("prof_*.npz")):
        done |= set(np.load(f)["acc"].tolist())
        nchunk = max(nchunk, int(f.stem.split("_")[1]) + 1)
    log = open(OUT / "scan.log", "a")
    rfout = open(OUT / "rf_hits.tsv", "a")
    rr = RangeReader(URL, k=k)
    tar = tarfile.open(fileobj=rr, mode="r|gz")
    pool = ProcessPoolExecutor(workers)
    dpool = ThreadPoolExecutor(1)
    pending, items, rfrows, dfuts = [], [], [], []
    t0, seen = time.time(), 0

    def take_one():
        nonlocal items, rfrows, nchunk
        acc, dom, p = pending.pop(0).result()
        items.append((acc, dom, p))
        rfrows += [(acc, g, s) for g, s in p["rf"]]
        if len(items) >= CHUNK:
            save_chunk(nchunk, items)
            nchunk += 1
            items = []
            rows, rfrows = rfrows, []
            dfuts.append(dpool.submit(diamond, rows))
            print(f"{seen} genomes, chunk {nchunk}, {time.time() - t0:.0f}s, {rr.cur_i * rr.size / 1e9:.1f} GB",
                  file=log, flush=True)
        while dfuts and dfuts[0].done():
            for r in dfuts.pop(0).result():
                rfout.write("\t".join(map(str, r)) + "\n")
            rfout.flush()

    for m in tar:
        if not m.isfile() or not m.name.endswith("_protein.fna.gz"):
            continue
        acc = Path(m.name).name.replace("_protein.fna.gz", "")
        dom = "archaea" if "/archaea/" in m.name else "bacteria"
        raw = tar.extractfile(m).read()
        seen += 1
        if acc in done:
            continue
        pending.append(pool.submit(work, acc, dom, raw))
        while len(pending) > 4 * workers:
            take_one()
        while pending and pending[0].done():
            take_one()
        if limit and seen >= limit:
            break
    while pending:
        take_one()
    if items:
        save_chunk(nchunk, items)
        dfuts.append(dpool.submit(diamond, rfrows))
    for f in dfuts:
        for r in f.result():
            rfout.write("\t".join(map(str, r)) + "\n")
    rfout.close()
    print(f"done: {seen} genomes, {time.time() - t0:.0f}s", file=log, flush=True)
    pool.shutdown()
    dpool.shutdown()
    os._exit(0)  # range threads may still hold sockets


if __name__ == "__main__":
    a = sys.argv[1:]
    main(int(a[0]) if a else 4, int(a[1]) if len(a) > 1 else 4, int(a[2]) if len(a) > 2 else None)
