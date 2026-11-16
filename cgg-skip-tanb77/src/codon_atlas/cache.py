"""Genome cache: one gzipped FASTA per accession, fetched once, shared by every project.

Layout: genomes/<shard>/<acc>.fna.gz, where <shard> is the first 8 characters of an NCBI accession
("GCA_0004" for GCA_000433515.1) and "_other" for anything else (MGnify, SPIRE, QTP bins...).
index/genomes.tsv lists each cached genome with its source and the project file it was first linked from.
Existing project files are hard-linked in, so linking costs no disk and leaves the original path working.
"""
import fcntl
import gzip
import io
import os
import re
import time
import zipfile
import zlib
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(os.environ.get("CODON_ATLAS", Path(__file__).resolve().parents[2] / "data" / "atlas"))
GENOMES = ROOT / "genomes"
PRODUCTS = ROOT / "products"
INDEX = ROOT / "index"
REF = ROOT / "ref"
NCBI = re.compile(r"^(?:RS_|GB_)?(GC[AF]_\d{9}\.\d+)$")
INDEX_COLS = ["acc", "path", "bytes", "source", "linked_from", "added"]


def norm(acc):
    """Strip GTDB's RS_/GB_ prefix and any file suffix: 'RS_GCF_000016525.1' -> 'GCF_000016525.1'."""
    acc = re.sub(r"\.(fna|fa|fasta)(\.gz)?$", "", Path(acc).name)
    m = NCBI.match(acc)
    return m.group(1) if m else acc


def shard(acc):
    acc = norm(acc)
    return acc[:8] if NCBI.match(acc) else "_other"


def genome_path(acc):
    acc = norm(acc)
    return GENOMES / shard(acc) / f"{acc}.fna.gz"


def product_dir(acc):
    acc = norm(acc)
    return PRODUCTS / shard(acc) / acc


@contextmanager
def lock(acc):
    """Per-genome exclusive lock (flock), so parallel agents never compute or write the same genome twice."""
    d = product_dir(acc)
    d.mkdir(parents=True, exist_ok=True)
    with open(d / ".lock", "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


@contextmanager
def index_lock():
    INDEX.mkdir(parents=True, exist_ok=True)
    with open(INDEX / ".genomes.lock", "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def _append_index(rows):
    fn = INDEX / "genomes.tsv"
    with index_lock():
        new = not fn.exists()
        with open(fn, "a") as f:
            if new:
                f.write("\t".join(INDEX_COLS) + "\n")
            for r in rows:
                f.write("\t".join(str(r[c]) for c in INDEX_COLS) + "\n")


def gzip_ok(fn):
    """True if fn is a complete gzip stream whose first line is a FASTA header."""
    try:
        with gzip.open(fn, "rb") as f:
            if not f.read(1) == b">":
                return False
            while f.read(1 << 24):
                pass
        return True
    except (OSError, EOFError, zlib.error):
        return False


def link(src, acc=None, source="", check=True):
    """Hard-link an existing gzipped genome FASTA into the cache (no copy). Returns the cache path.
    If the accession is already cached, nothing is done (the first copy wins)."""
    src = Path(src)
    acc = norm(acc or src.name)
    dst = genome_path(acc)
    if dst.exists():
        return dst
    if not src.name.endswith(".gz"):
        raise ValueError(f"{src}: cache holds gzipped FASTA only; gzip it first or use add_fasta()")
    if check and not gzip_ok(src):
        raise ValueError(f"{src}: not a complete gzipped FASTA")
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dst)
    except FileExistsError:
        return dst
    _append_index([dict(acc=acc, path=dst.relative_to(ROOT), bytes=dst.stat().st_size, source=source,
                        linked_from=src, added=time.strftime("%Y-%m-%dT%H:%M:%S"))])
    return dst


def add_fasta(records_or_path, acc, source=""):
    """Write a new genome into the cache from an uncompressed FASTA path or (name, seq) records."""
    acc = norm(acc)
    dst = genome_path(acc)
    if dst.exists():
        return dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(f".tmp{os.getpid()}")
    with gzip.open(tmp, "wt", compresslevel=6) as out:
        if isinstance(records_or_path, (str, Path)):
            with open(records_or_path) as f:
                for line in f:
                    out.write(line)
        else:
            for name, seq in records_or_path:
                out.write(f">{name}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")
    os.replace(tmp, dst)
    _append_index([dict(acc=acc, path=dst.relative_to(ROOT), bytes=dst.stat().st_size, source=source,
                        linked_from="", added=time.strftime("%Y-%m-%dT%H:%M:%S"))])
    return dst


def _ftp_dir(acc):
    """The assembly's NCBI FTP directory URL (over https), or None if NCBI lists none."""
    import requests
    d = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[:3]}/{acc[4:7]}/{acc[7:10]}/{acc[10:13]}/"
    r = requests.get(d, timeout=60)
    sub = sorted(set(re.findall(rf'href="({re.escape(acc)}_[^"/]+)/?"', r.text)))
    return f"{d}{sub[0]}/{sub[0]}" if r.status_code == 200 and sub else None


def _fetch_ftp(acc, tmp):
    """NCBI FTP site over https (fallback when the Datasets API is down); from stop-codon-census readthrough.py."""
    import requests
    for attempt in range(4):
        if attempt:
            time.sleep(30 * attempt)
        try:
            d = _ftp_dir(acc)
            if not d:
                continue
            r = requests.get(f"{d}_genomic.fna.gz", timeout=600)
            if r.status_code == 200 and r.content[:2] == b"\x1f\x8b":
                tmp.write_bytes(r.content)
                return
            if r.status_code == 404 and requests.get(f"{d}_assembly_stats.txt", timeout=60).status_code == 200:
                break  # NCBI keeps the reports but no longer serves the sequence: no point retrying
        except requests.RequestException:
            pass
    raise RuntimeError(f"download failed: {acc}")


def _fetch_ena(acc, tmp):
    """ENA browser API: assemblies submitted through ENA whose sequence NCBI no longer serves (GTDB's UHGG MGYG-HGUT
    representatives, GC[AF]_9023xxxxx: the NCBI directory keeps the reports but no _genomic.fna.gz). A GCF accession
    is fetched as its GCA twin. The total length must match NCBI's assembly_stats.txt; headers become the sequence
    accession (ENA|X|X.1 -> X.1)."""
    import requests
    d = _ftp_dir(acc)
    r = requests.get(f"{d}_assembly_stats.txt", timeout=60) if d else None
    m = re.search(r"^all\tall\tall\tall\ttotal-length\t(\d+)", r.text if r is not None else "", re.M)
    if not m:
        raise RuntimeError(f"download failed: {acc} (and no NCBI assembly_stats to check an ENA copy against)")
    r = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/GCA{acc[3:]}?download=true", timeout=600)
    r.raise_for_status()
    lines = [(">" + x[1:].split()[0].split("|")[-1]) if x.startswith(">") else x for x in r.text.splitlines()]
    n = sum(len(x) for x in lines if not x.startswith(">"))
    if n != int(m.group(1)):
        raise RuntimeError(f"download failed: {acc} (ENA GCA{acc[3:]} has {n} bp, NCBI says {m.group(1)})")
    tmp.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode(), compresslevel=6))


def fetch(acc):
    """Path to the cached genome, downloading it from NCBI (Datasets API, then FTP), else ENA, if it is not cached yet.
    Only NCBI accessions can be fetched; register other genomes with link() or add_fasta()."""
    acc = norm(acc)
    dst = genome_path(acc)
    if dst.exists():
        return dst
    if not NCBI.match(acc):
        raise FileNotFoundError(f"{acc} is not cached and is not an NCBI accession; register it with link()")
    import requests
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(f".tmp{os.getpid()}")
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
           "?include_annotation_type=GENOME_FASTA")
    source = "ncbi-datasets"
    try:
        r = requests.get(url, timeout=600)
        r.raise_for_status()
        z = zipfile.ZipFile(io.BytesIO(r.content))
        names = [n for n in z.namelist() if n.endswith(".fna")]
        if not names:
            raise RuntimeError(f"no fna for {acc}")
        tmp.write_bytes(gzip.compress(z.read(names[0]), compresslevel=6))
    except Exception:
        source = "ncbi-ftp"
        try:
            _fetch_ftp(acc, tmp)
        except RuntimeError:
            source = "ena"
            _fetch_ena(acc, tmp)
    if not gzip_ok(tmp):
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"{acc}: downloaded file is not a complete gzipped FASTA")
    os.replace(tmp, dst)
    _append_index([dict(acc=acc, path=dst.relative_to(ROOT), bytes=dst.stat().st_size, source=source,
                        linked_from="", added=time.strftime("%Y-%m-%dT%H:%M:%S"))])
    return dst


def read_fasta(fn):
    """Yield (name, upper-case sequence) from a gzipped or plain FASTA."""
    op = gzip.open if str(fn).endswith(".gz") else open
    name, buf = None, []
    with op(fn, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf).upper()


def cached():
    """Accessions in the cache (from the directory tree, so it never disagrees with the files)."""
    return sorted(p.name[:-7] for p in GENOMES.glob("*/*.fna.gz"))
