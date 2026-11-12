"""Extract selected members from a remote zip via HTTP Range requests (no full download), writing each
as gzip into outdir. Usage: url outdir list|NAME [NAME ...]   (NAME matched against member basename without extension)"""
import sys, os, io, gzip, zipfile, urllib.request
class HTTPFile(io.RawIOBase):
    def __init__(self, url):
        self.orig = url; self._resolve()
        self.pos = 0; self.fetched = 0
    def _resolve(self):
        req = urllib.request.Request(self.orig, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            self.url = r.geturl(); self.size = int(r.headers["Content-Length"])
    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos
    def read(self, n=-1):
        if n is None or n < 0: n = self.size - self.pos
        if n == 0 or self.pos >= self.size: return b""
        end = min(self.pos + n, self.size) - 1
        for attempt in range(3):
            req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
            try:
                with urllib.request.urlopen(req, timeout=300) as r:
                    data = r.read()
                break
            except urllib.error.HTTPError as e:
                if e.code != 403 or attempt == 2: raise
                self._resolve()  # presigned redirect URL expired; re-resolve
        self.pos += len(data); self.fetched += len(data)
        return data
    def readinto(self, b):
        d = self.read(len(b)); b[:len(d)] = d; return len(d)
url, outdir, names = sys.argv[1], sys.argv[2], sys.argv[3:]
f = HTTPFile(url)
z = zipfile.ZipFile(io.BufferedReader(f, buffer_size=1 << 20))
infos = z.infolist()
print("zip size", f.size, "members", len(infos), file=sys.stderr)
if names == ["list"]:
    for i in infos: print(i.filename, i.file_size)
    sys.exit()
want = set(names)
for i in infos:
    base = os.path.basename(i.filename)
    stem = base.rsplit(".", 1)[0] if "." in base else base
    for suf in (".fa", ".fna", ".fasta"):
        if base.endswith(suf): stem = base[: -len(suf)]
    if stem in want or base in want:
        out = os.path.join(outdir, os.path.basename(base) + ".gz")
        with z.open(i) as src, gzip.open(out, "wb") as dst:
            dst.write(src.read())
        print("wrote", out, i.file_size, file=sys.stderr)
        want.discard(stem); want.discard(base)
print("missing:", sorted(want), "bytes fetched", f.fetched, file=sys.stderr)
