"""RULES R2 species check: share of streamed reads mapping (MAPQ>=20, primary) at >=95% identity
(1 - NM / aligned query length). usage: python species_check.py prefix [reads_per_spot=2]   (reads prefix.bam, prefix.log)"""
import sys, re, subprocess
p = sys.argv[1]; rps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
spots = int(re.search(r'Read (\d+) spots', open(p + '.log').read()).group(1))
n = 0
for line in subprocess.run(['samtools', 'view', p + '.bam'], capture_output=True, text=True).stdout.splitlines():
    f = line.split('\t')
    nm = next(int(t[5:]) for t in f[11:] if t.startswith('NM:i:'))
    al = sum(int(k) for k, op in re.findall(r'(\d+)([MI])', f[5]))
    n += nm <= 0.05 * al
print(f'{p}\tspots {spots}\treads {rps * spots}\tmito_95 {n}\tshare {n / (rps * spots):.4f}\tpass {n / (rps * spots) >= 0.01}')
