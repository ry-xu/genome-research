"""Six-frame search (standard code, stops kept as *) for a peptide pattern; report the codon under each pattern position.
usage: python pep_in_genome.py REGEX genome.fna [...]"""
import sys, re
from Bio.Seq import Seq
pat = re.compile(sys.argv[1])
for path in sys.argv[2:]:
    name, chunks = None, []
    def flush():
        if name is None:
            return
        s = ''.join(chunks).upper()
        for strand, ss in (('+', s), ('-', str(Seq(s).reverse_complement()))):
            for f in range(3):
                sub = ss[f:]
                sub = sub[:len(sub) // 3 * 3]
                aa = str(Seq(sub).translate(table=1))
                for m in pat.finditer(aa):
                    cods = [sub[3 * i:3 * i + 3] for i in range(m.start(), m.end())]
                    print(path.rsplit('/', 1)[-1][:16], name[:70], len(s), strand, f, m.start(), m.group(), ' '.join(cods), sep='\t')
    for line in open(path):
        if line.startswith('>'):
            flush(); name = line[1:].strip(); chunks = []
        else:
            chunks.append(line.strip())
    flush()
