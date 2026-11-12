"""Arm B: which codon sits where homologs end, in genes that run on past their homolog's C-terminus?

If sense codon X became a stop, a gene caller that reads X as sense calls genes ending at X as running on to the next
stop of its table. Against Swiss-Prot, such a gene aligns to its homolog's full length and then keeps going; the codon
just past the homolog's last residue is X far more often than chance.

Per gene (called with the given table, best DIAMOND hit vs Swiss-Prot among hits < 95% identical, so a genome
that is itself in Swiss-Prot is scored like any other; audit 2026-11-10):
  early-end event  the alignment covers >= 60% of the Swiss-Prot protein, ends within 3 residues of its C-terminus,
                   and the query protein continues >= 5 residues past the alignment end
  e*               query residue aligned to the homolog's last residue: qe + (slen - se)
  window           codons at query positions e*-1 .. e*+3 (offset +1 = where the gene would stop; +-2 for C-terminal
                   length slack); the stop codon of the called gene counts if it falls there
Per codon X (sense under the calling table):
  O        events with X in the window (each event counts X once)
  O_fam    families among them: events joined when they share a Swiss-Prot mnemonic or their proteins are >= 50%
           identical over >= 50% of both (DIAMOND all-vs-all of the genome's event proteins; numbered paralogs and
           IS copies such as INSI1/INSI3 are one family)
  E_tail   expected from X's frequency at e*+8 .. end of the called gene (same events; the run-on)
  E_int    expected from X's frequency at positions 21 .. qe-21 (same events; ordinary coding sequence)
           E = N_events * (1 - (1 - f)^5) with the larger f of the two
  ratio = O / E, log10p = Poisson upper tail

usage: [THREADS=4] [BATCH=1] [FAMDEDUP=1] [ARMB_DB=<dmnd>] python src/term_codon.py <out_prefix> <fasta[.gz]>:<table> | inject:<fasta>:<X>:<s>:<q>:<seed>:<table> ...
       (table g = GTDB's 11-or-4 rule, or a number; inject: recodes the genome as inject.py does first;
        FAMDEDUP=1: one event per family, for pooled genomes;
        ARMB_DB: reference other than Swiss-Prot, e.g. data/giantvirus/gvref.dmnd; its ids are their own family keys)
  -> <out_prefix>.tsv (one row per genome x codon), <out_prefix>_events.tsv (one row per event)
"""
import gzip
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codonpos as cp
from inject import call, rc, recode, synonyms
from score import poisson_sf_log10

DB = Path(os.environ.get("ARMB_DB") or Path(__file__).resolve().parent.parent.parent / "tag-bacteria" / "data" /
          "sprot.dmnd")
COVER, SLACK, RUNON, W_LO, W_HI, TAIL0, INT_PAD = 0.6, 3, 5, -1, 3, 8, 20
SELF_PID, FAM_PID, FAM_COV, MIN_TRAIN = 95.0, 50, 50, 20000
# codon meanings that differ from table 11 (only what calling with these tables needs)
TABLE_AA = {4: {"TGA": "W"}, 25: {"TGA": "G"}, 16: {"TAG": "L"}, 15: {"TAG": "Q"}, 22: {"TAG": "L", "TCA": "*"},
            23: {"TTA": "*"}}


def read_any(fn):
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


def code(table):
    aa = {c: cp.AA[i] for i, c in enumerate(cp.CODONS)}
    aa.update(TABLE_AA.get(table, {}))
    return aa


def call_genes(contigs, table):
    """-> list of (gene_id, nt incl. stop, complete), and the table used. A forced table always means single mode
    trained with that table (meta mode's models would ignore it; audit 2026-11-10), so inputs < 20 kb are refused
    (no genes) rather than called with another code."""
    total = sum(len(s) for _, s in contigs)
    if table == "g":
        calls, used = call(contigs), None
    elif total < MIN_TRAIN:
        print(f"refused: {total} bp < {MIN_TRAIN} bp, too small to train table {table}", file=sys.stderr, flush=True)
        return [], int(table)
    else:
        f = pyrodigal.GeneFinder()
        f.train(*[s.encode() for _, s in contigs], translation_table=int(table))
        calls = [(n, g.begin, g.end, g.strand, g.partial_begin, g.partial_end) for n, s in contigs
                 for g in f.find_genes(s.encode())]
        used = f.training_info.translation_table
    seqs = dict(contigs)
    genes = []
    for n, b, e, st, pb, pe in calls:
        nt = seqs[n][b - 1:e]
        genes.append((f"{n}_{b}_{e}_{st}", nt if st == 1 else rc(nt), not (pb or pe)))
    if used is None:  # recover which table the GTDB rule picked from the stops genes end on
        ends = [g[1][-3:] for g in genes if g[2]]
        used = 4 if ends.count("TGA") < 0.002 * max(1, len(ends)) else 11
    return genes, used


def translate(nt, aa):
    return "".join(aa.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def events(genes, aa, threads):
    """-> list of (gene_id, codons list, ql, qe, e_star, sseqid, pident, family)."""
    return events_batch([(genes, aa)], threads)[0]


def events_batch(sets, threads):
    """sets = [(genes, aa), ...] -> one events list per set, from a single DIAMOND run over all of them (DIAMOND
    reports each query on its own, so batching only saves the per-run start-up; scratch/batch_check.sh)."""
    nts, prots = {}, {}
    with tempfile.TemporaryDirectory() as td:
        q, h = Path(td) / "q.faa", Path(td) / "h.tsv"
        with open(q, "w") as f:
            for i, (genes, aa) in enumerate(sets):
                for gid, nt, complete in genes:
                    if not complete or len(nt) % 3:
                        continue
                    p = translate(nt[:-3], aa)
                    if len(p) < 50:
                        continue
                    nts[f"{i}|{gid}"], prots[f"{i}|{gid}"] = nt, p
                    f.write(f">{i}|{gid}\n{p}\n")
        subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(DB), "-o", str(h), "--threads", str(threads),
                        "--evalue", "1e-5", "--max-target-seqs", "10", "--max-hsps", "1", "--quiet", "--outfmt", "6",
                        "qseqid", "sseqid", "pident", "qstart", "qend", "sstart", "send", "bitscore", "qlen", "slen"],
                       check=True)
        best = {}
        for line in open(h):
            r = line.rstrip("\n").split("\t")
            if float(r[2]) >= SELF_PID:
                continue
            if r[0] not in best or float(r[7]) > float(best[r[0]][7]):
                best[r[0]] = r
    out = [[] for _ in sets]
    for key, r in best.items():
        qs, qe, ss, se, ql, sl = (int(x) for x in (r[3], r[4], r[5], r[6], r[8], r[9]))
        if se < sl - SLACK or (se - ss + 1) < COVER * sl or ql - qe < RUNON:
            continue
        nt = nts[key]
        i, gid = key.split("|", 1)
        cods = [nt[k:k + 3] for k in range(0, len(nt), 3)]  # index k = query residue k+1; index ql = the stop
        out[int(i)].append([gid, cods, ql, qe, qe + (sl - se), r[1], float(r[2]), None])
    families(out, prots, threads)
    return [[tuple(ev) for ev in evs] for evs in out]


def families(out, prots, threads):
    """Set each event's family (in place): union of events sharing a Swiss-Prot mnemonic and of event proteins
    >= FAM_PID % identical over >= FAM_COV % of both, within one genome or pool (one DIAMOND all-vs-all per batch;
    hits between different sets are ignored). Family label = the smallest mnemonic in it."""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    for i, evs in enumerate(out):
        for ev in evs:
            union(("g", i, ev[0]), ("m", i, mnemonic(ev[5])))
    if any(out):
        with tempfile.TemporaryDirectory() as td:
            q, d, h = Path(td) / "e.faa", Path(td) / "e", Path(td) / "e.tsv"
            with open(q, "w") as f:
                for i, evs in enumerate(out):
                    for ev in evs:
                        f.write(f">{i}|{ev[0]}\n{prots[f'{i}|{ev[0]}']}\n")
            subprocess.run(["diamond", "makedb", "--in", str(q), "-d", str(d), "--quiet"], check=True)
            subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(d), "-o", str(h), "--threads", str(threads),
                            "--quiet", "--id", str(FAM_PID), "--query-cover", str(FAM_COV), "--subject-cover",
                            str(FAM_COV), "--max-target-seqs", "0", "--outfmt", "6", "qseqid", "sseqid"], check=True)
            for line in open(h):
                a, b = line.split()
                (ia, ga), (ib, gb) = a.split("|", 1), b.split("|", 1)
                if ia == ib and ga != gb:
                    union(("g", int(ia), ga), ("g", int(ib), gb))
    label = {}
    for i, evs in enumerate(out):
        for ev in evs:
            r = find(("g", i, ev[0]))
            label[r] = min(label.get(r, "~"), mnemonic(ev[5]))
    for i, evs in enumerate(out):
        for ev in evs:
            ev[7] = label[find(("g", i, ev[0]))]


def mnemonic(sid):
    """Swiss-Prot entry name mnemonic: sp|P27725|AK1H_SERMA -> AK1H. Families (see families()) join events by it and
    by protein identity, so a repeated element such as an IS transposase cannot make a codon pile up on its own.
    A non-UniProt id (ARMB_DB, e.g. genome|contig_7) is its own key."""
    if not sid.startswith(("sp|", "tr|")):
        return sid
    return sid.split("|")[-1].split("_")[0]


def tally(evs, sense):
    """-> {codon: (O, O_fam, E, E_tail, E_int, ratio, log10p)}, n_events."""
    n = len(evs)
    sset = set(sense)
    win = {c: 0 for c in sense}
    fams = {c: set() for c in sense}
    tail = {c: 0 for c in sense}
    inte = {c: 0 for c in sense}
    nt_tail = nt_int = 0
    for gid, cods, ql, qe, es, sid, pid, fam in evs:
        seen = set(cods[k - 1] for k in range(es + W_LO, es + W_HI + 1) if 1 <= k <= ql + 1)
        for c in seen & sset:
            win[c] += 1
            fams[c].add(fam)
        for k in range(es + TAIL0, ql + 1):
            c = cods[k - 1]
            if c in tail:
                tail[c] += 1
            nt_tail += 1
        for k in range(1 + INT_PAD, qe - INT_PAD):
            c = cods[k - 1]
            if c in inte:
                inte[c] += 1
            nt_int += 1
    w = W_HI - W_LO + 1
    res = {}
    for c in sense:
        ft = tail[c] / nt_tail if nt_tail else 0.0
        fi = inte[c] / nt_int if nt_int else 0.0
        et, ei = n * (1 - (1 - ft) ** w), n * (1 - (1 - fi) ** w)
        e = max(et, ei, 0.5)
        res[c] = (win[c], len(fams[c]), e, et, ei, win[c] / e, poisson_sf_log10(win[c], e))
    return res, n


def prepare(fn, table, inj=None):
    """-> (genes, table used). inj = (X, s, q, seed): recode the genome first, exactly as inject.py does."""
    contigs = list(read_any(fn))
    if inj:
        X, s, q, seed = inj
        base = call(contigs)
        contigs, _, _ = recode(contigs, base, X, s, q, random.Random(seed), synonyms(contigs, base))
    return call_genes(contigs, table)


def one_per_family(evs):
    """Pooled genomes (V-B): keep the first event of each family (families() was run on the whole pool), by contig
    name then position, so near-identical genomes in a pool cannot repeat an event (O then equals O_fam)."""
    def key(ev):
        contig, b, _, _ = ev[0].rsplit("_", 3)
        return contig, int(b)
    out, seen = [], set()
    for ev in sorted(evs, key=key):
        if ev[7] not in seen:
            seen.add(ev[7])
            out.append(ev)
    return out


def finish(label, table, genes, used, evs, famdedup=False):
    if famdedup:
        evs = one_per_family(evs)
    aa = code(used)
    sense = [c for c in cp.CODONS if aa[c] != "*"]
    res, n = tally(evs, sense)
    ncomp = sum(1 for g in genes if g[2])
    return dict(label=label, table=table, used=used, genes=len(genes), complete=ncomp, events=n, res=res, evs=evs)


def run(label, fn, table, threads=4, inj=None):
    genes, used = prepare(fn, table, inj)
    return finish(label, table, genes, used, events(genes, code(used), threads))


def parse(a):
    """'<fasta>:<table>' or 'inject:<fasta>:<X>:<s>:<q>:<seed>:<table>' -> (label, fasta, table, inj)."""
    if a.startswith("inject:"):
        _, fn, X, s, q, seed, table = a.split(":")
        inj = (X, float(s), float(q), int(seed))
    else:
        (fn, table), inj = a.rsplit(":", 1), None
    label = Path(fn).name.split(".fna")[0].split(".fa")[0] + (f"|{X}:{s}:{q}:{seed}" if inj else "")
    return label, fn, table, inj


def main(prefix, args, threads=4, batch=1, famdedup=False):
    rows = open(f"{prefix}.tsv", "w")
    rows.write("label\ttable\tused\tgenes\tcomplete\tevents\tcodon\taa\tO\tO_fam\tE\tE_tail\tE_int\tratio\tlog10p\t"
               "rank\n")
    evf = open(f"{prefix}_events.tsv", "w")
    evf.write("label\ttable\tgene\tql\tqe\te_star\tsseqid\tpident\twindow\tfamily\n")
    for b in range(0, len(args), batch):
        jobs = [parse(a) for a in args[b:b + batch]]
        prep = [prepare(fn, table, inj) for _, fn, table, inj in jobs]
        evss = events_batch([(genes, code(used)) for genes, used in prep], threads)
        for (label, fn, table, inj), (genes, used), evs in zip(jobs, prep, evss):
            r = finish(label, table, genes, used, evs, famdedup)
            aa = code(r["used"])
            order = sorted(r["res"], key=lambda c: -r["res"][c][5])
            for rank, c in enumerate(order, 1):
                O, Of, E, et, ei, ratio, lp = r["res"][c]
                rows.write(f"{label}\t{table}\t{r['used']}\t{r['genes']}\t{r['complete']}\t{r['events']}\t{c}\t"
                           f"{aa[c]}\t{O}\t{Of}\t{E:.2f}\t{et:.2f}\t{ei:.2f}\t{ratio:.3f}\t{lp:.2f}\t{rank}\n")
            for gid, cods, ql, qe, es, sid, pid, fam in r["evs"]:
                w = ",".join(cods[k - 1] for k in range(es + W_LO, es + W_HI + 1) if 1 <= k <= ql + 1)
                evf.write(f"{label}\t{table}\t{gid}\t{ql}\t{qe}\t{es}\t{sid}\t{pid}\t{w}\t{fam}\n")
            top = order[0]
            print(label, table, r["used"], "events", r["events"], "top", top, r["res"][top][0],
                  f"{r['res'][top][5]:.1f}", f"{r['res'][top][6]:.1f}", flush=True)
        rows.flush()
        evf.flush()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:], int(os.environ.get("THREADS", 4)), int(os.environ.get("BATCH", 1)),
         os.environ.get("FAMDEDUP") == "1")
