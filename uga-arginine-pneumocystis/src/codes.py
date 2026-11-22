"""Genetic codes and the "blinded" translation used throughout.

Blinded table t: the declared NCBI table's meanings, except that every codon that is sense in the standard code keeps
a sense meaning (its declared one, or its standard one if the table declares it a stop). So stop-to-sense knowledge
(mito UGA = Trp) is used, sense-to-stop knowledge (vertebrate AGR, Scenedesmus UCA, ...) is erased.
"""
import itertools
from Bio.Data import CodonTable

CODONS = [''.join(c) for c in itertools.product('TCAG', repeat=3)]
CI = {c: i for i, c in enumerate(CODONS)}
STD_STOPS = ('TAA', 'TAG', 'TGA')
STD = CodonTable.unambiguous_dna_by_id[1].forward_table
SENSE = [c for c in CODONS if c not in STD_STOPS]  # the 61 standard sense codons

# Codon meanings not in Biopython's tables (literature codes used as controls)
EXTRA = {}

_cache = {}


def blinded(table):
    """-> (dict codon -> aa or '*') with sense-to-stop erased"""
    if table in _cache:
        return _cache[table]
    t = CodonTable.unambiguous_dna_by_id[table]
    m = {}
    for c in CODONS:
        if c in t.forward_table:
            m[c] = t.forward_table[c]
        elif c in STD:
            m[c] = STD[c]  # declared stop, standard sense: erase
        else:
            m[c] = '*'
    _cache[table] = m
    return m


def canonical_stops(table):
    """standard stop codons that the declared table also treats as stops"""
    m = blinded(table)
    return {c for c in STD_STOPS if m[c] == '*'}


def translate_ext(nt, cds_len, table):
    """Translate CDS + extension. Canonical stops before the annotated last codon -> 'X' (read-through, errors);
    the first canonical stop at or after the annotated last codon ends the protein.
    -> (protein, cut_codon_index or None if no stop in the extension)"""
    m = blinded(table)
    last = max(cds_len // 3 - 1, 0)
    out = []
    for i in range(len(nt) // 3):
        c = nt[3 * i:3 * i + 3]
        a = m.get(c, 'X')
        if a == '*':
            if i >= last:
                return ''.join(out), i
            a = 'X'
        out.append(a)
    return ''.join(out), None
