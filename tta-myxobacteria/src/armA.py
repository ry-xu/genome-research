"""Arm A helpers shared by the calibration, control and GTDB analyses: build the null table from saved profiles,
score one profile under the table its genes were called with.
"""
from pathlib import Path

import numpy as np

import codonpos as cp
import score as sc

ROOT = Path(__file__).resolve().parent.parent


def sense_for(table):
    """Codon indices that are sense when genes were called with `table` (4: TGA sense; 16: TAG sense)."""
    stops = cp.TABLE_STOPS[table]
    return [i for i, c in enumerate(cp.CODONS) if c not in stops]


def gtdb_table(stops):
    """GTDB's calls: table 4 when (almost) no complete gene ends at TGA."""
    return 4 if stops[2] < 0.002 * max(1, stops.sum()) else 11


def load_gtdb(files):
    """-> acc, dom, pos [n,65,NB], stops [n,3], meta [n,3] concatenated over prof_*.npz files."""
    parts = [np.load(f) for f in files]
    cat = lambda k: np.concatenate([p[k] for p in parts])
    return cat("acc"), cat("dom"), cat("pos"), cat("stops"), cat("meta")


def null_from(gcs, poss):
    rho, irate = sc.null_table(list(zip(gcs, poss)))
    rho[:, cp.STOP_IDX] = np.nan  # a real stop that is sense under the calling table is scored with rho = 1
    irate[:, cp.STOP_IDX] = np.nan
    return rho, irate


def score(pos, gc, table, rho, irate):
    return sc.score(pos, gc, rho, irate, sense=sense_for(table))
