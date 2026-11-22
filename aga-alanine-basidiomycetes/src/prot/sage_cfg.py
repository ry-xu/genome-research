"""Write a Sage config. usage: python3 src/sage_cfg.py <out.json> <mode> [source_residue]
mode: 'std' (closed search) or 'sub' (closed search + every single-residue substitution FROM source_residue
as a variable mod, max one variable mod per peptide). Instrument: high-res MS2 (Orbitrap HCD) unless
PCC_LOWRES=1 in the environment. TMT=<mass> in the environment adds a static label on K and peptide N-termini
(229.162932 TMT6/10/11, 304.207146 TMTpro). Isotope errors are off so substitution deltas stay unambiguous."""
import json, os, sys

MASS = dict(G=57.021464, A=71.037114, S=87.032028, P=97.052764, V=99.068414, T=101.047679, C=103.009185,
            L=113.084064, N=114.042927, D=115.026943, Q=128.058578, K=128.094963, E=129.042593,
            M=131.040485, H=137.058912, F=147.068414, R=156.101111, Y=163.06332, W=186.079313)
CAM = 57.021464
AAS = 'GASPVTCLNDQKEMHFRYW'  # I is searched as L (same mass); I in the DB is converted to L by the tally
TMT = float(os.environ.get('TMT', 0))
STATIC = {'C': CAM, **({'K': TMT} if TMT else {})}  # residue mods carried by every copy, incl. substituted ones

def deltas(src):
    """substitution mass deltas from residue src (as stored, with its static mod) to every other residue (with its
    static mod: a substituted C is alkylated, a substituted K is TMT-labelled)."""
    m_src = MASS[src] + STATIC.get(src, 0)
    out = {}
    for y in AAS:
        if y == src: continue
        d = MASS[y] + STATIC.get(y, 0) - m_src
        out[y] = round(d, 6)
    return out

def cfg(mode, src=None):
    lowres = os.environ.get('PCC_LOWRES') == '1'
    var = {'M': [15.994915]}
    if mode == 'sub':
        d = list(deltas(src).values())
        var[src] = sorted(set(var.get(src, []) + d))
        if src == 'L': var['I'] = var['L']  # I and L share the deltas (same mass)
    return {
        'database': {
            'bucket_size': 16384,
            'enzyme': {'missed_cleavages': 2, 'min_len': 7, 'max_len': 30, 'cleave_at': 'KR', 'restrict': 'P',
                       'c_terminal': True},
            'fragment_min_mz': 150.0, 'fragment_max_mz': 2000.0,
            'peptide_min_mass': 600.0, 'peptide_max_mass': 4500.0,
            'ion_kinds': ['b', 'y'], 'min_ion_index': 2,
            'static_mods': {**STATIC, **({'^': TMT} if TMT else {})}, 'variable_mods': var, 'max_variable_mods': 1,
            'decoy_tag': 'rev_', 'generate_decoys': True,
        },
        'precursor_tol': {'ppm': [-10, 10]},
        'fragment_tol': {'da': [-0.4, 0.4]} if lowres else {'ppm': [-20, 20]},
        'precursor_charge': [2, 5], 'isotope_errors': [0, 0],
        'deisotope': not lowres, 'chimera': False, 'wide_window': False, 'predict_rt': True,
        'min_peaks': 15, 'max_peaks': 150, 'max_fragment_charge': 1 if lowres else 2,
        'min_matched_peaks': 4, 'report_psms': 1,
    }

if __name__ == '__main__':
    out, mode = sys.argv[1], sys.argv[2]
    src = sys.argv[3] if len(sys.argv) > 3 else None
    json.dump(cfg(mode, src), open(out, 'w'), indent=1)
