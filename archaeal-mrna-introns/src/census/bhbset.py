"""BHB label of a set of candidate removals, microhomology k, and same-shape decoys (RULES.md R3-R5)."""
import sys
sys.path.insert(0, '/Users/ryan/lab/work/archaeal-intron-census/src/checker')
from bhb import bhb_score, is_strict, is_relaxed

EDGE = 100

def slide_set(s, x, y):
    """All (x', y') giving the same spliced product as (x, y), and k = number of slides."""
    left = 0
    while x - 1 - left >= 0 and s[x - 1 - left] == s[y - 1 - left]:
        left += 1
    right = 0
    while y + right < len(s) and s[x + right] == s[y + right]:
        right += 1
    return [(x + t, y + t) for t in range(-left, right + 1)], left + right

def label(s, pairs):
    """(any strict-or-relaxed, any strict) over the removals."""
    anyb = anys = False
    for x, y in pairs:
        sc = bhb_score(s, x, y)
        if is_strict(sc):
            return True, True
        if is_relaxed(sc):
            anyb = True
    return anyb, anys

def ok_window(s, pairs):
    lo = min(x for x, _ in pairs); hi = max(y for _, y in pairs)
    return lo - 20 >= EDGE and hi + 20 <= len(s) - EDGE and 'N' not in s[lo - 20:hi + 20]

def decoys(s, pairs, L):
    """[(delta, bhb_any, bhb_strict)] for the kept shifted copies of the whole set."""
    out = []
    for m in (30, 90, 150, 210):
        for sign in (1, -1):
            dl = sign * (L + m)
            sh = [(x + dl, y + dl) for x, y in pairs]
            if ok_window(s, sh):
                a, b = label(s, sh)
                out.append((dl, a, b))
    return out
