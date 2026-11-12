# Shared ends: is TTA used like a stop across the clade, and is it old? (post hoc; Basil, 2026-11-11 10:55 EST)

Script: `shared_ends.py` (run: `uv run --no-project --with biopython python -I scratch/shared_ends/shared_ends.py
scratch/redteam scratch/shared_ends`), output `shared_ends.out`, rows `pairs.tsv`. Inputs are the red team's own gene
calls, RBH orthologs and sister-defined ends (`scratch/redteam/ortho/`); nothing new was frozen for this.

Method: an end of a clade (or sister) genome X is where its gene ends at the same place as its best sister ortholog
(TTA end: X runs on and the next codon is TTA; canonical end: X stops there too). For every other genome Y of the
same group with an RBH ortholog, the proteins are aligned (local, BLOSUM62) and Y's codon at the same place is read:
Y's stop if Y ends there, else Y's next codon. Pairs share genomes and ends, so the counts are descriptive, not
independent tests.

Clade X vs clade Y in another genus (resolved = X's last residue aligned to a residue of Y):

| X end | resolved | Y same codon | Y any canonical stop | Y TTA | Y other sense codon |
|---|---|---|---|---|---|
| TTA | 641 | 77 (12.0%) | 501 (78.2%) | 77 (12.0%) | 63 (9.8%) |
| TAA | 722 | 46 (6.4%) | 605 (83.8%) | 33 (4.6%) | 84 (11.6%) |
| TAG | 1,376 | 168 (12.2%) | 1,195 (86.8%) | 54 (3.9%) | 127 (9.2%) |
| TGA | 11,730 | 8,605 (73.4%) | 10,061 (85.8%) | 509 (4.3%) | 1,160 (9.9%) |

Same genus (JADLFU01 pair, JANWXX01 trio): TTA ends are TTA in Y at 79/239 (33.1%), a canonical stop at 155 (64.9%),
sense at 5 (2.1%); TAA ends stay TAA at 22.8%, TAG at 42.7%, TGA at 84.5%.

Sister X vs sister Y (negative control): Y has TTA at 3 of 38,417 resolved canonical ends (0.008%); in the clade,
777 of 18,007 (4.3%).

Reading:
- In the clade, TTA trades places with TAA, TAG and TGA at orthologous gene ends, at about the rate the canonical
  stops trade with each other; in the sisters it never does. TTA ends stay ends in the other genome as often as
  canonical ends do (about 90% across genera).
- 43 of 183 TTA-end genes with a resolved ortholog in another genus (23%) have TTA there too (corrected 11:37
  after fact-check 2: I first wrote 43 of 225, because pairs.tsv has CRLF line ends and my shell filter on the last
  column never matched "unresolved"; recount in `scratch/factcheck2/shared_ends_check.out`), and TTA ends are shared
  between every pair of the five genera at least once except JAXLYQ01 with CAYNEK01 and DASTPP01. That fits TTA ends
  being present in the clade's common ancestor; parallel TGA->TTA or TAA->TTA changes (one base each) at the same
  end cannot be ruled out site by site.
