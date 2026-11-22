# Tweet fact-check: posts/roundup-1, tweet 2 (AGA = Ala in basidiomycete mitochondria)

Thread: ~/.local/share/lab-dashboard-v3/posts/roundup-1/thread.md (21 Nov 2026). Line numbers are of post.md in this folder on 21 Nov.
Figure: figs/roundup-1/public/aga_alanine.png = this folder's image.png, unchanged.

| Claim | Source | Status |
|---|---|---|
| about 30 basidiomycete fungi | l.1 ("about 30"), l.5 (32 genomes, about 29 species) | OK |
| Trichosporon, Tremella, Tilletia, Cystobasidium | l.5 | OK |
| AGA codes for alanine instead of arginine | l.1 ("most likely means alanine") | OK; post says "most likely". Protein evidence is one peptide (l.13, l.27) |
| first known fungal case | l.1 ("No known fungal genetic code reads AGA as alanine"); l.19 | OK |
| at conserved sites AGA aligns to alanine | l.5 (78% at well-conserved positions), l.7 (selection-free test: 38 Ala vs 4 Arg positions) | OK |
| 28 of the 32 have a tRNA for AGA with the alanine identity mark G3:U70 | l.5, l.9 | OK |
| none of the 4,322 arginine-reading genomes with an AGA tRNA has that mark | l.9 (4,322 put AGA at >= 3 Arg positions, none at Ala, and have a UCU tRNA; none has the pair) | OK |
| some green algae use the same route (Barcyte et al. 2025) | l.19 (Marsupiomonadales; Turmel 2020, Barcyte 2025; same tRNA route) | OK |

Status: ready; Mintaka 19:01 20 Nov recomputed all numbers; Sabik 21:44 20 Nov the peptide (README).

## Fresh-subagent recheck (21 Nov)
Applied: "codes for" -> "most likely codes for" (source l.1; one peptide); "the first known fungal case" -> "the first known case
of AGA as alanine in fungi" (Ashbya CUU/CUA = Ala is a prior fungal alanine reassignment, l.19). Figure numbers agree (4,387 in the
figure vs 4,322 in the tweet are different denominators, both correct).
