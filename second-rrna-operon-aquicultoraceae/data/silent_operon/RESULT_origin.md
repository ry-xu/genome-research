<!-- Copied from the lab record (work/silent-operon/, 22 Nov 2026). Paths such as out/, src/, recheck/ or frozen-core/ refer to the original working folder; the small tables from out/ are in this folder (see ../README.md), the rest is not included. -->
# T1 Origin: gene trees of divergent and normal 16S/23S (B_origin, 22 Nov)

**Finding.** The divergent operon is one old lineage that is older than the Aquicultoraceae genera and is also
present in the sister family JAMCWL01. Inside the family it follows the species tree. No donor is visible: every
outside sequence is about as far from it as its own genome's normal copy is, and where its long branch attaches is
unresolved.

Numbers (rules in `RULES.md`, T1):
- **One clade.** Divergent copies are monophyletic in both trees (UFBoot 100; 23S 13 copies / 5 genera, 16S 13 copies / 6 genera).
  The stem branch is long: 0.34 subs/site (23S), 0.23 (16S); within the clade, root-to-tip 0.13 vs 0.25 for the family normal copies' tree.
- **Present in the sister family.** GCA_023481285.1 (JAMCWL01, the only other family of order Aquicultorales in GTDB R232)
  has two 16S: one normal (sits with the family's normal copies) and one (JAMCWL010000101.1:766-2426 -) that is sister to the
  Aquicultoraceae divergent clade (UFBoot 100). It is 74-75% identical to the divergent 16S vs 70.5% to normal 16S, and
  deviates at 13.4% of census-invariant columns (divergent copies 13-14%, normal 0-2%). Recheck (fresh agent, own code):
  it carries 24/41 (cmalign) or 29/52 (mafft) of the changes at >= 99%-conserved columns that the family's divergent
  copies share, against 1-6 for the long-branched Tenuifilum divergent copy, so this is shared ancestry, not long-branch
  attraction. The contig is not contamination (coverage 3.9x vs genome 4.1x, GC 55.4 vs 55.9%). A 23S fragment on the
  same contig that I had counted as support was **killed** by the recheck: over aligned columns it is not closer to the
  divergent type, and it sits on the opposite strand (likely a misjoin). The JAMCWL01 copy also does not share the
  family's divergent-locus flanking genes. So the divergent lineage was in the common ancestor of both families, or it
  moved between their ancestors before the Aquicultoraceae genera split; these two cannot be told apart.
- **No donor.** Best outside hits of the divergent 23S: 73.5-77.6% over 1.3-2.4 kb, scattered over Bacillota, Sysuimicrobiota,
  Armatimonadota, Hydrothermota, Chloroflexota, Pseudomonadota, Actinomycetota (no lineage recurs). Divergent 16S: 71-75%
  (Desulfobacterota, Actinomycetes, Bacillota...). Local blast favours short conserved stretches, so I compared over
  aligned positions (cmalign columns): each divergent 23S is 0.708-0.717 identical to the best family normal copy,
  0.706-0.719 to the best other-family normal of the class, and 0.721-0.731 to the best outside sequence (a Sysuimicrobiota
  MAG, 4 of 5). 16S: 0.695-0.711 to family normal, 0.702-0.720 to the best outside. The divergent copy is about equally
  far from everything, so no lineage stands out as a donor. The HGT rule (outside >= own normal + 3 points) is not met;
  the 23S excess is 0.4-1.5 points. Normal copies' best outside hits are, as expected, Actinomycetota
  (Geothermincolia, Coriobacteriia, 82-89%). Phage databases were not searched at nucleotide level (no hit at protein level in T2).
- **Deep placement unresolved.** 23S ML tree: the divergent clade attaches among long-branched outside hits (UFBoot 29-43),
  not inside class Aquicultoria; AU test rejects "divergent + Aquicultorales normals form a clade" (p = 0.012). The
  "+ all Aquicultoria normals" constrained search scored worse than the narrower constraint, so that search failed and its
  p-value (0.0007) is not trusted. 16S: "inside Aquicultorales" p = 0.042, "inside class Aquicultoria" not rejected (p = 0.21).
  A long, GC-shifted branch pulled toward other fast-evolving outside hits (which were picked *because* they hit it) is the
  likely cause of the 23S placement; the two genes disagree, so by the pre-written rule this is **unresolved**.
- **Vertical inheritance inside the family.** Among genera, the divergent trees share the strongly supported split
  {Aquicultor_A, CADDYS01, DASXVE01 | DASQQM01, Kmv33} with the GTDB species tree (div 23S 98, div 16S 100; GTDB 99).
  No divergent-tree split with UFBoot >= 95 conflicts with the species tree or the normal-copy trees (the only conflicts are
  at 59-92 support). Combined with the fixed locus (T2), this is vertical inheritance, not transfer between genera.
- **Controls behave.** E. coli and Vibrio copies cluster with Gammaproteobacteria; T. bispora's two 16S types (7% apart)
  sit together inside Actinomycetes; the Tenuifilum divergent operon (12%) sits with its own genome's normal copy.
- A third 23S in DASXVE01 (DASXVE010000021.1, 78% to normal) groups with Clostridia: a binned contaminant, not a third type.

**Mechanism.** A single origin (a duplication, or an acquisition too old to trace) at or before the base of
Aquicultorales, then inheritance at a fixed locus through the Aquicultoraceae radiation, with losses in some lineages. It has been kept separate from the normal operon (no gene conversion)
for the whole history of the order, which needs something to keep it.

Files: `src/extract.py`, `src/t1_hits.py`, `src/t1_treeset.py`, `src/t1_tree.sh`, `src/t1_read.py`, `src/t1_au.py`,
`src/t1_vertical.py`; trees `out/t1_23S.contree`, `out/t1_16S.contree`; AU `out/t1_au_*.txt`.
