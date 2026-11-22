## Draft 6. Trixagus symbionts, TGA = Trp and RF2 loss (`outbox/2026-11-11-trixagus-tga-trp/post.md`)

**Verdict: STANDS.** I found no report of a code change in Shikimatogenerans or in *Nardonella*. The throscid
symbiont paper named in the GenBank records is still not published or preprinted anywhere I could search, and
the records themselves use table 11. The draft's account of Choi et al. 2026, Parks et al. 2026 and gTranslate
2026 is accurate.

One sentence I could not check: "a 2025 conference abstract on the work does not mention the code." I did not
find that abstract.

Optional addition: the *Nardonella* genome has a paper (Palmieri, Pavarini & Sharma 2022, *Microbiol Resour
Announc*), which the draft does not cite. It says nothing about the code.

### Verified

- **GenBank CP157895.2, CP158689.2, CP157896.2** (efetch): "MAG: Candidatus Shikimatogenerans sp. Tder
  chromosome" and the Ttur and Tcar equivalents. Reference: "Wierz,J.C., Malm,H., Kaltenpoth,M. and Engl,T.",
  "Diversity, functionality, and evolutionary history of bacterial symbionts in false click beetles (Coleoptera,
  Throscidae)", "Unpublished"; submitted 04-JUN-2024; BioProject PRJNA1062330. Every CDS is `transl_table=11`.
  The PGAP comment gives "Pseudo Genes (internal stop) :: 25 of 44" for Tder and "19 of 24" for Ttur. Meaning:
  the depositors annotated with the standard bacterial code, and PGAP flagged the internal stops as pseudogenes,
  as the draft says.
- **Choi et al. 2026, *ISME J*** (PMC13712303;
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13712303/fullTextXML): "To date, only six endosymbionts,
  including this first case in Bacteroidota identified here, have been found to use the genetic code 4 (Table
  1)."; "Except for Brownia EUSC, all the other Bacteroidota symbionts retain both prfA and prfB in their
  genomes". It names Shikimatogenerans only as a bostrichid beetle symbiont. Meaning: a 2026 paper that
  surveyed exactly this question did not know of the throscid cases. If the *Trixagus* result holds, Choi's
  "only six" and "all the other Bacteroidota symbionts retain ... prfB" are both out of date, which the draft
  could say.
- **Parks et al. 2026, *Microb Genom*** (PMC13462651;
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13462651/fullTextXML): "suggesting that at least two
  independent reassignment events have occurred in the MGR lineage"; "Notably, many species that use GC4 fail
  the GTDB QC criteria as they lack a sufficient number of phylogenetic marker genes due to their highly reduced
  genomes." Its list of known UGA = Trp lineages has no Shikimatogenerans. Meaning: matches the draft; the
  second quote explains why GTDB-based screens would miss these tiny genomes.
- **Chaumeil et al. 2026, *PeerJ*** (gTranslate; PMC13695967;
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13695967/fullTextXML): "we identify a basal lineage of Ca.
  Stammera capleta that uses the standard bacterial genetic code"; "the first instances of UGA-to-tryptophan
  reassignment in the Patescibacteriota class Minisyncoccia". No Shikimatogenerans. Matches the draft.
- **Kiefer et al. 2023, *ISME J*** (PMC10284843) and **Kiefer et al. 2021, *Commun Biol*** (PMC8113238), both
  Europe PMC full text: no match for UGA, TGA, RF2, prfB, "release factor", "genetic code" or "translation
  table". Meaning: the papers that named and surveyed Shikimatogenerans report no code change.
- **Wierz et al. 2024, *ISME J*** (PMC11322605; Europe PMC full text): "under SRP482428 as part of the BioProject
  PRJNA1062330 for the Throscidae samples". No code mention. Meaning: the throscid reads were deposited with
  this earlier paper; the genome analysis is the unpublished one.
- **Lukeš et al. 2025, Table S1**: bacterial UGA = Trp is listed for "Mycoplasmatales and several separate
  lineages of bacterial endosymbionts"; no lineage of Shikimatogenerans or *Nardonella* is named.
- **Palmieri, Pavarini & Sharma 2022, *Microbiol Resour Announc*** (PMC9670985;
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9670985/fullTextXML): "The final genome was annotated
  using PGAP version 6.0"; "Although some genes are incomplete". No code mention. Assembly record
  GCA_023651715.1 (NCBI Datasets): University of Wisconsin-Madison, 2022-06-06, PRJNA798699.

### Listing only (not opened)

- Krüsemer et al. 2026, *Nat Commun* (PMC13043703): opened only far enough to see that it covers Shikimatogenerans
  in the beetle *Oryzaephilus*, not throscids.
- A bioRxiv preprint titled "A mitochondrion-like ribosome in an insect endosymbiont" (abstract only; it is
  about *Hodgkinia*).

### Searches run

Europe PMC:

- `Shikimatogenerans`
- `Throscidae AND (symbiont* OR endosymbiont* OR Bacteroidota)`
- `Trixagus`
- `(Bostrichidae OR Throscidae OR Silvanidae OR Nosodendridae) AND (symbiont OR endosymbiont) AND ("UGA" OR "TGA" OR "genetic code" OR "translation table" OR prfB OR "release factor 2")`
- `"false click beetles"`
- `AUTH:"Wierz JC" AND PUB_YEAR:[2024 TO 2026]`
- `AUTH:"Engl T" AND (symbiont OR endosymbiont) AND PUB_YEAR:[2025 TO 2026]`
- `(Throscidae OR "throscid" OR Eucnemidae OR Elateroidea) AND (symbio* OR Bacteroidota OR Flavobacteri*) AND PUB_YEAR:[2023 TO 2026]`
- `Nardonella AND Metamasius`
- `Nardonella AND ("UGA" OR "TGA" OR "genetic code" OR "translation table" OR prfB OR "release factor")`

Web searches:

- `"Diversity, functionality, and evolutionary history of bacterial symbionts in false click beetles" Throscidae Wierz Engl Kaltenpoth`
- `Throscidae Shikimatogenerans symbiont genome "genetic code" OR "UGA" OR "stop codon" OR "translation table 4" Trixagus endosymbiont tyrosine Bacteroidota beetle`
- `Wierz dissertation OR thesis Throscidae symbiont Shikimatogenerans "Trixagus" genome reduction tyrosine Kaltenpoth Engl 2025`

"Shikimatogenerans" returned 16 records in Europe PMC and "Trixagus" 36; none reports a code change, and none
of the *Trixagus* records is about symbionts. The exact title of the unpublished paper returned nothing.

### Could not check

- The 2025 conference abstract the draft cites.
- A Wierz thesis, if one exists (the Mainz repository surfaced in a search; I did not find a throscid chapter).
- The Edmond (Max Planck) data repository.
- Whether the Jena group has a manuscript under review. Only they can say; the draft already recommends
  writing to them first.

---

