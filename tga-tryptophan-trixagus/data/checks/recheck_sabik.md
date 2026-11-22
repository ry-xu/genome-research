# trixagus-tga-trp: second check of the parts Mintaka did not recount (Sabik, 21:36 20 Nov). Own code on NCBI efetch GenBank files.
- Pseudogenes (biopython on CP157893.2 Tser, CP157895.2 Tder, CP158689.2 Ttur, CP157896.2 Tcar): CDS with protein 174/158/153/166;
  pseudo CDS 27/44/24/0; with PGAP "internal stop" note 21/25/19/0 = post's "19-25". Share 21/195, 25/183, 19/172 = 10.8/13.7/11.0%
  = post's "11-14%" (denominator CDS + internal-stop pseudogenes; 10.8 rounds to 11). My own frame-intact pseudo CDS whose only
  in-frame stops are TGA: 21/23/19/0. All records /transl_table=11.
- RF2: every PGAP release-factor CDS (inference HMM:NF015427). Readers have exactly one each, "PCRF domain-containing protein",
  GGQ at 216/218/207, carrying RF1's PxT motif (PKT/PIT/PKT) ~40-75 aa before GGQ: RF1 kept, no RF2. Tcar has two: RF1
  (ABNO60_00780, PKT) and "peptide chain release factor-like protein" ABNO60_00745 (345 aa, GGQ@231, no SPF and no PxT in the
  window) = the post's divergent RF2 that lost SPF. Consistent with the post. Name searches for "release factor 2"/prfB find
  nothing in any record (PGAP names neither), so the RF2 call needs a profile or motif, as the project did.
- Novelty: one web search (Shikimatogenerans Throscidae Trixagus Wierz Engl 2026) finds no published throscid paper yet.
- Not checked: raw-read RF2 scan, tree/AU tests, Codetta re-run, Nardonella.
Verdict: checked numbers reproduce; the post's own advice (write to the Jena group first) stands.
