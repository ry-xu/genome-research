# Blinded AGG panel: pass rule (fixed 2026-11-11 15:00 EST, before launch)

Panel: 9 mitogenomes coded R01-R09 (key: blind_key.tsv, not shown to the agent).
- Tests: NC_062302.1 and MT483995.1 (*Porphyridium aerugineum*, two isolates).
- Negatives (red algae, same family/order, AGG expected Arg): MT483996.1, NC_062400.1, NC_062401.1, NC_062303.1, NC_062397.1.
- Positive control, known AGG = Met: NC_024758.1 (*Chromochloris zofingiensis*, Noutahi et al. 2019).
- Different-meaning control, AGG = Ala expected: CM026549.1 (*Scenedesmus*).

PASS: the agent's final call is AGG = Met in both tests, and in none of the 5 negatives (Arg or "no call" is fine there).
FAIL: either test not called AGG = Met (Arg, another amino acid or no call), or any negative called AGG = Met.
Sensitivity: if the positive control is not called AGG = Met, a missed test counts as "uninformative", not as evidence against.
Calls on other codons are recorded but do not change PASS/FAIL. "Called" means the agent's own final, frozen call; its tiers are read
as it defines them, and I report both readings if a tier is ambiguous.
