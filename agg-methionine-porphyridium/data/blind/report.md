# Blinded AGG panel R01–R09: unblinded 2026-11-11 15:22 EST

- **Agent:** a fresh subagent (a97848f0bfd348ff2), told only to infer each genome's mitochondrial code.
- **Methods:** Codetta v2, plus its own reference-protein likelihood using Swiss-Prot mitochondrial proteins.
- **Calls:** `calls.md`. **Scripts:** `agent/scripts/`.
- **Rule:** `../blind_prereg.md`, fixed 15:00 before launch.

## Verdict: PASS

| Code | Accession | Organism | Role | Agent's final call on AGG | Score |
|---|---|---|---|---|---|
| R06 | NC_062302.1 | *Porphyridium aerugineum* SAG 1380-2 | test | **Met**. Method 2: n 52, LLR vs Arg 73.5. Codetta M: 26 columns, gap 85.7 vs Arg | hit |
| R02 | MT483995.1 | *P. aerugineum* CCMP1948 | test | **Met**. Method 2: n 50, LLR 72.3. Codetta M: 25 columns, gap 68.8 | hit |
| R05 | MT483996.1 | *P. purpureum* | negative | not Met (table 1 + TGA = Trp) | ok |
| R03 | NC_062401.1 | *Erythrolobus coxiae* | negative | not Met | ok |
| R07 | NC_062303.1 | *Rhodella violacea* | negative | not Met | ok |
| R08 | NC_062397.1 | *Neorhodella cyanea* | negative | not Met | ok |
| R09 | NC_062400.1 | *Dixoniella grisea* | negative | not Met | ok |
| R01 | NC_024758.1 | *Chromochloris zofingiensis* | positive control | **Met**. Method 2: n 105, LLR 75.3. Codetta M: 48 columns | sensitive |
| R04 | CM026549.1 | *Scenedesmus* sp. | AGG = Ala control | **Ala**. Method 2: n 92. Codetta A: 40 columns | ok |

**No tier ambiguity.** Every AGG call is a "call" under both the agent's frozen v2 rule and its post-hoc v3 filter.
The v3 filter can only remove calls; it removed one, R05 TAA → Met.

## Other calls (recorded; they don't change the verdict)

- **R04, *Scenedesmus*.** TAG = Leu is called. TCA is flagged as a "weak, possible stop": it is the first stop-like
  codon after 10 of 12 gene ends. Together these are the known Scenedesmus code (NCBI table 22), recovered blind.
- **R01, *Chromochloris*.** CGG = Leu is called by Method 2 alone; Codetta's top residue is also Leu but below its
  threshold.
- **TGA = Trp** in the red algae, as expected from their code.
