DS02 (GCA_002251995.2) mito contigs vs Wickerhamiella mitogenomes, 2026-11-12.
- ds02_vs_bk.paf, ds02_vs_wick.paf: fact-checker's minimap2 -x asm20 -c alignments of the 15 DS02 mito contigs
  (all hit BK074532 W. infanticola only). Identity outside gaps (1-de, weighted by block length): 96.8% over
  295,627 aligned bp; blocks >= 5 kb range 95.3-97.7%; longest block 28,163 bp of alignment.
- wick_gb024.fa: the 20 Wickerhamiella mitogenomes from data/gb_mito/gb_024.gb.gz.
- bk074532_vs_wick.paf: BK074532 vs the other 19 (wick_gb024_noBK074532.fa), minimap2 -t 2 -x asm20 -c.
  Only W. verensis BK074542 aligns: 92.0% outside gaps over 26,444 bp.
