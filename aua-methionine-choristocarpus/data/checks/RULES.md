# Choristocarpus AUA = Met: independent test on Phaeoexplorer WGS assemblies (written 09:04 14 Nov, before download)
Question: does an independently sequenced and assembled Choristocarpus tenellus mito (Genoscope GCA_964200185.2,
not the Choi et al. record NC_085326) read ATA at Met positions, and does the sister genus Discosporangium
mesarthrocarpum (GCA_964200235.2, Discosporangiaceae) do the same?

Mito contigs: blastn (megablast off, -task dc-megablast, evalue 1e-20) of NC_085326.1 against each assembly; a contig is
mito if aligned length >= 2 kb in total and >= 70% identity. Report contig ids, lengths, and covered share of NC_085326.

Reading (no table assumed for ATA): references = mito CDS translations (as annotated, all read ATA as Ile) of one
genome per brown-algal order other than Sphacelariales from data/records_ranked + data/cds_mito (first by accession),
plus Phaeostrophion ON552993.1 separately reported. tblastn (-db_gencode 1, -seg no, evalue 1e-5) each reference
protein against the mito contigs; per contig codon keep the query residue from each reference's best HSP covering it.
A codon position is M (or I) when >= 2 references align there and >= 2/3 of them agree. Skip contig codons aligned to
reference residues 1-10 (start codons). Count ATA, ATG, ATT, ATC at M and I positions.

Call ATA = Met for an assembly if ATA_M >= 10, ATA_M >= 2 * ATA_I, and ATT_M <= 0.1 * (ATT_M + ATT_I).
Controls run the same way: NC_085326.1 itself (expected M-lean, K1 gave 39/2) and Phaeostrophion ON552993.1 +
Ectocarpales/Laminariales records (expected I). If a control misbehaves, report and do not call.
Reported as "mine only"; then kill subagent before any claim.

Amendment A (09:05, before any read of the downloads): NCBI files Choristocarpus under Sphacelariales, so the reference set is
data/ first-by-accession per order except Sphacelariales: Cutleria NC_085327, Desmarestia NC_007684, Dictyota NC_007685,
Pylaiella NC_003055, Fucus NC_007683, Laminaria NC_004024, Scytothamnus NC_085336, Sporochnus NC_085337.
tblastn uses -db_gencode 4 (TGA=Trp, ATA still Ile) so TGA does not break HSPs. For a control from a reference order the
reference of that order is dropped. Extra negative control: Protohalopteris MZ156064 (Sphacelariales).

Amendment B (09:07, after contig finding, before any codon read): contig 1 of each assembly (cten 120,357 nt; dmes 131,317 nt)
aligns to the Choristocarpus plastid NC_085297 over most of its length (131,232 and 96,674 nt aligned) and is the plastid; it
is dropped. Mito = cten CAXNAD020000071.1 (32,342 nt, 99.9% to NC_085326); dmes CAXNAE020000396.1 (34,744 nt). Dmes
CAXNAE020000402.1 is the same molecule (rotated, 99.99-100% over 34.7 kb) and is dropped so no codon counts twice.
Discosporangium mito aligns to NC_085326 at 76% identity: a different genome, not a copy.

Amendment C (09:13, before any RNA read is looked at): RNA-editing kill test. Choristocarpus RNA-seq ERR12954301 (Phaeoexplorer,
TRANSCRIPTOMIC): the first ~40 MB of read 1 from ENA (gz prefix, streamed). minimap2 -ax sr to the WGS mito contig
CAXNAD020000071.1; primary alignments, MAPQ >= 20, base quality >= 20. At each site in read/cten_WGS_sites.tsv, count read
bases at codon position 3 (strand-aware). Groups: ATA at ref-Met, ATA at ref-Ile, ATG at ref-Met (control).
Editing explains AUA=Met if, among ATA-at-Met sites with depth >= 10, the median share of G (on the coding strand) is >= 0.5.
AUA is read unedited (supports a code change) if that median G share is <= 0.1 and the ATA-at-Ile median is similar.
Fewer than 10 ATA-at-Met sites with depth >= 10: no call, report depth.

Amendment D (09:16): the 40 MB prefix (457,550 reads) gave depth 1-3 at the sites. I saw those few bases before writing this:
the ATA-at-Met sites with reads showed A. Same test and thresholds, on read 1 bytes 0-799,999,999 (~9M reads), streamed
through minimap2; only mapped SAM lines are kept on disk.
Amendment D2 (09:22): the 800 MB stream returned only 79 mapped reads (ENA stopped the stream without an error). I'm refetching
the read-1 prefix as checked 40 MB byte ranges up to 400 MB, then running the same test. No new results were seen.
