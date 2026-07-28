# **Molecular consequences of chromosome and chromatid type breakage fusion bridge cycles in maize**

This repository contains all the code used in "Molecular consequences of chromosome and chromatid type breakage fusion bridge cycles in maize". We utilized an inducible centromere system in Maize to generate 11 chromosome BFB derived neochromosomes and 11 chromatid BFB derived neochromosomes. Nanopore long-read sequencing was used to study the structural variants and neotelomeres formed. Chromatid BFB strucutral variants were limited to foldback duplications (also known as foldback inversions and inverted duplications). Chromosome type BFB was found to produce tandem duplications and internal deletions as well as complex rearrangements resembling chromothripsis. Neotelomere formation sites largely mirrored neoteolomere sites in humans and wheat where formation typically occurs over regions with homology to the telomerase RNA template. These addition sites occasionally contained imperfect telomere motifs at addition sites. Several addition sites showed an unreported pattern where new sequence was inserted between the telomere motif and the DNA break site. The inserted sequence was found to contain reverse homology just upstream of the addition site.

# Repositories

-   [Genome Assembly](#genome-assembly)
-   [Dorado Scripts](#dorado-scripts)
-   [BFB Pipelines](#bfb-pipelines)
-   [SV Plots](#sv-plots)
-   [Telomere Analysis](#telomere-analysis)
-   [Breaksite Analysis](#breaksite-analysis)
-   [Coverage Calculations](#coverage-calculations)

## Genome Assembly

#### This section contains all the code needed to generate the c2 assembly used in the paper.

1)  yak.sh: Count K-mers from the illumina reads for triobinning
2)  15kb_Filter.sh: Filter Nanopore reads to remove reads under 15kb
3)  ONT_cat.sh: Combine the three filtered fastq.gz files into a single fastq.gz
4)  Hifiasm.sh: Assembly step. HiFi reads from 4b(4), ONT reads from step 3, and the .yak files from step 1 are used to generate a triobinned assembly. .gfa files from each parent are converted to fasta.
5)  Ragtag.sh: Scaffold the contigs from Hifiasm using Mo17 genome.

###### Note: If attempting this in the future and a new W22 genome is released scaffolding on that will likely work better because c2 is derived from W22. When this work was done, the available W22 reference was from 2017 and contained 68,134 contigs. Doing this will change coordinates of breakpoints later!

6)  Command_Line_Final.sh: Commands run directly on terminal to combined the new reference and ABS to create a diploid reference and remove unplaced contigs.\
7)  mask_telo.py: This hard masks (replace with N) all telomere repeats. This was used on the diploid refernce for the telomere analysis.\

