# **Molecular consequences of chromosome and chromatid type breakage fusion bridge cycles in maize**

This repository contains all the code used in "Molecular consequences of chromosome and chromatid type breakage fusion bridge cycles in maize". We utilized an inducible centromere system in Maize to generate 11 chromosome BFB derived neochromosomes and 11 chromatid BFB derived neochromosomes. Nanopore long-read sequencing was used to study the structural variants and neotelomeres formed. Chromatid BFB structural variants were limited to foldback duplications (also known as foldback inversions and inverted duplications). Chromosome type BFB was found to produce tandem duplications and internal deletions as well as complex rearrangements resembling chromothripsis. Neotelomere formation sites largely mirrored neotelomere sites in humans and wheat where formation typically occurs over regions with homology to the telomerase RNA template. These addition sites occasionally contained imperfect telomere motifs at addition sites. Several addition sites showed an unreported pattern where new sequence was inserted between the telomere motif and the DNA break site. The inserted sequence was found to contain reverse homology just upstream of the addition site.

All raw reads used in this analysis can be found under Bioproject Accession Number PRJNA1493139

# Directories

-   [Genome Assembly](#genome-assembly)
-   [Dorado Scripts](#dorado-scripts)
-   [BFB Pipelines](#bfb-pipelines)
-   [SV Plots](#sv-plots)
-   [Telomere Analysis](#telomere-analysis)
-   [Breaksite Analysis](#breaksite-analysis)
-   [Coverage Calculations](#coverage-calculations)

## Genome Assembly

#### This directory section contains all the code needed to generate the c2 assembly used in the paper.

1.  yak.sh: Count K-mers from the illumina reads for triobinning
2.  15kb_Filter.sh: Filter Nanopore reads to remove reads under 15kb
3.  ONT_cat.sh: Combine the three filtered fastq.gz files into a single fastq.gz
4.  Hifiasm.sh: Assembly step. HiFi reads from 4b(4), ONT reads from step 3, and the .yak files from step 1 are used to generate a triobinned assembly. .gfa files from each parent are converted to fasta.
5.  Ragtag.sh: Scaffold the contigs from Hifiasm using Mo17 genome.

###### Note: If attempting this in the future and a new W22 genome is released scaffolding on that will likely work better because c2 is derived from W22. When this work was done, the available W22 reference was from 2017 and contained 68,134 contigs. Doing this will change coordinates of breakpoints later!

6.  Command_Line_Final.sh: Commands run directly on terminal to combined the new reference and ABS to create a diploid reference and remove unplaced contigs.\
7.  mask_telo.py: This hard masks (replace with N) all telomere repeats. This was used on the diploid refernece for the telomere analysis.

## Dorado Scripts

#### This directory contains all the code to process the raw pod5 and multiplexed files prior to analysis. The outputs from this are on the SRA.

###### Note: These scripts use dorado version 1.1.1 in a singularity container. The only reason a singularity container was used is because at the time this analysis was done Dorado 1.1.1 was the newest version on our HPC and it was only installed in a singularity container.

1.  dorado_basecall.sh: Basecalls .pod5 files.
2.  dorado_demux.sh: Demultiplexes the libraries generated with SQK-NBD114.24

## BFB Pipelines

#### This directory contains the scripts used to do alignment, initial neotelomere detection, and variant calling. There are three scripts here, but all do nearly the same thing. The only difference is in Step 1. This is the most computationally intense step outside of basecalling and should be done an an HPC. 

Running these scripts in a directory with the name for your sample will result in the easiest analysis. The initial command uses a \* wildcard to just grab the unaligned .bam file (or .fq), so be sure there are not other files in the directory.

### Step 1:

This step does the initial genome alignment, sorting, and indexing. `Dorado aligner` is used for the standard pipeline and the multi input pipeline. Minimap2 is used for for Hifi reads. For the multi input pipeline there are several lines for different dorado alignments. These are then combined with `samtools merge`. Lines can be added or removed for samples with more or less than 3 input bam files. After this all the scripts are identical. Reads are then sorted and indexed.

### Step 2:

This step extracts reads with a minimum of 4 intact telomere 7-mers for initial neotelomere identification. These alignments are not used for the later telomere analysis.

### Step 3:

This step creates filtered alignments that remove multimapping reads and limit alignment to chromosome 4. A SM tag is also added here to ensure Freebayes can handle the files in the next step. This file `chr4_reads\_\${PREFIX}.bam` is used for all the subsequent analyses. Adaptive sampling primarily enriched chromosome 4 and all SVs should be on this chromosome.

### Step 4:

This step calls SNPs on chromosome 4 from both genomes for step 5.

### Step 5:

This uses the SNPs called in step 4 and the chromosome 4 alignment to remove reads that contain called SNPs. This helps to remove residual background from small unassembled or unscaffolded parts of the genome. The SNP filtered alignment is very similar, but these background reads can create a lot of false structural variation.

### Step 6:

This creates a new alignment that has foldback candidates. Nanopore R10.4 chemistry has a frequent false foldback signature. This has been speculated to arise from when both strands of a read are read in quick succession and not split by the basecaller (failed duplex reads). This unfortunately looks like a foldback inversion. This script gets rid of these false signatures by finding reads with a primary and supplementary alignment on the same chromosome and in the same location but on opposite strands. It then checks if one alignment is at least 1.5 times longer than the other. This gets rid of a lot of the false alignments because the positive and negative alignment frequently have similar lengths. Manual verification is still essential because a lot of false foldbacks still make it through.

### Step 7:

This step runs Sniffles and Severus on the chromosome 4 and SNP filtered chromosome 4 alignments. Manual variant calling with by coverage analysis and split read identification in IGV is still incredibly useful. Calling large and nested SVs in a heterogeneous and repetative genome like Maize is quite challenging!

## SV Plots

#### This directory contains the code and input files used to create the coverage plots in figure 2. 

pyGenomeTracks was installed in a conda environment and run on command line to generate each plot where 4b# is replaced with 4, 5, or 10 to generate each plot.

`pyGenomeTracks --tracks tracks_4b#.ini --region chr4_RagTag:110000000-190000000 --out BFB_variant_landscape_region_4b#.pdf`

## Telomere Analysis

#### This directory contains the code used for the neotelomere analysis. 

telomere_realign.sh: Realign all neotelomere reads to telomere masked reference genomes and runs `modkit pileup`. Realigning to a telomere masked reference gives more coverage where methylation state is callable. Reads where the primary alignment is at a canonical telomere and the supplementary alignment are at the neotelomere site are not counted by `modkit pileup`. Realigning to the telomere masked reference makes the neotelomere site the primary alignment.

site_extraction.sh: Extracts neotelomere sites from .bedmethyl and removes 5hm calls.

Metaplot_position_based.py: Plots neotelomere methylation. This script requires .bedmethyl to be in a specific directory structure. Two initial sub directories called "Telo_Left" and "Telo_Right". Each of these directories contains directories with the sample name. The sample name directories should contain .bedmethyl files for the control and neotelomere samples.

## Breaksite Analysis

1.  w22_c2_blast_loop.sh: Finds c2 breaksite locations in the W22 genome. The regions 10 kb upstream and downstream the breaksite is extracted and blasted to the W22 genome. The original position is then converted to a coordinate in the W22 genome
2.  breaksite_enrichment_ABS.py and w22_feature_annotate.py: Creates a file that contains information on breaksite coordinates and genomic feature. The outputs from these files are also included in the Breaksite_Analysis directory because these scripts require the .gff3 files for W22 and ABS.
3.  piechart_combined.py and piechart_combined_ABS.py: These create the piechart seen in Figure 4 A. The only difference is that piechart_combined.py does not include ABS as a category.
4.  combined_ABS_categoy.py: This creates figure 4B and 4C. These figures show the distribution of neoteolomeres and their orientation and plot the location of breaks along with genomic context.

## Coverage Calculations

These scripts were run to calculate coverage. the .sh scripts (bam_cover_loop.sh and fastq_coverage.sh) were used to calculate coverage for table S3. 4b4_Coverage were run from command line to calculate coverage of the 4b(4) neochromosome to compare the ONT and HiFi data.
