#!/bin/bash
#SBATCH --job-name=sq
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=02:00:00
#SBATCH --output=sq.%j.out
#SBATCH --error=sq.%j.err

###ONT Reads were first filtered to remove reads under 15kb

ml SeqKit/2.9.0
seqkit seq -j 8 -m 15000 barcode05.fastq.gz -o barcode5_15kb_ONT.fq.gz
seqkit seq -j 8 -m 15000 barcode06.fastq.gz -o barcode6_15kb_ONT.fq.gz
seqkit seq -j 8 -m 15000 c2_control_reads.fasta -o c2_control_15kb_ONT.fq.gz
