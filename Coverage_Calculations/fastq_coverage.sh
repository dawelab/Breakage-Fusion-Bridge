#!/bin/bash
#SBATCH --job-name=count_bases
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=count_%j.out
#SBATCH --error=count_%j.err


echo "=== Processing file (can add file name here for clean output) ==="
date
pigz -p 4 -dc ###fastq.gz of fq.gz### | awk 'NR%4==2 {total += length($0)} END {printf "Total bases: %'\''d\n", total}'
