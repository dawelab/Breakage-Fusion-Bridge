#!/bin/bash
#SBATCH --job-name=bamcoverage
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5gb
#SBATCH --time=00:10:00
#SBATCH --output=bamcoverage.%j.out
#SBATCH --error=bamcoverage.%j.err

replace 4b# with neochromosome of interest. chr4_reads_###.bam is one of the output files from the BFB_pipeline scripts. 

ml deepTools/3.5.5-gfbf-2023a
bamCoverage -b chr4_reads_4b#.bam -o 4b#_coverage_100.bw --binSize 100 --normalizeUsing None -p max
