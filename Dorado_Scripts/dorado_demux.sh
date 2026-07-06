#!/bin/bash
#SBATCH --job-name=dorado
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10gb
#SBATCH --time=12:00:00
#SBATCH --output=dorado.%j.out
#SBATCH --error=dorado.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=### Email ###


SAMPLE="Sample_Name"
singularity exec --nv /apps/singularity-images/dorado_1.1.1.sif dorado demux -o Demux_${SAMPLE} --no-classify -t 10 --sample-sheet #Path to sample_sheet/${SAMPLE}.sample_sheet ${SAMPLE}_unaligned.bam
