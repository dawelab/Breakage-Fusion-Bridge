#!/bin/bash
#SBATCH --job-name=dorado
#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:A100:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50gb
#SBATCH --time=72:00:00
#SBATCH --output=dorado.%j.out
#SBATCH --error=dorado.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=### Email ###

### Replace Sample_Name with the directory name containing the Pod5 files. Keep this SAMPLE variable consistent between downstream scripts.
SAMPLE="Sample_Name"
singularity exec /apps/singularity-images/dorado_1.1.1.sif dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.2.0 --kit-name SQK-NBD114-24 ### Path to Pod5 Directory ###/${SAMPLE} --modified-bases-models dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mC_5hmC@v2 > ${SAMPLE}_unaligned.bam
