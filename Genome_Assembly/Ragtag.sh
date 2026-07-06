#!/bin/bash
#SBATCH --job-name=ragtag
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128gb
#SBATCH --time=12:00:00
#SBATCH --output=ragtag.%j.out
#SBATCH --error=ragtag.%j.err

module load RagTag/2.1.0-foss-2022a
ragtag.py scaffold /home/kbb45638/references/GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic.fasta hap2_maternal.fasta -t 12
