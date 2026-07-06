#!/bin/bash
#SBATCH --job-name=sq
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=02:00:00
#SBATCH --output=sq.%j.out
#SBATCH --error=sq.%j.err

### the c2 gneome was sequenced 3 seperate times using ONT so these sequencing runs were all combined for assembly

cat barcode5_15kb_ONT.fq.gz barcode6_15kb_ONT.fq.gz c2_control_15kb_ONT.fq.gz > c2_ONT_15kb.fq.gz
