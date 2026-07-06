#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=160gb
#SBATCH --time=60:00:00
#SBATCH --output=hifiasm.%j.out
#SBATCH --error=hifiasm.%j.err

ml hifiasm/0.25.0
hifiasm -o trio_integrated_asm -t 36 \
  -1 /scratch/kbb45638/c2_tester_assembly/yak/k31/pat.yak -2 /scratch/kbb45638/c2_tester_assembly/yak/k31/mat.yak \
  --ul /scratch/kbb45638/c2_tester_assembly/Reads_for_Assembly/c2_ONT_15kb.fq.gz \
  /scratch/kbb45638/c2_tester_assembly/HiFi_Reads/L4_Combined.fq
