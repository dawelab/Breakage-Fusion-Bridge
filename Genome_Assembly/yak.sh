#!/bin/bash
#SBATCH --job-name=yak
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=64gb
#SBATCH --time=04:00:00
#SBATCH --output=yak.%j.out
#SBATCH --error=yak.%j.err

#Generate parental K-mers
/home/kbb45638/software+environments/yak/yak count -k31 -b37 -t36 -o pat.yak /scratch/kbb45638/c2_tester_assembly/parental_reads/4330-1_R1_001.fastq.gz /scratch/kbb45638/c2_tester_assembly/parental_reads/4330-1_R2_001.fastq.gz

/home/kbb45638/software+environments/yak/yak count -k31 -b37 -t36 -o mat.yak /scratch/kbb45638/c2_tester_assembly/parental_reads/4319-1_R1_001.fastq.gz /scratch/kbb45638/c2_tester_assembly/parental_reads/4319-1_R2_001.fastq.gz
