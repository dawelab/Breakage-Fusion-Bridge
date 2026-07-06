#!/bin/bash
#SBATCH --job-name=awk
#SBATCH --partition=batch
#SBATCH --mail-type=END
#SBATCH --mail-user=kempton@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --output=awk.%j.out
#SBATCH --error=awk.%j.err

#L or R added to indicate direction of telomere. Samples with R later have their modkit calls put in reverse order for plotting

#4a1 L
awk -v START=178101388 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4a1b1_methylation_NewGenome.bed > 4a1_telo_modkit.bed
awk -v START=178101388 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4a1_ABS_modkit.bed

#4b1 L 
awk -v START=161404502 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4a1b1_methylation_NewGenome.bed > 4b1_telo_modkit.bed
awk -v START=161404502 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4b1_ABS_modkit.bed

#4b2 R
awk -v START=150402744 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4b2_methylation_NewGenome.bed > 4b2_telo_modkit.bed
awk -v START=150402744 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' c2_methylation.bed > 4b2_c2_modkit.bed

#4b3 L Pulling 4b3 from 120-12 becasue it had the highest coverage of 4b3. Extreme resolution would be higher by concatenating all 4b3 samples
awk -v START=151760101 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 120-12_methylation_NewGenome.bed > 4b3_telo_modkit.bed
awk -v START=151760101 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4b3_ABS_modkit.bed

#4b4 L
awk -v START=116023409 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4b4_ONT_methylation_NewGenome.bed > 4b4_telo_modkit.bed
awk -v START=116023409 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4b4_ABS_modkit.bed

#4b5 L
awk -v START=130346197 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4b5_methylation_NewGenome.bed > 4b5_telo_modkit.bed
awk -v START=130346197 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4b5_ABS_modkit.bed

#4b7 L
awk -v START=172716549 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4b7_methylation_NewGenome.bed > 4b7_telo_modkit.bed
awk -v START=172716549 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4b7_ABS_modkit.bed

#4b8 L
awk -v START=127343146 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4b8_methylation_NewGenome.bed > 4b8_telo_modkit.bed
awk -v START=127343146 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4b8_ABS_modkit.bed

#4a9 R
awk -v START=163361467 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4a9b9_methylation_NewGenome.bed > 4a9_telo_modkit.bed
awk -v START=163361467 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4a9_ABS_modkit.bed

#4b9 R
awk -v START=163414315 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4a9b9_methylation_NewGenome.bed > 4b9_telo_modkit.bed
awk -v START=163414315 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' c2_methylation.bed > 4b9_c2_modkit.bed

#4b10 R
awk -v START=125923390 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4b10_methylation_NewGenome.bed > 4b10_telo_modkit.bed
awk -v START=125923390 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4b10_ABS_modkit.bed

#4a3b3 homo 4b3.d1 L
awk -v START=168797655 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4a3b3_homo_methylation_NewGenome.bed > 4b3.d1_telo_modkit.bed
awk -v START=168797655 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' c2_methylation.bed > 4b3.d1_c2_modkit.bed

#4a3b3 homo 4a3 R
awk -v START=186519890 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4a3b3_homo_methylation_NewGenome.bed > 4a3_telo_modkit.bed
awk -v START=186519890 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4a3_ABS_modkit.bed

#4a3b3 sup 4a3.d1 L
awk -v START=171651351 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4a3b3_super_methylation_NewGenome.bed > 4a3.d1_telo_modkit.bed
awk -v START=171651351 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' c2_methylation.bed > 4a3.d1_c2_modkit.bed

#121-5 4a3.d2 L
awk -v START=141286410 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 121-5_methylation_NewGenome.bed > 4a3.d2_telo_modkit.bed
awk -v START=141286410 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' c2_methylation.bed > 4a3.d2_c2_modkit.bed

#120-12 4a3.d3 Staggered Edge. Ignore
#awk -v START=141286410 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 120-12_methylation_NewGenome.bed > 4b3.d3_telo_modkit.bed

#4307x4309-1 C2-10 4a3.d4 R
awk -v START=114927974 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4309_1_C2_10_methylation_NewGenome.bed > 4a3.d4_telo_modkit.bed
awk -v START=114927974 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' c2_methylation.bed > 4a3.d4_c2_modkit.bed

#4307x4309-3 C2-10 4a3.d5 R
awk -v START=149905150 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4309_3_C2_10_methylation_NewGenome.bed > 4a3.d5_telo_modkit.bed
awk -v START=149905150 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' c2_methylation.bed > 4a3.d5_c2_modkit.bed

#4307x4309-3 C2-3 4a3.d6 L
awk -v START=162247299 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 4309_3_C2_3_methylation_NewGenome.bed > 4a3.d6_telo_modkit.bed
awk -v START=162247299 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' c2_methylation.bed > 4a3.d6_c2_modkit.bed

#112-5 4a8.d1 R
awk -v START=151397899 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 112-5_methylation_NewGenome.bed > 4a8.d1_telo_modkit.bed
awk -v START=151397899 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' c2_methylation.bed > 4a8.d1_c2_modkit.bed

#115-2 4a8.d2 R
awk -v START=131347204 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 115-2_methylation_NewGenome.bed > 4a8.d2_telo_modkit.bed
awk -v START=131347204 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 4a8.d2_ABS_modkit.bed

#128-1 174,660 ABS R
awk -v START=174653704 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 128-1_methylation_NewGenome.bed > 128-1_telo_modkit.bed
awk -v START=174653704 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 128-1_ABS_modkit.bed

#128-1 166,995 c2 Ignore because only 1 read and still unconfirmed

#137-4 L
awk -v START=186426084 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' 137-4_methylation_NewGenome.bed > 137-4_telo_modkit.bed
awk -v START=186426084 '$4!="h" && ($1=="chr4_RagTag" || $1=="CM039153.1_RagTag") && $2>=START && $2<=(START+10000)' ABS_methylation.bed > 137-4_ABS_modkit.bed