#!/bin/bash
#SBATCH --job-name=count_total_bases
#SBATCH --partition=batch
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --time=12:00:00
#SBATCH --output=count_bases.%j.out
#SBATCH --error=count_bases.%j.err

ml SAMtools/1.23.1-GCC-13.3.0

# Path to main samples directory
# Upload was a directory containing several subdirectories. These subdirectories each contained files that were being uploaded to the SRA.

SAMPLES_DIR="/scratch/kbb45638/telomere_final_analysis/ncbi_upload/upload"

echo -e "Sample_ID\tTotal_Bases"
echo -e "--------------------------"

#
# Loop
#

echo -e "File\tTotal_Bases"

for parent_dir in "$SAMPLES_DIR"/*; do
    if [ -d "$parent_dir" ]; then

        for bam in "$parent_dir"/*.bam; do
            if [ -f "$bam" ]; then

                file_name=$(basename "$bam")

                bases=$(samtools stats -@ 4 "$bam" 2>/dev/null \
                    | grep "total length:" | cut -f3)

                echo -e "${file_name}\t${bases}"

            fi
        done

    fi
