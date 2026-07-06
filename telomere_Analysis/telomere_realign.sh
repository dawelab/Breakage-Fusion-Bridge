#!/bin/bash
#SBATCH --job-name=maize_bfb_pipeline
#SBATCH --partition=batch
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40gb
#SBATCH --time=4:00:00
#SBATCH --output=telo_realign.%j.out
#SBATCH --error=telo_realign.%j.err

ml SAMtools/1.21-GCC-13.3.0
ml modkit/0.5.0-GCCcore-13.3.0
# ==========================================
# CONFIGURATION
# ==========================================
DIR_LIST="directory_list.txt"                 # Text file with directory names for samples that have a telomere addition site
PARENT_DIR="/scratch/kbb45638/telomere_final_analysis/analysis"        # Subdirectory where all of the directories in DIR_LIST are located
NEW_REF="/scratch/kbb45638/telomere_final_analysis/reference/ABSxc2v2.1_masked.fasta"     # The new reference genome
THREADS=12
OUTPUT_DIR=$(pwd)

# ==========================================
# PIPELINE EXECUTION
# ==========================================

# Read the text file line by line
while IFS= read -r DIR_NAME; do
    
    # Skip empty lines in the text file just in case
    [ -z "$DIR_NAME" ] && continue 
    
    TARGET_DIR="$PARENT_DIR/$DIR_NAME"
    
    # 1. Sanity Check: Does the directory exist?
    if [ ! -d "$TARGET_DIR" ]; then
        echo "WARNING: Directory not found -> $TARGET_DIR. Skipping..."
        continue
    fi

    # 2. Find the input BAM using an array to safely handle the wildcard
    BAM_FILES=("$TARGET_DIR"/*_telo_sorted.bam)
    INPUT_BAM="${BAM_FILES[0]}"

    # Sanity Check: Does the file actually exist?
    if [ ! -e "$INPUT_BAM" ]; then
        echo "WARNING: No *_telo_sorted.bam found in $TARGET_DIR. Skipping..."
        continue
    fi

    # Extract the prefix to dynamically name the output files
    FILENAME=$(basename "$INPUT_BAM")
    PREFIX="${FILENAME%_telo_sorted.bam}"
    
    # Define the output file names
    NEW_BAM="$OUTPUT_DIR/${PREFIX}_realigned_NewGenome.bam"
    BEDMETHYL="$OUTPUT_DIR/${PREFIX}_methylation_NewGenome.bed"
    
    echo "--------------------------------------------------------"
    echo "Processing: $PREFIX in $DIR_NAME"
    
    # 3. Re-alignment & Sorting
    echo "   -> Aligning to new reference and sorting..."
    singularity exec /apps/singularity-images/dorado_1.1.1.sif dorado aligner "$NEW_REF" "$INPUT_BAM" -t "$THREADS" | samtools sort -@ "$THREADS" -o "$NEW_BAM"
    
    # 4. Indexing (Required before Modkit can run)
    echo "   -> Indexing new BAM..."
    samtools index -@ "$THREADS" "$NEW_BAM"
    
    # 5. Modkit bedMethyl Extraction
    echo "   -> Generating bedMethyl table for Cytosines..."
    modkit pileup --ref "$NEW_REF" --threads "$THREADS" "$NEW_BAM" "$BEDMETHYL"
    
    echo "Finished: $PREFIX"
    
done < "$DIR_LIST"

echo "--------------------------------------------------------"
echo "All directories processed successfully!"
