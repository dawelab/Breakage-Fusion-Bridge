#!/bin/bash
#SBATCH --job-name=maize_bfb_pipeline
#SBATCH --partition=batch
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50gb
#SBATCH --time=14:00:00
#SBATCH --output=BFB_pipeline.%j.out
#SBATCH --error=BFB_pipeline.%j.err

#Run this script in a directory ONLY the unaligned .bam from dorado_demux.sh. Name this directory the prefix that will signify the sample. 

PREFIX=$(basename "$PWD")
REFERENCE="Path to Reference"
POS_TELO="TTTAGGGTTTAGGGTTTAGGGTTTAGGG"
NEG_TELO="CCCTAAACCCTAAACCCTAAACCCTAAA"

echo "Starting pipeline for ${PREFIX}..."

# ---------------------------------------------------------
# 1. Global Alignment with Dorado
# ---------------------------------------------------------
ml SAMtools/1.21-GCC-13.3.0
echo "Aligning reads..."
singularity exec /apps/singularity-images/dorado_1.1.1.sif dorado aligner -t 20 $REFERENCE *.bam > ${PREFIX}_aligned.bam

echo "Sorting and indexing global BAM..."
samtools sort -@ 20 -o ${PREFIX}_sorted.bam ${PREFIX}_aligned.bam
samtools index ${PREFIX}_sorted.bam

# ---------------------------------------------------------
# 2. Extract Neotelomere Reads
# ---------------------------------------------------------
echo "Extracting neotelomere reads..."
samtools view -h ${PREFIX}_sorted.bam \
  | awk -v p1="$POS_TELO" -v p2="$NEG_TELO" '
      $1 ~ /^@/ || $10 ~ p1 || $10 ~ p2 {print}
    ' \
  | samtools view -b -o ${PREFIX}_telo_reads.bam
samtools sort -@ 20 -o ${PREFIX}_telo_sorted.bam ${PREFIX}_telo_reads.bam
samtools index ${PREFIX}_telo_sorted.bam

# ---------------------------------------------------------
# 3. Extract Chr4, Filter MapQ, & Fix Read Groups
# ---------------------------------------------------------
#echo "Extracting Chromosome 4 reads and filtering MapQ 0..."
# Output to a temporary file first
samtools view -b -q 1 ${PREFIX}_sorted.bam chr4_RagTag CM039153.1_RagTag > tmp_chr4_${PREFIX}.bam

echo "Adding missing SM (Sample) tag for Freebayes..."
# Assign a clean read group with the SM tag to all reads
samtools addreplacerg -r "@RG\tID:1\tSM:${PREFIX}" -o chr4_reads_${PREFIX}.bam tmp_chr4_${PREFIX}.bam

# Clean up the temporary file and index the new, fixed BAM
rm tmp_chr4_${PREFIX}.bam
samtools index chr4_reads_${PREFIX}.bam

# ---------------------------------------------------------
# 4. Variant Calling (Freebayes)
# ---------------------------------------------------------
echo "Calling SNPs on Chr4..."
module load freebayes/1.3.7-gfbf-2024a-R-4.4.2
freebayes -f $REFERENCE chr4_reads_${PREFIX}.bam > ${PREFIX}_SNPs.vcf

# ---------------------------------------------------------
# 5. Read Filtering (Remove Alternate Alleles)
# ---------------------------------------------------------
echo "Filtering reads with alternate alleles..."
ml Python
ml Pysam
python <<EOF
import pysam

prefix = "${PREFIX}"
bam_file = f"chr4_reads_{prefix}.bam"
vcf_file = f"{prefix}_SNPs.vcf"
output_bam_file = f"{prefix}_SNP_filtered.bam"

def load_snps(vcf_file):
    snp_positions = {}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            if len(record.alts) == 1:
                snp_positions[(record.chrom, record.pos)] = (record.ref, record.alts[0])
    return snp_positions

def contains_snp(read, snp_positions):
    if read.is_unmapped or not read.query_sequence:
        return False
    
    # get_aligned_pairs handles CIGAR parsing automatically (skips indels/clipping)
    aligned_pairs = read.get_aligned_pairs(matches_only=True)
    read_seq = read.query_sequence
    ref_name = read.reference_name
    
    for q_pos, r_pos in aligned_pairs:
        chrom_pos = (ref_name, r_pos + 1) # pysam is 0-based, VCF is 1-based
        if chrom_pos in snp_positions:
            ref_base, alt_base = snp_positions[chrom_pos]
            if read_seq[q_pos] == alt_base:
                return True
    return False

snp_positions = load_snps(vcf_file)
with pysam.AlignmentFile(bam_file, "rb") as bam_in, \
     pysam.AlignmentFile(output_bam_file, "wb", header=bam_in.header) as bam_out:
    
    for read in bam_in:
        if read.is_unmapped or not read.query_sequence:
            continue
        if not contains_snp(read, snp_positions):
            bam_out.write(read)

print("SNP filtering complete.")
EOF

samtools index ${PREFIX}_SNP_filtered.bam

# ---------------------------------------------------------
# 6. Isolate Foldback Inversion Candidates
# ---------------------------------------------------------
echo "Detecting foldback inversions..."
python <<EOF
import pysam
from collections import defaultdict

# FIX: Moved the $ to properly expand the PREFIX variable
bam_file = "chr4_reads_${PREFIX}.bam"
output_bam = "${PREFIX}_foldbacks_unsorted.bam" # Writing as unsorted first
MIN_LENGTH_RATIO = 1.5

def is_real_foldback(alignments):
    if len(alignments) != 2:
        return False
    
    aln1, aln2 = alignments
    if aln1.is_reverse == aln2.is_reverse:
        return False
    if aln1.reference_name != aln2.reference_name:
        return False

    len1 = aln1.reference_end - aln1.reference_start
    len2 = aln2.reference_end - aln2.reference_start
    
    ratio = max(len1, len2) / max(min(len1, len2), 1)
    
    return ratio > MIN_LENGTH_RATIO

reads_by_name = defaultdict(list)
with pysam.AlignmentFile(bam_file, "rb") as bam:
    for read in bam:
        if read.is_supplementary or read.has_tag('SA'):
            reads_by_name[read.query_name].append(read)

with pysam.AlignmentFile(bam_file, "rb") as bam_in, \
     pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:
    
    for qname, alignments in reads_by_name.items():
        if is_real_foldback(alignments):
            for aln in alignments:
                bam_out.write(aln)

print("Foldback inversion isolation complete.")
EOF

# Coordinate-sort the foldbacks so IGV can calculate coverage
echo "Sorting and indexing foldback BAM..."
samtools sort -@ 20 -o ${PREFIX}_foldbacks.bam ${PREFIX}_foldbacks_unsorted.bam
samtools index ${PREFIX}_foldbacks.bam
rm ${PREFIX}_foldbacks_unsorted.bam

# ---------------------------------------------------------
# 7. Structural Variant Calling
# ---------------------------------------------------------
echo "Running SV callers on pre-filtered and post-filtered BAMs..."

ml Sniffles/2.4-GCC-13.3.0
# Sniffles on pre-SNP-filtered (MapQ >= 1 only)
sniffles --input chr4_reads_${PREFIX}.bam --vcf ${PREFIX}_chr4_sniffles.vcf --threads 20
# Sniffles on post-SNP-filtered
sniffles --input ${PREFIX}_SNP_filtered.bam --vcf ${PREFIX}_SNP_filtered_sniffles.vcf --threads 20

ml Miniforge3/24.11.3-0
source activate severus_env
# Severus on pre-SNP-filtered (MapQ >= 1 only)
severus --target-bam chr4_reads_${PREFIX}.bam --out-dir sev_${PREFIX}_chr4_out --threads 20
# Severus on post-SNP-filtered
severus --target-bam ${PREFIX}_SNP_filtered.bam --out-dir sev_${PREFIX}_SNP_filtered_out --threads 20

echo "Pipeline finished successfully!"
