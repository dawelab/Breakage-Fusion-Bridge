import pandas as pd
import matplotlib.pyplot as plt
import glob

# ==========================================
# 1. SETUP: Provide paths and hardcoded sites
# ==========================================
# 1a. BEDPE file(s) - use glob pattern if you have multiple, or direct path for one
BEDPE_FILE = 'cat.bedpe' 

# 1b. Hardcoded breaks
my_hardcoded_breaks = [
    {'chrom': 'chr4_RagTag', 'pos': 178101389},
    {'chrom': 'chr4_RagTag', 'pos': 161404502},
    {'chrom': 'chr4_RagTag', 'pos': 151760101},
    {'chrom': 'chr4_RagTag', 'pos': 186529890},
    {'chrom': 'chr4_RagTag', 'pos': 116023409},
    {'chrom': 'chr4_RagTag', 'pos': 130346197},
    {'chrom': 'chr4_RagTag', 'pos': 172716549},
    {'chrom': 'chr4_RagTag', 'pos': 127343146},
    {'chrom': 'chr4_RagTag', 'pos': 131357204},
    {'chrom': 'chr4_RagTag', 'pos': 163371467},
    {'chrom': 'chr4_RagTag', 'pos': 174663704},
    {'chrom': 'chr4_RagTag', 'pos': 174542727}, #4b9.d1 foldback
    {'chrom': 'chr4_RagTag', 'pos': 178122315}, #4a1 foldback
    {'chrom': 'chr4_RagTag', 'pos': 125933390},
    {'chrom': 'chr4_RagTag', 'pos': 136438699},
    {'chrom': 'chr4_RagTag', 'pos': 186426084}
]

# 1c. Your feature files
GENES_GFF = '/Users/kempton/Desktop/IGV_Resources/Genomes/ABSxc2/AbsGenomePBHIFI_version_1_liftoffA188.sorted.gff3'
TE_GFF = '/Users/kempton/Desktop/IGV_Resources/Genomes/ABSxc2/AbsGenomePBHIFI_version_1.fa.mod.EDTA.intact.gff3'

# ==========================================
# 2. PARSE AND COMBINE ALL BREAKPOINTS
# ==========================================
print("Loading BEDPE files and hardcoded sites...")

# --- Process BEDPE ---
# If using multiple bedpe files with glob, you can loop this, but for one file:
bedpe = pd.read_csv(BEDPE_FILE, sep='\t', header=None, comment='#')

# Extract breakpoints (adding +1 to convert 0-based to 1-based coordinates)
breaks1 = pd.DataFrame({'chrom': bedpe[0], 'pos': bedpe[1] + 1})
breaks2 = pd.DataFrame({'chrom': bedpe[3], 'pos': bedpe[4] + 1})

# Combine into a single dataframe for BEDPE breaks
bedpe_breaks = pd.concat([breaks1, breaks2])

# --- Process Hardcoded ---
hardcoded_breaks = pd.DataFrame(my_hardcoded_breaks)

# --- Merge Everything ---
# Combine BEDPE and hardcoded breaks into one master list
all_breaks = pd.concat([bedpe_breaks, hardcoded_breaks]).drop_duplicates().reset_index(drop=True)

print(f"Total unique breakpoints to analyze: {len(all_breaks)}")

# ==========================================
# 3. PARSE THE GFF3 FILES
# ==========================================
gff_cols = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

print("Loading GFF3 files...")
genes_df = pd.read_csv(GENES_GFF, sep='\t', comment='#', names=gff_cols)
tes_df = pd.read_csv(TE_GFF, sep='\t', comment='#', names=gff_cols)

# ==========================================
# 4. INTERSECT BREAKS WITH FEATURES
# ==========================================
def check_overlap(breaks, features):
    """Checks if each breakpoint falls within any feature interval."""
    overlaps = []
    for _, row in breaks.iterrows():
        # Find any feature on the same chrom where start <= break_pos <= end
        match = features[
            (features['seqid'] == row['chrom']) & 
            (features['start'] <= row['pos']) & 
            (features['end'] >= row['pos'])
        ]
        overlaps.append(not match.empty)
    return overlaps

print("Intersecting with genes...")
all_breaks['hit_gene'] = check_overlap(all_breaks, genes_df)

print("Intersecting with transposons...")
all_breaks['hit_te'] = check_overlap(all_breaks, tes_df)

# Categorize each breakpoint
def categorize(row):
    if row['hit_gene'] and row['hit_te']:
        return 'Both (Gene & TE)'
    elif row['hit_gene']:
        return 'Gene Only'
    elif row['hit_te']:
        return 'TE Only'
    else:
        return 'Intergenic / Neither'

all_breaks['category'] = all_breaks.apply(categorize, axis=1)

# ==========================================
# 5. VISUALIZE THE RESULTS
# ==========================================
# Count how many breaks fall into each category
counts = all_breaks['category'].value_counts()

# Create a bar chart
plt.figure(figsize=(8, 6))
bars = counts.plot(kind='bar', color=['#4C72B0', '#DD8452', '#55A868', '#C44E52'][:len(counts)])

plt.title('All BFB Breakpoint Overlaps with Genomic Features', fontsize=14)
plt.ylabel('Number of Breakpoints', fontsize=12)
plt.xticks(rotation=45, ha='right')

# Add the exact numbers on top of the bars
for p in bars.patches:
    bars.annotate(str(int(p.get_height())), 
                  (p.get_x() + p.get_width() / 2., p.get_height()), 
                  ha='center', va='bottom', xytext=(0, 5), 
                  textcoords='offset points')

plt.tight_layout()
plt.show()

# Save the final annotated table
all_breaks.to_csv('annotated_all_breaks.csv', index=False)
print("Finished! Saved raw overlap data to 'annotated_all_breaks.csv'")