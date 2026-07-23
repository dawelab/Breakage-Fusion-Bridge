import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import chisquare

# Force publication-standard sans-serif font
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'

# ==========================================
# 1. USER CONFIGURATION
# ==========================================
CSV_FILE_1 = 'annotated_all_breaks.csv' 
CSV_FILE_2 = 'w22_feature_annotations.csv'

# Background region files
GENES_GFF = '/Users/kempton/Desktop/IGV_Resources/Genomes/ABSxc2/AbsGenomePBHIFI_version_1_liftoffA188.sorted.gff3'
TE_GFF = '/Users/kempton/Desktop/IGV_Resources/Genomes/ABSxc2/AbsGenomePBHIFI_version_1.fa.mod.EDTA.intact.gff3'

CHROM = 'chr4_RagTag'
REGION_START = 110000000
REGION_END = 190000000
REGION_LENGTH = REGION_END - REGION_START + 1

# --- ABS REPEAT REGION ---
ABS_REPEAT_START = 189006860
ABS_REPEAT_END = 189588809

HARDCODED_BREAKS = []

# Category Order and Color Mapping (4 Categories)
STANDARD_ORDER = [
    'Gene',
    'TE Only',
    'ABS Repeat',
    'Intergenic / Neither'
]

COLOR_MAP = {
    'Gene': '#DD8452',                 # Orange
    'TE Only': '#55A868',              # Green
    'ABS Repeat': '#8172B3',           # Purple
    'Intergenic / Neither': '#C44E52'  # Red
}

# ==========================================
# 2. DATA PROCESSING FUNCTIONS
# ==========================================
def simplify_category(cat_str, pos_bp=None):
    if pd.isna(cat_str) or str(cat_str).upper() == "UNMAPPED":
        return None
    
    cat_str = str(cat_str)
    
    if "Gene" in cat_str or "Both" in cat_str:
        return "Gene"
    if "TE Only" in cat_str or "TE Intact" in cat_str:
        return "TE Only"
    
    if pos_bp is not None and pd.notna(pos_bp):
        if ABS_REPEAT_START <= pos_bp <= ABS_REPEAT_END:
            return "ABS Repeat"
            
    return "Intergenic / Neither"

def get_combined_sample_counts():
    print(f"Loading sample data from {CSV_FILE_1} and {CSV_FILE_2}...")
    
    # --- Load CSV 1 ---
    df1 = pd.read_csv(CSV_FILE_1)
    pos_col1 = 'pos' if 'pos' in df1.columns else ('w22_location' if 'w22_location' in df1.columns else 'position')
    df1['pos_bp'] = pd.to_numeric(df1[pos_col1], errors='coerce')
    df1['simplified_cat'] = df1.apply(lambda r: simplify_category(r['category'], r['pos_bp']), axis=1)
    counts1 = df1['simplified_cat'].value_counts()
    
    # --- Load CSV 2 (W22 Format) ---
    df2 = pd.read_csv(CSV_FILE_2)
    pos_col2 = 'w22_location' if 'w22_location' in df2.columns else ('pos' if 'pos' in df2.columns else 'position')
    df2['pos_bp'] = pd.to_numeric(df2[pos_col2], errors='coerce')
    df2['simplified_cat'] = df2.apply(lambda r: simplify_category(r['genomic_context'], r['pos_bp']), axis=1)
    counts2 = df2['simplified_cat'].value_counts()
    
    # --- Hardcoded Breaks ---
    hardcoded_mapped = [simplify_category(cat, pos) for pos, cat in HARDCODED_BREAKS]
    counts_hardcoded = pd.Series(hardcoded_mapped).value_counts()
    
    combined_series = counts1.add(counts2, fill_value=0).add(counts_hardcoded, fill_value=0)
    
    return combined_series.reindex(STANDARD_ORDER, fill_value=0)

def paint_features_on_mask(gff_file, mask, chrom, reg_start, reg_end):
    gff_cols = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    df = pd.read_csv(gff_file, sep='\t', comment='#', names=gff_cols)
    
    df = df[df['seqid'] == chrom]
    
    for _, row in df.iterrows():
        if row['end'] >= reg_start and row['start'] <= reg_end:
            eff_start = max(reg_start, row['start']) - reg_start
            eff_end = min(reg_end, row['end']) - reg_start + 1
            mask[eff_start:eff_end] = True
            
    return mask

def get_control_counts():
    print(f"Generating digital chromosome for background calculations ({CHROM}:{REGION_START}-{REGION_END})...")
    gene_mask = np.zeros(REGION_LENGTH, dtype=bool)
    te_mask = np.zeros(REGION_LENGTH, dtype=bool)
    abs_repeat_mask = np.zeros(REGION_LENGTH, dtype=bool)

    abs_start_idx = max(REGION_START, ABS_REPEAT_START) - REGION_START
    abs_end_idx = min(REGION_END, ABS_REPEAT_END) - REGION_START + 1
    abs_repeat_mask[abs_start_idx:abs_end_idx] = True

    print("Painting genes onto mask...")
    gene_mask = paint_features_on_mask(GENES_GFF, gene_mask, CHROM, REGION_START, REGION_END)

    print("Painting transposons onto mask...")
    te_mask = paint_features_on_mask(TE_GFF, te_mask, CHROM, REGION_START, REGION_END)

    bp_gene_total = np.sum(gene_mask)
    bp_te_only = np.sum(~gene_mask & te_mask)
    bp_abs_repeat = np.sum(~gene_mask & ~te_mask & abs_repeat_mask)
    bp_neither = np.sum(~gene_mask & ~te_mask & ~abs_repeat_mask)

    counts = pd.Series({
        'Gene': bp_gene_total,
        'TE Only': bp_te_only,
        'ABS Repeat': bp_abs_repeat,
        'Intergenic / Neither': bp_neither
    })
    
    return counts.reindex(STANDARD_ORDER, fill_value=0)

# ==========================================
# 3. VISUALIZATION FUNCTION
# ==========================================
def plot_pies(sample_counts, control_counts):
    print("Generating 4-category pie charts...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    
    # Chi-Square goodness-of-fit test
    total_samples = sample_counts.sum()
    control_proportions = control_counts / control_counts.sum()
    expected_frequencies = control_proportions * total_samples
    
    valid_idx = expected_frequencies > 0
    chi2_stat, p_val = chisquare(f_obs=sample_counts[valid_idx], f_exp=expected_frequencies[valid_idx])
    
    datasets = [
        (sample_counts, "Combined SV Breakpoints\n(Observed)", "A"),
        (control_counts, f"Expected Background Composition\n({CHROM}: {REGION_START:,} - {REGION_END:,})", "B")
    ]
    
    for ax, (data, title, label) in zip(axes, datasets):
        sizes = [val for val in data.values if val > 0]
        plot_colors = [COLOR_MAP[cat] for cat, val in zip(data.index, data.values) if val > 0]
        
        total = sum(sizes)
        if "Combined" in title:
            outside_labels = [f"{(val/total)*100:.1f}%\n({int(val)})" for val in sizes]
        else:
            outside_labels = [f"{(val/total)*100:.1f}%" for val in sizes]
        
        ax.pie(
            sizes, 
            labels=outside_labels, 
            colors=plot_colors, 
            startangle=140, 
            wedgeprops={'edgecolor': 'white', 'linewidth': 2},
            textprops={'fontsize': 12, 'fontweight': 'bold'}
        )
        ax.set_title(title, fontsize=16, pad=20)
        ax.text(-0.1, 1.05, label, transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

    # DIRECT UNROUNDED SCIENTIFIC NOTATION ON PLOT
    p_str = f"p = {p_val:.10e} (Chi2 = {chi2_stat:.2f})"
    fig.text(0.5, 0.90, f"Chi-Square Test vs. Expected: {p_str}", ha='center', fontsize=12, fontstyle='italic')

    # Legend
    legend_elements = [mpatches.Patch(facecolor=COLOR_MAP[cat], edgecolor='white', label=cat) for cat in STANDARD_ORDER]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=13, bbox_to_anchor=(0.5, 0.02))
        
    plt.tight_layout(rect=[0, 0.10, 1, 0.88])
    
    plt.savefig("all_samples_combined_abs_repeat_piecharts.png", dpi=300)
    plt.savefig("all_samples_combined_abs_repeat_piecharts.pdf", format='pdf')
    print("Grind complete! Saved as 'all_samples_combined_abs_repeat_piecharts.png' and '.pdf'.")

# ==========================================
# 4. EXECUTION
# ==========================================
if __name__ == "__main__":
    s_counts = get_combined_sample_counts()
    c_counts = get_control_counts()
    plot_pies(s_counts, c_counts)