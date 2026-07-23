import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import chisquare

# Force publication-standard sans-serif font
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'

# ==========================================
# 1. USER CONFIGURATION
# ==========================================
ANNOTATED_CSV = "w22_feature_annotations.csv"
GENE_GFF3 = "/Users/kempton/Desktop/Telomere_Paper/breaksite_enrichment/w22/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3"
TE_GFF3 = "/Users/kempton/Desktop/Telomere_Paper/breaksite_enrichment/w22/W22.structuralTEv2.disjoined.2018-09-22.gff3"

CHROM = "chr4" 
CONTROL_START = 110000000
CONTROL_END = 190000000

# Strict Category Order and Color Mapping (Identical to previous plot)
STANDARD_ORDER = [
    'Gene Only',
    'TE Only',
    'Both (Gene & TE)',
    'Intergenic / Neither'
]

COLOR_MAP = {
    'Gene Only': '#DD8452',            # Orange
    'TE Only': '#55A868',              # Green
    'Both (Gene & TE)': '#4C72B0',     # Blue
    'Intergenic / Neither': '#C44E52'  # Red
}

# ==========================================
# 2. DATA PROCESSING FUNCTIONS
# ==========================================
def load_annotations():
    print("Loading gene and TE annotations for W22...")
    genes_pr = pr.read_gff3(GENE_GFF3)
    genes_pr = genes_pr[genes_pr.Feature == "gene"]
    
    tes_pr = pr.read_gff3(TE_GFF3)
    te_df = tes_pr.df
    
    def fix_chrom(c):
        c_str = str(c)
        return f"chr{c_str}" if c_str.isdigit() else c_str
        
    te_df["Chromosome"] = te_df["Chromosome"].apply(fix_chrom)
    fixed_tes_pr = pr.PyRanges(te_df)
    
    return genes_pr, fixed_tes_pr

def get_sample_counts():
    print(f"Loading sample data from {ANNOTATED_CSV}...")
    df = pd.read_csv(ANNOTATED_CSV)
    
    def simplify_context(context):
        if pd.isna(context) or context == "UNMAPPED":
            return None
        if "Gene & TE Intact" in context:
            return "Both (Gene & TE)"
        if "Gene" in context:
            return "Gene Only"
        if "TE Intact" in context:
            return "TE Only"
        return "Intergenic / Neither"

    df["broad_context"] = df["genomic_context"].apply(simplify_context)
    counts = df["broad_context"].value_counts()
    return counts.reindex(STANDARD_ORDER, fill_value=0)

def get_control_counts(genes_pr, tes_pr):
    print(f"Calculating background composition for W22 ({CHROM}:{CONTROL_START}-{CONTROL_END})...")
    window_size = CONTROL_END - CONTROL_START
    window_pr = pr.from_dict({"Chromosome": [CHROM], "Start": [CONTROL_START], "End": [CONTROL_END]})
    
    te_df = tes_pr.df
    if "intact" in te_df.columns:
        intact_df = te_df[te_df["intact"].astype(str).str.upper() == "TRUE"]
        tes_intact_pr = pr.PyRanges(intact_df) if not intact_df.empty else pr.PyRanges()
    else:
        tes_intact_pr = tes_pr

    genes_in_window = genes_pr.intersect(window_pr).merge()
    tes_intact = tes_intact_pr.intersect(window_pr).merge() if not tes_intact_pr.empty else pr.PyRanges()
    
    # Calculate overlaps
    overlap_pr = genes_in_window.intersect(tes_intact)
    overlap_bp = sum(overlap_pr.lengths()) if not overlap_pr.empty else 0
    
    genes_only = genes_in_window.subtract(tes_intact) if not overlap_pr.empty else genes_in_window
    gene_only_bp = sum(genes_only.lengths()) if not genes_only.empty else 0
    
    tes_only = tes_intact.subtract(genes_in_window) if not overlap_pr.empty else tes_intact
    te_only_bp = sum(tes_only.lengths()) if not tes_only.empty else 0
            
    intergenic_bp = window_size - overlap_bp - gene_only_bp - te_only_bp
    
    counts = pd.Series({
        "Gene Only": gene_only_bp, 
        "TE Only": te_only_bp, 
        "Both (Gene & TE)": overlap_bp,
        "Intergenic / Neither": intergenic_bp
    })
    return counts.reindex(STANDARD_ORDER, fill_value=0)

# ==========================================
# 3. VISUALIZATION FUNCTION
# ==========================================
def plot_pies(sample_counts, control_counts):
    print("Generating side-by-side pie charts...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    
    # Run Chi-Square goodness-of-fit test
    total_samples = sample_counts.sum()
    control_proportions = control_counts / control_counts.sum()
    expected_frequencies = control_proportions * total_samples
    
    # Calculate p-value (avoid zero division)
    valid_idx = expected_frequencies > 0
    chi2_stat, p_val = chisquare(f_obs=sample_counts[valid_idx], f_exp=expected_frequencies[valid_idx])
    
    datasets = [
        (sample_counts, "Sample SV Breakpoints\n(Observed)", "C"),
        (control_counts, f"Expected Background Composition\n({CHROM}: {CONTROL_START:,} - {CONTROL_END:,})", "D")
    ]
    
    for ax, (data, title, label) in zip(axes, datasets):
        sizes = [val for val in data.values if val > 0]
        plot_colors = [COLOR_MAP[cat] for cat, val in zip(data.index, data.values) if val > 0]
        
        # Calculate outside text labels
        total = sum(sizes)
        if "Sample" in title:
            outside_labels = [f"{(val/total)*100:.1f}%\n({val})" for val in sizes]
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
        
        # Add Panel Letter (A / B) for publication layout
        ax.text(-0.1, 1.05, label, transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

    # Add Chi-square p-value text onto figure
    p_str = "p < 0.001" if p_val < 0.001 else f"p = {p_val:.4f}"
    fig.text(0.5, 0.90, f"Chi-Square Test vs. Expected: {p_str}", ha='center', fontsize=12, fontstyle='italic')

    # Create the custom universal legend
    legend_elements = [mpatches.Patch(facecolor=COLOR_MAP[cat], edgecolor='white', label=cat) for cat in STANDARD_ORDER]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=14, bbox_to_anchor=(0.5, 0.02))
        
    plt.tight_layout(rect=[0, 0.10, 1, 0.88])
    
    # Save both PNG (for presentations) and PDF (for journal submission)
    plt.savefig("w22_combined_breaksite_piecharts.png", dpi=300)
    plt.savefig("w22_combined_breaksite_piecharts.pdf", format='pdf')
    print("Complete! Saved as 'w22_combined_breaksite_piecharts.png' and '.pdf'.")

# ==========================================
# 4. EXECUTION
# ==========================================
if __name__ == "__main__":
    s_counts = get_sample_counts()
    genes, tes = load_annotations()
    c_counts = get_control_counts(genes, tes)
    
    plot_pies(s_counts, c_counts)