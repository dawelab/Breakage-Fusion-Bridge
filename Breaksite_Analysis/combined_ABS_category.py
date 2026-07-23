import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Force publication-standard sans-serif font
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'

# ==========================================
# 1. USER CONFIGURATION & FILE PATHS
# ==========================================
BREAKPOINTS_CSV_1 = 'annotated_all_breaks.csv'
BREAKPOINTS_CSV_2 = 'w22_feature_annotations.csv'
TELOMERES_CSV = 'telomere_sites.csv'

CENTROMERES_MB = [110.0, 189.5]

# --- ABS REPEAT REGION COORDINATES (in bp) ---
ABS_REPEAT_START = 189006860
ABS_REPEAT_END = 189588809

# Set fixed, symmetrical coordinate boundaries (105 Mb to 194 Mb)
SHARED_X_MIN = 105.0 
SHARED_X_MAX = 194.0

# --- Color & Shape Definitions (Now 4 Categories) ---
BREAKPOINT_COLORS = {
    'Gene': '#DD8452',                 # Orange
    'TE Only': '#55A868',              # Green
    'ABS Repeat': '#8172B3',           # Purple (New Category!)
    'Intergenic / Neither': '#C44E52'  # Red
}

TELOMERE_COLORS = {
    'ABS': '#C44E52', # Red
    'c2': '#55A868'   # Green
}

TELOMERE_SHAPES = {
    'R': '>',
    'L': '<'
}

# ==========================================
# 2. LOAD & PREPROCESS DATA
# ==========================================
def simplify_category(cat_str, pos_bp=None):
    """Maps genomic contexts and reclassifies intergenic hits in ABS Repeat region."""
    if pd.isna(cat_str) or str(cat_str).upper() == "UNMAPPED":
        return None
    
    cat_str = str(cat_str)
    
    if "Gene" in cat_str or "Both" in cat_str:
        return "Gene"
    if "TE Only" in cat_str or "TE Intact" in cat_str:
        return "TE Only"
    
    # Check if Intergenic hit falls inside ABS Repeat region
    if pos_bp is not None and pd.notna(pos_bp):
        if ABS_REPEAT_START <= pos_bp <= ABS_REPEAT_END:
            return "ABS Repeat"
            
    return "Intergenic / Neither"

# Load Breakpoints File 1
df_b1 = pd.read_csv(BREAKPOINTS_CSV_1)
pos_col1 = 'pos' if 'pos' in df_b1.columns else ('w22_location' if 'w22_location' in df_b1.columns else 'position')
cat_col1 = 'category' if 'category' in df_b1.columns else 'genomic_context'

df_b1['pos_bp'] = pd.to_numeric(df_b1[pos_col1], errors='coerce')
df_b1['pos_Mb'] = df_b1['pos_bp'] / 1_000_000
df_b1['clean_cat'] = df_b1.apply(lambda r: simplify_category(r[cat_col1], r['pos_bp']), axis=1)
df_b1 = df_b1.dropna(subset=['pos_Mb', 'clean_cat'])[['pos_Mb', 'clean_cat']]

# Load Breakpoints File 2
df_b2 = pd.read_csv(BREAKPOINTS_CSV_2)
pos_col2 = 'w22_location' if 'w22_location' in df_b2.columns else ('pos' if 'pos' in df_b2.columns else 'position')
cat_col2 = 'genomic_context' if 'genomic_context' in df_b2.columns else 'category'

df_b2['pos_bp'] = pd.to_numeric(df_b2[pos_col2], errors='coerce')
df_b2['pos_Mb'] = df_b2['pos_bp'] / 1_000_000
df_b2['clean_cat'] = df_b2.apply(lambda r: simplify_category(r[cat_col2], r['pos_bp']), axis=1)
df_b2 = df_b2.dropna(subset=['pos_Mb', 'clean_cat'])[['pos_Mb', 'clean_cat']]

df_breaks = pd.concat([df_b1, df_b2], ignore_index=True)

# Load Telomeres
df_telo = pd.read_csv(TELOMERES_CSV)
pos_col_t = 'position' if 'position' in df_telo.columns else 'pos'
df_telo['pos_Mb'] = pd.to_numeric(df_telo[pos_col_t], errors='coerce') / 1_000_000

# Staggering logic for telomeres
df_telo = df_telo.sort_values('pos_Mb').reset_index(drop=True)
overlap_threshold = 2.5  
y_step_size = 0.08      

offsets = []
levels = {}
for pos in df_telo['pos_Mb']:
    lvl = 0
    while True:
        if lvl not in levels:
            break
        if abs(pos - levels[lvl]) > overlap_threshold:
            break
        lvl += 1
    levels[lvl] = pos
    if lvl == 0:
        y_off = 0
    elif lvl % 2 == 1:
        y_off = ((lvl + 1) // 2) * y_step_size
    else:
        y_off = -(lvl // 2) * y_step_size
    offsets.append(y_off)

df_telo['y_offset'] = offsets

# ==========================================
# 3. FIGURE 1: SV BREAKPOINT LOLLIPOPS (4 CATEGORIES)
# ==========================================
fig1, ax1 = plt.subplots(figsize=(13, 4))

# 1. Symmetric Backbone (zorder=1)
ax1.hlines(y=0, xmin=SHARED_X_MIN, xmax=SHARED_X_MAX, color='#E0E0E0', linewidth=8, zorder=1)

# 2. Lollipops (zorder=3 and 4)
lollipop_top = 0.2
stem_width = 2.0
dot_size = 55

for _, row in df_breaks.iterrows():
    cat = row['clean_cat']
    x_val = row['pos_Mb']
    color = BREAKPOINT_COLORS.get(cat, '#8172B3')
    ax1.vlines(x=x_val, ymin=0, ymax=lollipop_top, color=color, linewidth=stem_width, alpha=0.8, zorder=3)
    ax1.scatter(x_val, lollipop_top, color=color, s=dot_size, zorder=4, edgecolors='white', linewidth=0.8)

# 3. Centromeres (zorder=5)
ax1.scatter(CENTROMERES_MB, [0, 0], color='white', edgecolors='black', marker='o', s=180, zorder=5, linewidth=1.5)

# Formatting
ax1.set_yticks([])  
ax1.set_xlim(SHARED_X_MIN, SHARED_X_MAX)
ax1.set_ylim(-0.1, lollipop_top + 0.2)

ax1.set_xlabel('Genomic Position (Mb)', fontsize=16, labelpad=10)
ax1.tick_params(axis='x', labelsize=14)

ax1.set_title('Structural Variant Breakpoint Contexts', fontsize=15, pad=15)

for spine in ['top', 'right', 'left']:
    ax1.spines[spine].set_visible(False)
ax1.tick_params(axis='y', length=0)

# 4-Category Legend
legend_breaks = [Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10, label=cat) for cat, c in BREAKPOINT_COLORS.items()]
ax1.legend(handles=legend_breaks, title="Breakpoint Context", bbox_to_anchor=(1.02, 1.0), loc='upper left')

fig1.subplots_adjust(left=0.05, right=0.78, top=0.85, bottom=0.20)

fig1.savefig("figure1_sv_lollipops.png", dpi=300, bbox_inches='tight')
fig1.savefig("figure1_sv_lollipops.pdf", format='pdf', bbox_inches='tight')

# ==========================================
# 4. FIGURE 2: TELOMERE SITES
# ==========================================
fig2, ax2 = plt.subplots(figsize=(13, 3.5))

# 1. Background stems (zorder=0)
for _, row in df_telo.iterrows():
    x_val = row['pos_Mb']
    y_off = row['y_offset']
    if y_off != 0:
        ax2.vlines(x=x_val, ymin=0, ymax=y_off, color='#A0A0A0', linewidth=1.2, zorder=0)

# 2. Symmetric Thick Backbone (zorder=1)
ax2.hlines(y=0, xmin=SHARED_X_MIN, xmax=SHARED_X_MAX, color='#E0E0E0', linewidth=12, zorder=1)

# 3. Telomere Triangles (zorder=4)
for _, row in df_telo.iterrows():
    genome = row['genome']
    direction = row['direction']
    x_val = row['pos_Mb']
    y_off = row['y_offset']
    
    if pd.isna(direction) or direction not in TELOMERE_SHAPES:
        continue
    color = TELOMERE_COLORS.get(genome, 'black')
    marker = TELOMERE_SHAPES[direction]
    ax2.scatter(x_val, y_off, color=color, marker=marker, s=140, zorder=4, edgecolors='white', linewidth=0.8)

# 4. Centromeres (zorder=5)
ax2.scatter(CENTROMERES_MB, [0, 0], color='white', edgecolors='black', marker='o', s=180, zorder=5, linewidth=1.5)

# Formatting
ax2.set_yticks([])  
ax2.set_xlim(SHARED_X_MIN, SHARED_X_MAX)
ax2.set_ylim(-0.35, 0.35)

ax2.set_xlabel('Genomic Position (Mb)', fontsize=16, labelpad=10)
ax2.tick_params(axis='x', labelsize=14)

ax2.set_title('Telomere Formation Sites', fontsize=15, pad=15)

for spine in ['top', 'right', 'left']:
    ax2.spines[spine].set_visible(False)
ax2.tick_params(axis='y', length=0)

legend_telo = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='white', markeredgecolor='black', markersize=10, label='Centromere'),
    Line2D([0], [0], marker='>', color='w', markerfacecolor=TELOMERE_COLORS['ABS'], markersize=10, label='ABS Genome'),
    Line2D([0], [0], marker='>', color='w', markerfacecolor=TELOMERE_COLORS['c2'], markersize=10, label='c2 Genome'),
    Line2D([0], [0], marker='>', color='w', markerfacecolor='grey', markersize=10, label='Direction (R)'),
    Line2D([0], [0], marker='<', color='w', markerfacecolor='grey', markersize=10, label='Direction (L)')
]
ax2.legend(handles=legend_telo, title="Telomere & Feature Details", bbox_to_anchor=(1.02, 1.0), loc='upper left')

fig2.subplots_adjust(left=0.05, right=0.78, top=0.85, bottom=0.20)

fig2.savefig("figure2_telomeres.png", dpi=300, bbox_inches='tight')
fig2.savefig("figure2_telomeres.pdf", format='pdf', bbox_inches='tight')

print("Execution complete. Both plots rendered with 4 categories for breakpoints.")
