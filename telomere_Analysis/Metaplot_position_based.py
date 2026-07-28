import os
import pandas as pd
from pathlib import Path

base_dir = Path(".")           # where the script is run
dirs_to_search = [base_dir / "Telo_Left", base_dir / "Telo_Right"]

dfs = {}   # dictionary to store all dataframes

for parent in dirs_to_search:
    for subdir in parent.iterdir():
        if not subdir.is_dir():
            continue

        dirname = subdir.name  # this is the unique identifier like "4a3_345"

        # Find the files inside this subdirectory
        for file in subdir.iterdir():
            if not file.name.endswith("_modkit.bed"):
                continue

            fname = file.name

            # telo file
            if "_telo_" in fname:
                key = f"{dirname}_telo"
                dfs[key] = pd.read_csv(file, sep="\t", header=None)

            # cont file (ABS or c2)
            elif "_ABS_" in fname or "_c2_" in fname:
                key = f"{dirname}_cont"
                dfs[key] = pd.read_csv(file, sep="\t", header=None)

cleaned_dfs = {}

for name, df in dfs.items():
    pos = df[2]
    coverage = df[4]
    frac_mod = df[10]

    cleaned = (
        pd.DataFrame({
            "pos": pos,
            "coverage": coverage,
            "frac_mod": frac_mod
        })
        .sort_values("pos")
        .reset_index(drop=True)
    )

    cleaned_dfs[name] = cleaned

merged_dfs = {}


# group dataframe names by their directory prefix
from collections import defaultdict

groups = defaultdict(dict)

for name, df in cleaned_dfs.items():
    dirname, kind = name.rsplit("_", 1)  # splits into e.g. "4a3_345" and "telo" or "cont"
    groups[dirname][kind] = df

# now merge telo + cont for each directory
for dirname, sub in groups.items():
    telo = sub.get("telo")
    cont = sub.get("cont")

    if telo is None:
        print(f"Warning: no telo file for {dirname}")
        continue

    if cont is None:
        print(f"Warning: no control file for {dirname}")
        continue

    merged = pd.merge(
        telo, cont,
        on="pos",
        how="left",
        suffixes=("_telo", "_cont")
    )

    merged_dfs[dirname] = merged

print(merged_dfs.keys())

#I need to fix areas where some rows get duplicated.Here I am dropping whichever row has lower coverage for the telomere and if telomere is equal then I drop whichever has lower coverage for control. This seems to be an error when a strand there is a misscalled based on the opposite strand. For example if I have 8 methylated cytosines at position 5 on strand + but a miscalled cytosine at position 5 on strand - it will end up making two rows here. The lower covergae rule should be almost if not perfect at removing all the samples with this type of problem

deduped_merged_dfs = {}

for dirname, df in merged_dfs.items():
    # Sort so the "best" row per collapsed site is FIRST:
    df_sorted = df.sort_values(
        by=["coverage_telo", "coverage_cont"],
        ascending=[False, False]
    )

    # Drop duplicates by position (keep the best row)
    df_dedup = df_sorted.drop_duplicates(subset=["pos"], keep="first")

    # Optional: sort output back by position
    df_dedup = df_dedup.sort_values("pos").reset_index(drop=True)

    deduped_merged_dfs[dirname] = df_dedup

merged_dfs = deduped_merged_dfs

from pathlib import Path

#This is where I will fill in all the empty rows. So currently I only have rows for genomic positions where there is a cytosine. This part will add the position and fill in coverage and fraction modified with NaN

for key, df in merged_dfs.items():

    # Ensure sorted by genomic position
    df = df.sort_values("pos")

    # Build full continuous position range
    full_range = range(df["pos"].min(), df["pos"].max() + 1)

    # Set pos as index → reindex → fill missing with NaN
    df_filled = (
        df.set_index("pos")
          .reindex(full_range)
          .reset_index()
          .rename(columns={"index": "pos"})
    )

    # Save back into the dictionary
    merged_dfs[key] = df_filled

# Get the list of directory names under A/
R_dirs = {d.name for d in Path("Telo_Right").iterdir() if d.is_dir()}

# Reverse merged dataframes that belong to directories in R_dir/
#This makes the telomere end go to the bottom of the plot. 
# I had to switch this to R_dir. When it was left dir the shorter reads would end in the middle of the plot which messed it up. Now telomeres begin on the left side so the top of my dataframe is directly at the addition site. 
for dirname, df in merged_dfs.items():
    if dirname in R_dirs:
        merged_dfs[dirname] = df.iloc[::-1].reset_index(drop=True)
        
first100_dfs = {}

###Modify the number in df.iloc to change length for analysis

for name, df in merged_dfs.items():
    first100_dfs[name] = df.iloc[:10000].reset_index(drop=True)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", context="talk")

plt.figure(figsize=(14, 8))

window = 100
smoothed_arrays = {}
max_len = 0

# ---- Smooth curves and track max length ----
for name, df in first100_dfs.items():
    delta = df["frac_mod_telo"] -  df["frac_mod_cont"]
    delta_smooth = pd.Series(delta).rolling(window=window, center=True, min_periods=1).mean()

    smoothed_arrays[name] = delta_smooth.values
    max_len = max(max_len, len(delta_smooth))

# ---- Pad everything with NaN so all arrays have equal length ----
padded = []
for name, arr in smoothed_arrays.items():
    pad_len = max_len - len(arr)
    arr_padded = np.pad(arr, (0, pad_len), constant_values=np.nan)
    padded.append(arr_padded)

    # plot individual curves (truncate to actual length)
    plt.plot(range(len(arr)), arr, linewidth=1, alpha=1, label=name)

# ---- Compute mean ignoring NaN ----
padded = np.vstack(padded)
mean_curve = np.nanmean(padded, axis=0)

plt.plot(
    range(max_len),
    mean_curve,
    linewidth=4,
    color="black",
    label="Mean Δ (cont − telo, smoothed)"
)

plt.xlabel("Nucleotides after Telomere Site")
plt.xlim(0,5000)
plt.ylabel("Δ methylation (telo - cont)")
plt.title(f"Smoothed Methylation Difference in {window} bins")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
#plt.show()
plt.savefig("All_Nucleotides.png", dpi=1200)