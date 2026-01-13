
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# ---------------------------
# Step 0: Setup
# ---------------------------
bench_dir = "bench_results"  # folder containing all benchmark files
num_map = {'G': 1e9, 'M': 1e6, 'K': 1e3}

def parse_number(s):
    """Parse numbers with G/M/K suffix or GB/s"""
    s = s.strip().replace(",", "").replace("GB/s", "")
    match = re.match(r"([\d\.]+)([GMK]?)", s)
    if match:
        val, suffix = match.groups()
        val = float(val)
        if suffix in num_map:
            val *= num_map[suffix]
        return val
    return float(s)

def extract_dna_len(dataset_name):
    """Extract numeric DNA length from synthetic datasets"""
    if dataset_name.startswith("synt") and dataset_name.endswith(".fa"):
        return int(dataset_name[4:-3])
    return np.nan

# ---------------------------
# Step 1: Parse all files
# ---------------------------
all_rows = []

for fname in os.listdir(bench_dir):
    
    dataset_name = fname.replace(".txt","")
    dna_length = extract_dna_len(dataset_name)
    
    with open(os.path.join(bench_dir,fname), "r") as f:
        content = f.read()
    
    # Split into experiment blocks
    blocks = re.split(r'\n(?=[^\s].*?:\n)', content)
    
    for block in blocks:
        lines = block.strip().splitlines()
        if not lines:
            continue
        impl_name = lines[0].strip().replace(":", "")
        metrics = {}
        for line in lines[1:]:
            if ':' in line:
                key, val = line.split(":",1)
                key = key.strip()
                val = val.strip().replace("%","")
                try:
                    metrics[key] = parse_number(val)
                except:
                    metrics[key] = np.nan
        row = {
            "dataset": dataset_name,
            "dna_length": dna_length,
            "implementation": impl_name,
            "throughput": metrics.get("throughput", np.nan),
            "cycles": metrics.get("cycles", np.nan),
            "cycles/byte": metrics.get("cycles/byte", np.nan),
            "instructions": metrics.get("instructions", np.nan),
            "instr/byte": metrics.get("instr/byte", np.nan),
            "branches": metrics.get("branches", np.nan),
            "branches/byte": metrics.get("branches/byte", np.nan),
            "branch_miss_percent": metrics.get("% branch miss", np.nan),
            "branch_misses": metrics.get("branch misses", np.nan)
        }
        all_rows.append(row)

# Convert to DataFrame
df = pd.DataFrame(all_rows)
df['dataset'] = df['dataset'].str.replace(r'\.bench$', '', regex=True)
# ---------------------------
# Step 2: Plotting
# ---------------------------
sns.set(style="whitegrid", font_scale=1.1)
palette = sns.color_palette("tab10")

# --- 1) Throughput by dataset & implementation ---
throughput_impls = [
    "Needletail (reader)",
    "Paraseq (reader)",
    "DNA string (slice)",
    "DNA packed (slice)",
    "DNA columnar (slice)"
]

df_throughput = df[df['implementation'].isin(throughput_impls)]
df_throughput = df_throughput[
    (~df_throughput['dataset'].str.contains(".gz")) & 
    (~df_throughput['dataset'].str.startswith("synt"))
]

# Optional: sort datasets alphabetically
df_throughput = df_throughput.sort_values(by='dataset')



# ---------------------------
# Compute speedup relative to Needletail (reader)
# ---------------------------

# First, filter as before (real datasets only, relevant implementations)
throughput_impls = [
    "Needletail (reader)",
    "Paraseq (reader)",
    "DNA string (slice)",
    "DNA packed (slice)",
    "DNA columnar (slice)"
]

df_throughput = df[df['implementation'].isin(throughput_impls)]
df_throughput = df_throughput[
    (~df_throughput['dataset'].str.endswith(".gz")) &
    (~df_throughput['dataset'].str.startswith("synt"))
]

# Pivot to compute speedup easily
df_pivot = df_throughput.pivot(index='dataset', columns='implementation', values='throughput')

# Compute speedup relative to Needletail (reader)
for col in df_pivot.columns:
    df_pivot[col + " speedup"] = df_pivot[col] / df_pivot["Needletail (reader)"]

# Melt back to long format for plotting
speedup_cols = [col for col in df_pivot.columns if "speedup" in col]
df_speedup = df_pivot[speedup_cols].reset_index()
df_speedup = df_speedup.melt(id_vars='dataset', var_name='implementation', value_name='speedup')
df_speedup['implementation'] = df_speedup['implementation'].str.replace(" speedup", "")

# ---------------------------
# Plot: Throughput + Speedup on same subplot
# ---------------------------

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# ---------------------------
# Filter for real datasets and selected implementations
# ---------------------------
throughput_impls = [
    "Needletail (reader)",
    "Paraseq (reader)",
    "DNA string (slice)",
    "DNA packed (slice)",
    "DNA columnar (slice)"
]

df_throughput = df[df['implementation'].isin(throughput_impls)]
df_throughput = df_throughput[
    (~df_throughput['dataset'].str.contains(".gz")) &
    (~df_throughput['dataset'].str.startswith("synt"))
]

# Sort datasets alphabetically for aligned plotting
df_throughput = df_throughput.sort_values(by='dataset')

# ---------------------------
# Compute speedup relative to Needletail
# ---------------------------
breakpoint()
df_pivot = df_throughput.pivot(index='dataset', columns='implementation', values='throughput')
df_speedup = df_pivot.divide(df_pivot["Needletail (reader)"], axis=0)
df_speedup = df_speedup.reset_index().melt(id_vars='dataset', var_name='implementation', value_name='speedup')

# ---------------------------
# Plotting two subplots
# ---------------------------

# Assume df_throughput and df_speedup are already prepared
sns.set(style="whitegrid", font_scale=1.1)
palette = sns.color_palette("tab10")
hue_order = ["Needletail (reader)", "Paraseq (reader)", "DNA string (slice)", "DNA packed (slice)", "DNA columnar (slice)"]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14,8), sharex=True)

# --- Top subplot: absolute throughput ---
sns.barplot(
    data=df_throughput,
    x='dataset',
    y='throughput',
    hue='implementation',
    palette=palette,
    hue_order=hue_order,
    ax=ax1,
    dodge=True
)
ax1.set_ylabel("Throughput (GB/s)")
ax1.set_title("Throughput by Dataset and Implementation")
ax1.legend_.remove()  # remove legend from top plot

# --- Bottom subplot: speedup relative to Needletail ---
sns.barplot(
    data=df_speedup,
    x='dataset',
    y='speedup',
    hue='implementation',
    palette=palette,
    hue_order=hue_order,
    ax=ax2,
    dodge=True
)
ax2.set_ylabel("Speedup vs Needletail")
ax2.axhline(1.0, color='gray', linestyle='--', linewidth=1)
ax2.legend_.remove()  # remove legend from bottom plot
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')

# --- Create a single shared legend above the figure ---
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=len(hue_order), frameon=False, title="Implementation")

plt.tight_layout(rect=[0,0,1,0.95])  # leave space on top for the legend
plt.savefig("plot_throughput_and_speedup_aligned_legend.png")
plt.show()

exit()



# --- 2) Throughput vs DNA length (syntX only) ---
df_synt = df.dropna(subset=['dna_length'])
plt.figure(figsize=(10,6))
sns.lineplot(data=df_synt, x='dna_length', y='throughput', hue='implementation', marker="o", palette=palette)
plt.xlabel("DNA length (syntX)")
plt.ylabel("Throughput (GB/s)")
plt.title("Throughput Scaling vs DNA Length (syntX)")
plt.xticks(sorted(df_synt['dna_length'].unique()))
plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
plt.tight_layout()
plt.savefig("bench_plots/plot2_throughput_vs_dna_length.png")
plt.close()

# --- 3) Cycles per byte ---
plt.figure(figsize=(14,6))
sns.barplot(data=df, x='dataset', y='cycles/byte', hue='implementation', palette=palette)
plt.xticks(rotation=45, ha='right')
plt.ylabel("Cycles per byte")
plt.title("Cycles per Byte by Dataset and Implementation")
plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
plt.tight_layout()
plt.savefig("bench_plots/plot3_cycles_per_byte.png")
plt.close()

# --- 4) Instructions per byte ---
plt.figure(figsize=(14,6))
sns.barplot(data=df, x='dataset', y='instr/byte', hue='implementation', palette=palette)
plt.xticks(rotation=45, ha='right')
plt.ylabel("Instructions per byte")
plt.title("Instructions per Byte by Dataset and Implementation")
plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
plt.tight_layout()
plt.savefig("bench_plots/plot4_instr_per_byte.png")
plt.close()

print("All 4 plots generated successfully.")
