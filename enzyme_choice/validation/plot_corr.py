import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

# Load data
file_path = sys.argv[1]
df = pd.read_csv(file_path)

# Keep relevant columns and ensure numeric
cols = [
    'Centromere Mean Spacing (bp)', 'mean_monomer_length_CEN',
    'Arm Mean Spacing (bp)', 'mean_monomer_length_ARMS'
]
df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')
df = df.dropna(subset=cols + ['Enzyme'])

# Rename for plotting
df_cen = df.rename(columns={
    'Centromere Mean Spacing (bp)': 'In silico mean monomer length (bp)',
    'mean_monomer_length_CEN': 'Experimental Mean Monomer Length'
})

df_arm = df.rename(columns={
    'Arm Mean Spacing (bp)': 'In silico mean monomer length (bp)',
    'mean_monomer_length_ARMS': 'Experimental Mean Monomer Length'
})

# Plot settings
sns.set(style="whitegrid")
point_size = 100

# Function to compute correlation labels
def get_corr_labels(df, region_label):
    label_map = {}
    for enzyme in sorted(df['Enzyme'].unique()):
        sub = df[df['Enzyme'] == enzyme]
        if len(sub) >= 2:
            x = sub['In silico mean monomer length (bp)']
            y = sub['Experimental Mean Monomer Length']
            spearman_corr, _ = spearmanr(x, y)
            pearson_corr, _ = pearsonr(x, y)
            label = f"{enzyme}\nPearson: {pearson_corr:.2f}\nSpearman: {spearman_corr:.2f}"
        else:
            label = f"{enzyme}\n(n=1)"
        label_map[enzyme] = label
    df['Enzyme_Label'] = df['Enzyme'].map(label_map)
    return df

# Apply correlation labels
df_cen = get_corr_labels(df_cen, "Centromere")
df_arm = get_corr_labels(df_arm, "Arm")

# Plot for CEN
plt.figure(figsize=(9, 7))
sns.scatterplot(
    data=df_cen,
    x='In silico mean monomer length (bp)',
    y='Experimental Mean Monomer Length',
    hue='Enzyme_Label',
    s=point_size
)
lims = [
    min(df_cen['In silico mean monomer length (bp)'].min(), df_cen['Experimental Mean Monomer Length'].min()),
    max(df_cen['In silico mean monomer length (bp)'].max(), df_cen['Experimental Mean Monomer Length'].max())
]
plt.plot(lims, lims, 'k--', alpha=0.7)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("In silico mean monomer length (bp)")
plt.ylabel("Experimental Mean Monomer Length (CEN)")
plt.title("Centromere: In Silico vs Experimental")
plt.legend(title="Enzyme (with correlations)", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

# Plot for ARMS
plt.figure(figsize=(9, 7))
sns.scatterplot(
    data=df_arm,
    x='In silico mean monomer length (bp)',
    y='Experimental Mean Monomer Length',
    hue='Enzyme_Label',
    s=point_size
)
lims = [
    min(df_arm['In silico mean monomer length (bp)'].min(), df_arm['Experimental Mean Monomer Length'].min()),
    max(df_arm['In silico mean monomer length (bp)'].max(), df_arm['Experimental Mean Monomer Length'].max())
]
plt.plot(lims, lims, 'k--', alpha=0.7)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("In silico mean monomer length (bp)")
plt.ylabel("Experimental Mean Monomer Length (ARMS)")
plt.title("Arms: In Silico vs Experimental")
plt.legend(title="Enzyme (with correlations)", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

