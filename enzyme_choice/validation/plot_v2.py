import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

# Read input file from sys.argv
file_path = sys.argv[1]
df = pd.read_csv(file_path)

# Ensure numeric
cols = [
    'Centromere Mean Spacing (bp)', 'mean_monomer_length_CEN',
    'Arm Mean Spacing (bp)', 'mean_monomer_length_ARMS'
]
df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')
df = df.dropna(subset=cols + ['Enzyme'])

# Rename columns for plotting
df_cen = df.rename(columns={
    'Centromere Mean Spacing (bp)': 'In silico mean monomer length (bp)',
    'mean_monomer_length_CEN': 'Experimental Mean Monomer Length (CEN)'
})

df_arm = df.rename(columns={
    'Arm Mean Spacing (bp)': 'In silico mean monomer length (bp)',
    'mean_monomer_length_ARMS': 'Experimental Mean Monomer Length (ARMS)'
})

# Plot style
sns.set(style="whitegrid")

# Correlation text function
def get_corr_text(x, y, label):
    spearman_corr, _ = spearmanr(x, y)
    pearson_corr, _ = pearsonr(x, y)
    return f"{label}\nSpearman: {spearman_corr:.2f}\nPearson: {pearson_corr:.2f}"

# Plot for Centromeres
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=df_cen,
    x='In silico mean monomer length (bp)',
    y='Experimental Mean Monomer Length (CEN)',
    hue='Enzyme',
    s=100
)
lims = [
    min(df_cen['In silico mean monomer length (bp)'].min(), df_cen['Experimental Mean Monomer Length (CEN)'].min()),
    max(df_cen['In silico mean monomer length (bp)'].max(), df_cen['Experimental Mean Monomer Length (CEN)'].max())
]
plt.plot(lims, lims, 'k--', alpha=0.7)
plt.xlim(0,10000)
plt.ylim(bottom=0)
plt.title("Centromere: In Silico vs. Experimental Monomer Length")
plt.xlabel("In silico mean monomer length (bp)")
plt.ylabel("Experimental Mean Monomer Length (CEN)")

# Add correlation text
for enzyme in df_cen['Enzyme'].unique():
    sub = df_cen[df_cen['Enzyme'] == enzyme]
    if len(sub) >= 2:
        x = sub['In silico mean monomer length (bp)']
        y = sub['Experimental Mean Monomer Length (CEN)']
        txt = get_corr_text(x, y, enzyme)
        plt.annotate(txt, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')

plt.legend(title='Enzyme')
plt.tight_layout()
plt.show()

# Plot for Arms
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=df_arm,
    x='In silico mean monomer length (bp)',
    y='Experimental Mean Monomer Length (ARMS)',
    hue='Enzyme',
    s=100
)
lims = [
    min(df_arm['In silico mean monomer length (bp)'].min(), df_arm['Experimental Mean Monomer Length (ARMS)'].min()),
    max(df_arm['In silico mean monomer length (bp)'].max(), df_arm['Experimental Mean Monomer Length (ARMS)'].max())
]
plt.plot(lims, lims, 'k--', alpha=0.7)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.title("Arms: In Silico vs. Experimental Monomer Length")
plt.xlabel("In silico mean monomer length (bp)")
plt.ylabel("Experimental Mean Monomer Length (ARMS)")

# Add correlation text
for enzyme in df_arm['Enzyme'].unique():
    sub = df_arm[df_arm['Enzyme'] == enzyme]
    if len(sub) >= 2:
        x = sub['In silico mean monomer length (bp)']
        y = sub['Experimental Mean Monomer Length (ARMS)']
        txt = get_corr_text(x, y, enzyme)
        plt.annotate(txt, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')

plt.legend(title='Enzyme')
plt.tight_layout()
plt.show()

