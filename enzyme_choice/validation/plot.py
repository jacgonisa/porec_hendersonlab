import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read file from command-line argument
file_path = sys.argv[1]
df = pd.read_csv(file_path, sep='\t')

# Ensure relevant columns are numeric
cols = [
    'Centromere Mean Spacing (bp)', 'Centromere Spacing Std Dev',
    'mean_monomer_length_CEN', 'sd_monomer_length_CEN',
    'Arm Mean Spacing (bp)', 'Arm Spacing Std Dev',
    'mean_monomer_length_ARMS', 'sd_monomer_length_ARMS'
]
df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')

# Filter out rows with missing values in key columns
df = df.dropna(subset=cols + ['Enzyme'])

# Set plot style
sns.set(style="whitegrid")

# Centromere plot
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=df,
    x='Centromere Mean Spacing (bp)',
    y='mean_monomer_length_CEN',
    hue='Enzyme',
    s=100
)
# Error bars
#for _, row in df.iterrows():
#    plt.errorbar(
#        row['Centromere Mean Spacing (bp)'],
#        row['mean_monomer_length_CEN'],
#        xerr=row['Centromere Spacing Std Dev'],
#        yerr=row['sd_monomer_length_CEN'],
#        fmt='none', ecolor='gray', alpha=0.7, capsize=3
#    )

plt.title("Centromere: Predicted Spacing vs. Experimental Monomer Length")
plt.xlabel("Centromere Mean Spacing (bp)")
plt.ylabel("Mean Monomer Length (CEN)")
plt.xlim(0,10000)
plt.ylim(bottom=0)
plt.legend(title='Enzyme')
plt.tight_layout()
plt.show()

# Arm plot
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=df,
    x='Arm Mean Spacing (bp)',
    y='mean_monomer_length_ARMS',
    hue='Enzyme',
    s=100
)
# Error bars
#for _, row in df.iterrows():
#    plt.errorbar(
#        row['Arm Mean Spacing (bp)'],
#        row['mean_monomer_length_ARMS'],
#        xerr=row['Arm Spacing Std Dev'],
#        yerr=row['sd_monomer_length_ARMS'],
#        fmt='none', ecolor='gray', alpha=0.7, capsize=3
#    )

plt.title("Arms: Predicted Spacing vs. Experimental Monomer Length")
plt.xlabel("Arm Mean Spacing (bp)")
plt.ylabel("Mean Monomer Length (ARMS)")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend(title='Enzyme')
plt.tight_layout()
plt.show()

