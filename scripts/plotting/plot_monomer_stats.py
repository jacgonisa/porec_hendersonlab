import sys
import matplotlib.pyplot as plt
import numpy as np
from statistics import median
import seaborn as sns

# Use Seaborn style
sns.set(style="whitegrid")
plt.rcParams.update({'font.size': 14, 'figure.figsize': (8, 6)})

# Centromere coordinates
centromeres = {
    "Chr1": (14841147, 17216861),
    "Chr2": (4621558, 6841935),
    "Chr3": (13596351, 15826119),
    "Chr4": (5208113, 7982091),
    "Chr5": (12402000, 15178500),
}

def in_centromere(chrom, start, end):
    if chrom not in centromeres:
        return False
    cen_start, cen_end = centromeres[chrom]
    return not (end < cen_start or start > cen_end)

# Read SAM and collect monomer data
in_sam = sys.argv[1]
name = sys.argv[2]

cen_read_monomers = {}
arm_read_monomers = {}
cen_monomer_lengths = []
arm_monomer_lengths = []

with open(in_sam, "r") as f:
    for line in f:
        if line.startswith("@") or "Xw:Z:" not in line:
            continue
        fields = line.strip().split("\t")
        read_name = fields[0].split(":")[0]
        xw_field = [x for x in fields if x.startswith("Xw:Z:")][0]
        mappings = xw_field.replace("Xw:Z:", "").split(";")
        
        for m in mappings:
            if m.startswith("*"):
                continue
            parts = m.split(":")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            gcoords = parts[3]
            try:
                gstart, gend = map(int, gcoords.split("-"))
            except:
                continue
            mono_len = abs(gend - gstart)
            if in_centromere(chrom, gstart, gend):
                cen_read_monomers.setdefault(read_name, []).append(mono_len)
                cen_monomer_lengths.append(mono_len)
            else:
                arm_read_monomers.setdefault(read_name, []).append(mono_len)
                arm_monomer_lengths.append(mono_len)

# ------- Plot 1: Monomers per read -------
cen_counts = [len(v) for v in cen_read_monomers.values()]
arm_counts = [len(v) for v in arm_read_monomers.values()]
bins = np.arange(1, max(cen_counts + arm_counts) + 2)

plt.figure()
plt.hist(cen_counts, bins=bins, alpha=0.6, label=f"Centromere (median={median(cen_counts)})", color="firebrick", edgecolor="black")
plt.hist(arm_counts, bins=bins, alpha=0.6, label=f"Arms (median={median(arm_counts)})", color="steelblue", edgecolor="black")
plt.xlabel("Number of monomers per read")
plt.ylabel("Read count")
plt.title(f"{name}: Monomers per read")
plt.legend()
plt.xlim(0, 15)
plt.tight_layout()
plt.savefig("monomers_per_read_hist.png", dpi=300)
plt.close()

# ------- Plot 2: Monomer lengths -------
length_bins = np.arange(0, 2001, 25)  # bin width = 25 bp

plt.figure()
plt.hist(cen_monomer_lengths, bins=length_bins, alpha=0.6, label=f"Centromere (median={median(cen_monomer_lengths)} bp)", 
         color="firebrick", edgecolor="black")
plt.hist(arm_monomer_lengths, bins=length_bins, alpha=0.6, label=f"Arms (median={median(arm_monomer_lengths)} bp)", 
         color="steelblue", edgecolor="black")
plt.xlabel("Monomer length (bp)")
plt.ylabel("Count")
plt.title(f"{name}: Monomer lengths")
plt.legend()
plt.xlim(0, 2000)
plt.tight_layout()
plt.savefig("monomer_length_hist.png", dpi=300)
plt.close()

