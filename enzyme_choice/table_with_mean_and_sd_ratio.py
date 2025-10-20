import sys
import re
import pandas as pd
import numpy as np

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse(dna):
    return "".join(complement.get(base, base) for base in reversed(dna))

def read_sequence(fasta_f):
    raw_data = {}
    with open(fasta_f, "r") as f:
        line = f.readline()
        while line:
            if line.startswith(">"):
                this_name = line.strip().split(">")[1]
                raw_data[this_name] = ""
            else:
                raw_data[this_name] += line.strip()
            line = f.readline()
    return raw_data

def read_bed(bed_file):
    bed_data = {}
    with open(bed_file, "r") as f:
        for line in f:
            chrom, start, end = line.strip().split("\t")
            if chrom not in bed_data:
                bed_data[chrom] = []
            bed_data[chrom].append((int(start), int(end)))
    return bed_data

def find_sites(sequence, enzyme):
    r_enzyme = reverse(enzyme)
    cut_sites = [m.start() for m in re.finditer(enzyme, sequence)]
    cut_sites += [m.start() for m in re.finditer(r_enzyme, sequence)]
    return sorted(cut_sites)

def compute_density_and_spacing(chrom, cut_sites, centromere_regions):
    centromere_sites = []
    arm_sites = []
    
    for site in cut_sites:
        if any(start <= site <= end for start, end in centromere_regions.get(chrom, [])):
            centromere_sites.append(site)
        else:
            arm_sites.append(site)
    
    centromere_size = sum(end - start for start, end in centromere_regions.get(chrom, []))
    total_size = len(raw_data[chrom])
    arm_size = total_size - centromere_size
    
    centromere_density = (len(centromere_sites) / centromere_size * 1000) if centromere_size > 0 else 0
    arm_density = (len(arm_sites) / arm_size * 1000) if arm_size > 0 else 0
    
    centromere_spacings = np.diff(centromere_sites) if len(centromere_sites) > 1 else []
    arm_spacings = np.diff(arm_sites) if len(arm_sites) > 1 else []
    
    centromere_mean_spacing = np.mean(centromere_spacings) if len(centromere_spacings) > 0 else 0
    centromere_std_spacing = np.std(centromere_spacings) if len(centromere_spacings) > 0 else 0
    arm_mean_spacing = np.mean(arm_spacings) if len(arm_spacings) > 0 else 0
    arm_std_spacing = np.std(arm_spacings) if len(arm_spacings) > 0 else 0
    
    return (len(centromere_sites), len(arm_sites), centromere_density, arm_density,
            centromere_mean_spacing, centromere_std_spacing, arm_mean_spacing, arm_std_spacing)

# Input files
fasta_f = sys.argv[1]
bed_file = sys.argv[2]

# Read data
raw_data = read_sequence(fasta_f)
bed_data = read_bed(bed_file)

# Enzymes to analyze
enzymes = ["GATC", "GGCC", "AAGCTT", "CTGCAG"]  # Add more as needed

# Collect results
results = []
epsilon = 1e-6  # small constant to avoid division by zero

for enzyme in enzymes:
    for chrom, sequence in raw_data.items():
        cut_sites = find_sites(sequence, enzyme)
        (centromere_sites,
         arm_sites,
         centromere_density,
         arm_density,
         centromere_mean_spacing,
         centromere_std_spacing,
         arm_mean_spacing,
         arm_std_spacing) = compute_density_and_spacing(chrom, cut_sites, bed_data)
        
        # Compute homogeneity metrics
        if centromere_mean_spacing > 0:
            cen_homogeneity = centromere_mean_spacing / (centromere_std_spacing + epsilon)
        else:
            cen_homogeneity = np.nan
        
        if arm_mean_spacing > 0:
            arm_homogeneity = arm_mean_spacing / (arm_std_spacing + epsilon)
        else:
            arm_homogeneity = np.nan
        
        # Difference between centromere and arm homogeneity
        homogeneity_difference = cen_homogeneity - arm_homogeneity
        
        results.append({
            "Chromosome": chrom,
            "Enzyme": enzyme,
            "Total Sites": len(cut_sites),
            "Centromere Sites": centromere_sites,
            "Arm Sites": arm_sites,
            "Centromere Density (cuts/kb)": centromere_density,
            "Arm Density (cuts/kb)": arm_density,
            "Centromere/Arm Ratio": centromere_density / arm_density if arm_density > 0 else float('inf'),
            "Centromere Mean Spacing (bp)": centromere_mean_spacing,
            "Centromere Spacing Std Dev": centromere_std_spacing,
            "Arm Mean Spacing (bp)": arm_mean_spacing,
            "Arm Spacing Std Dev": arm_std_spacing,
            "Centromere Homogeneity": cen_homogeneity,
            "Arm Homogeneity": arm_homogeneity,
            "Homogeneity Difference (Cen - Arm)": homogeneity_difference
        })

# Convert results to DataFrame and print
df_results = pd.DataFrame(results)
print(df_results.to_csv(index=False))

