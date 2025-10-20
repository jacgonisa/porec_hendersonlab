import sys
import re
import pandas as pd

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

def compute_density(chrom, cut_sites, centromere_regions):
    centromere_sites = 0
    arm_sites = 0
    
    for site in cut_sites:
        if any(start <= site <= end for start, end in centromere_regions.get(chrom, [])):
            centromere_sites += 1
        else:
            arm_sites += 1
    
    centromere_size = sum(end - start for start, end in centromere_regions.get(chrom, []))
    total_size = len(raw_data[chrom])
    arm_size = total_size - centromere_size
    
    centromere_density = (centromere_sites / centromere_size * 1000) if centromere_size > 0 else 0
    arm_density = (arm_sites / arm_size * 1000) if arm_size > 0 else 0
    
    return centromere_sites, arm_sites, centromere_density, arm_density

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
for enzyme in enzymes:
    for chrom, sequence in raw_data.items():
        cut_sites = find_sites(sequence, enzyme)
        centromere_sites, arm_sites, centromere_density, arm_density = compute_density(chrom, cut_sites, bed_data)
        
        results.append({
            "Chromosome": chrom,
            "Enzyme": enzyme,
            "Total Sites": len(cut_sites),
            "Centromere Sites": centromere_sites,
            "Arm Sites": arm_sites,
            "Centromere Density (cuts/kb)": centromere_density,
            "Arm Density (cuts/kb)": arm_density,
            "Centromere/Arm Ratio": centromere_density / arm_density if arm_density > 0 else float('inf')
        })

# Convert results to DataFrame and print
df_results = pd.DataFrame(results)
print(df_results.to_csv(index=False))

