#!/bin/bash
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p YOUR_PARTITION
#SBATCH --job-name=pore_c_analysis
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00

# Example SLURM submission script for Pore-C analysis
# Adjust account, partition, and resource allocations for your HPC cluster

# Load required modules (adjust for your HPC environment)
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl  # or your system's default modules

# Set variables
ENZYME=$1       # e.g., AlwI, NlaIII, DpnII
BARCODE=$2      # e.g., 01, 02, etc.
FASTQ_DIR="/path/to/fastq/data"
REFERENCE="/path/to/reference/genome.fa"
OUTPUT_DIR="results_${ENZYME}_barcode${BARCODE}"

# Run Pore-C pipeline
nextflow run epi2me-labs/wf-pore-c \
   -profile singularity \
   --fastq ${FASTQ_DIR}/barcode${BARCODE}.fastq \
   --chunk_size 20000 \
   --cutter ${ENZYME} \
   --hi_c \
   --mcool \
   --pairs \
   --coverage \
   --ref ${REFERENCE} \
   --threads ${SLURM_CPUS_PER_TASK} \
   --out_dir ${OUTPUT_DIR}

# Optional: Convert BAM to SAM for monomer analysis
if [ -f "${OUTPUT_DIR}/bams/*.cs.bam" ]; then
    echo "Converting BAM to SAM for monomer statistics..."
    samtools view -h ${OUTPUT_DIR}/bams/*.cs.bam > ${OUTPUT_DIR}/aligned.sam

    # Generate monomer plots
    python scripts/plotting/plot_monomer_stats.py \
        ${OUTPUT_DIR}/aligned.sam \
        ${ENZYME}_barcode${BARCODE}
fi

echo "Analysis complete: ${OUTPUT_DIR}"
