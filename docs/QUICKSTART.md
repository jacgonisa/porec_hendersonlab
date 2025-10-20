# Quick Start Guide

## 1. Clone the Repository

```bash
git clone https://github.com/jacgonisa/porec_hendersonlab.git
cd porec_hendersonlab
```

## 2. Install Python Dependencies

```bash
pip install -r requirements.txt
```

Or with conda:
```bash
conda install matplotlib numpy seaborn
```

## 3. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

## 4. Run Your First Analysis

### Pore-C Quick Test

```bash
# Prepare your data
fastq_file="your_reads.fastq"
reference="your_reference.fa"
enzyme="NlaIII"

# Run the pipeline
nextflow run epi2me-labs/wf-pore-c \
   -profile singularity \
   --fastq $fastq_file \
   --cutter $enzyme \
   --ref $reference \
   --out_dir results_test
```

### On SLURM Cluster

```bash
# Edit the script with your parameters
nano scripts/porec/run_porec_analysis.sh

# Submit the job
sbatch scripts/porec/sub_porec.sh AlwI 01
```

## 5. Analyze Results

After the pipeline completes:

```bash
# Convert BAM to SAM
samtools view -h results_test/bams/sample.cs.bam > sample.cs.sam

# Generate plots
python scripts/plotting/plot_monomer_stats.py sample.cs.sam test_analysis
```

This creates:
- `monomers_per_read_hist.png`
- `monomer_length_hist.png`

## Common Workflows

### Full Pore-C Workflow with QC

```bash
# 1. QC raw reads
NanoPlot --fastq raw_reads.fastq -o qc_raw/

# 2. Filter reads (optional)
filtlong --min_length 1000 --min_mean_q 10 \
   raw_reads.fastq > filtered_reads.fastq

# 3. Run Pore-C
nextflow run epi2me-labs/wf-pore-c \
   -profile singularity \
   --fastq filtered_reads.fastq \
   --cutter AlwI \
   --ref reference.fa \
   --out_dir results_AlwI \
   --hi_c --mcool --pairs

# 4. Plot statistics
samtools view -h results_AlwI/bams/*.cs.bam > aligned.sam
python scripts/plotting/plot_monomer_stats.py aligned.sam AlwI_sample
```

### Hi-C Workflow

```bash
# Prepare samplesheet (CSV format):
# sample,fastq_1,fastq_2
# sample1,reads_R1.fastq.gz,reads_R2.fastq.gz

nextflow run nf-core/hic \
   -profile singularity \
   --input samplesheet.csv \
   --fasta reference.fa \
   --digestion dpnii \
   --outdir results_hic
```

## Next Steps

- See [README.md](../README.md) for detailed documentation
- Check [SLURM_GUIDE.md](SLURM_GUIDE.md) for cluster-specific info
- Browse `examples/` for more submission scripts

## Troubleshooting

**Can't find Nextflow pipelines?**
```bash
# Pull pipelines manually
nextflow pull epi2me-labs/wf-pore-c
nextflow pull nf-core/hic
```

**Permission denied on scripts?**
```bash
chmod +x scripts/porec/*.sh
chmod +x scripts/hic/*.sh
```

**Python module not found?**
```bash
# Install with pip
pip install matplotlib numpy seaborn

# Or use conda
conda install -c conda-forge matplotlib numpy seaborn
```
