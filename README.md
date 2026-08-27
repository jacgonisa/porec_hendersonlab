# PoreC Analysis Pipeline

> **A joint contribution by Hanwen Cao** (hc617@cam.ac.uk), **[Jacob Gonzalez Isa](https://github.com/jacgonisa)** (jg2070@cam.ac.uk), **and [Katie Jenike](https://github.com/kjenike)** (kj436@cam.ac.uk) — with special thanks to Katie Jenike for the design of the cut-density track.

A 3C genomics project to understand centromeres better through Hi-C and Pore-C chromatin conformation capture analysis at Henderson Lab.

## Overview

This repository contains workflows for:
- **Restriction enzyme selection** - In silico analysis to choose optimal enzymes (avoiding DpnII's centromere overdigestion). Also available as a no-install browser app: [`webapp/`](webapp/) — upload a genome FASTA + centromere BED and get the enzyme ranking, metrics, and plots (hostable on GitHub Pages).
- **Hi-C analysis** using nf-core/hic pipeline
- **Pore-C analysis** using epi2me-labs/wf-pore-c pipeline
- **Monomer statistics** plotting and visualization for centromeric regions

## Repository Structure

```
.
├── enzyme_choice/            # In silico enzyme selection analysis (with PDF report)
├── webapp/                   # Browser app for enzyme selection (FASTA + BED → plots)
├── reference_data/           # Documentation for published Hi-C datasets
├── scripts/
│   ├── hic/                  # Hi-C pipeline scripts
│   ├── porec/                # Pore-C pipeline scripts
│   ├── plotting/             # Data visualization scripts
│   └── utils/                # Utility scripts (filtlong, nanoplot, etc.)
├── examples/
│   ├── plots/                # Example monomer analysis plots
│   └── slurm_submission_example.slurm
├── config/                   # Configuration files
└── docs/                     # Additional documentation
```

## Requirements

### Software Dependencies
- [Nextflow](https://www.nextflow.io/) (>=21.04.0)
- [Singularity](https://sylabs.io/singularity/) or Docker
- Python 3.7+ with:
  - matplotlib
  - numpy
  - seaborn
- [Samtools](http://www.htslib.org/)
- [Filtlong](https://github.com/rrwick/Filtlong) (optional, for read filtering)
- [NanoPlot](https://github.com/wdecoster/NanoPlot) (optional, for QC)

### HPC Environment
These scripts are designed for SLURM-based HPC clusters. Adjust resource allocations in submission scripts as needed.

## Usage

### 0. Enzyme Selection (Recommended First Step)

Before running experiments, analyze which restriction enzyme is best suited for your genome:

```bash
# Generate comprehensive enzyme statistics
python enzyme_choice/table_with_mean_and_sd_ratio.py \
    reference.fa \
    centromeres.bed \
    > enzyme_comparison.csv

# Or use the interactive visualizer
python enzyme_choice/plot_re_v2.py reference.fa centromeres.bed
```

**Key finding for Arabidopsis:** DpnII (GATC) overdigests centromeric regions (Centromere/Arm ratio >1.5). We recommend **AlwI** or **NlaIII** instead for more balanced cutting.

See [`enzyme_choice/README.md`](enzyme_choice/README.md) for detailed guidance.

### 1. Hi-C Analysis

Run the nf-core/hic pipeline with paired-end short reads:

```bash
nextflow run nf-core/hic \
   -profile singularity \
   --input samplesheet.csv \
   --fasta reference.fa \
   --outdir results_hic \
   --digestion "dpnii"
```

See `scripts/hic/run_hic_nextflow.sh` for a working example.

### 2. Pore-C Analysis

Run the epi2me-labs/wf-pore-c pipeline with Oxford Nanopore long reads:

```bash
nextflow run epi2me-labs/wf-pore-c \
   -profile singularity \
   --fastq reads.fastq \
   --chunk_size 20000 \
   --cutter NlaIII \
   --hi_c \
   --mcool \
   --pairs \
   --coverage \
   --ref reference.fa \
   --threads 19 \
   --out_dir results_porec
```

**Common restriction enzymes:**
- `NlaIII` - 4-cutter (CATG^)
- `DpnII` - 4-cutter (^GATC)
- `AlwI` - 4-cutter (GGATC)
- `BglII` - 6-cutter (A^GATCT)
- `MluI` - 6-cutter (A^CGCGT)

See `scripts/porec/run_porec_analysis.sh` for a working example.

#### SLURM Submission

For cluster job submission:

```bash
sbatch scripts/porec/sub_porec.sh <enzyme> <barcode>
```

Example:
```bash
sbatch scripts/porec/sub_porec.sh AlwI 01
```

### 3. Monomer Statistics Plotting

After Pore-C analysis, visualize monomers per read and monomer length distributions:

**Step 1: Convert BAM to SAM**
```bash
samtools view -h results_porec/bams/sample.cs.bam > sample.cs.sam
```

**Step 2: Generate plots**
```bash
python scripts/plotting/plot_monomer_stats.py sample.cs.sam output_prefix
```

This generates two plots:
- `monomers_per_read_hist.png` - Distribution of contact monomers per read
- `monomer_length_hist.png` - Distribution of monomer lengths

**Example plots:** See [`examples/plots/`](examples/plots/) for example outputs

The script automatically separates statistics for:
- **Centromeric regions** (firebrick)
- **Chromosome arms** (steelblue)

**Centromere coordinates (Arabidopsis thaliana Col-0):**
- Chr1: 14,841,147 - 17,216,861
- Chr2: 4,621,558 - 6,841,935
- Chr3: 13,596,351 - 15,826,119
- Chr4: 5,208,113 - 7,982,091
- Chr5: 12,402,000 - 15,178,500

## Utilities

### Quality Filtering with Filtlong

Filter reads by quality and length:

```bash
filtlong --min_length 1000 --min_mean_q 10 input.fastq > output.filtered.fastq
```

See `scripts/utils/filtlong.sh`

### QC with NanoPlot

Generate quality control plots:

```bash
NanoPlot --fastq input.fastq -o nanoplot_output/
```

See `scripts/utils/nanoplot.sh`

## Reference Hi-C Data

We process published Hi-C datasets for comparison with our Pore-C data:

- **Sakamoto et al.** - Arabidopsis Hi-C (DRR327470, DRR327471)
- **Teano et al.** - Arabidopsis Hi-C replicates (SRR13739209, SRR13739213)

These datasets serve as:
- Validation controls for pipeline correctness
- Comparison baseline for centromere interaction patterns
- Reference for evaluating DpnII overdigestion vs our AlwI/NlaIII approach

See [`reference_data/README.md`](reference_data/README.md) for detailed information on downloading and processing these datasets.

## Example Workflows

### Complete Pore-C Analysis Pipeline

```bash
# 1. Optional: Filter reads
filtlong --min_length 1000 --min_mean_q 10 raw_reads.fastq > filtered_reads.fastq

# 2. Optional: QC
NanoPlot --fastq filtered_reads.fastq -o qc_results/

# 3. Run Pore-C pipeline
nextflow run epi2me-labs/wf-pore-c \
   -profile singularity \
   --fastq filtered_reads.fastq \
   --cutter AlwI \
   --ref reference.fa \
   --out_dir results_AlwI

# 4. Convert BAM to SAM
samtools view -h results_AlwI/bams/sample.cs.bam > sample.cs.sam

# 5. Generate monomer statistics plots
python scripts/plotting/plot_monomer_stats.py sample.cs.sam AlwI_analysis
```

## Output Files

### Pore-C Pipeline Outputs
- `bams/` - Aligned reads with contact annotations
- `pairs/` - Contact pairs in .pairs format
- `mcool/` - Multi-resolution cooler files for visualization
- `hic/` - .hic format files (if `--hi_c` enabled)
- `coverage/` - Genomic coverage tracks

### Monomer Plots
- `monomers_per_read_hist.png` - Histogram showing distribution of monomers detected per read
- `monomer_length_hist.png` - Histogram showing distribution of individual monomer lengths

## Configuration Notes

### Memory and CPU Requirements
- **Hi-C pipeline**: 64 GB RAM, 32 CPUs recommended
- **Pore-C pipeline**: 128 GB RAM, 64 CPUs recommended (adjust `--chunk_size` if memory limited)

### Nextflow Profiles
- `singularity` - Uses Singularity containers (recommended for HPC)
- `docker` - Uses Docker containers
- `standard` - Local execution without containers

## Troubleshooting

**Issue: Nextflow fails with "Cannot find revision"**
- Solution: Ensure internet connection for pipeline download, or use `-offline` with cached pipelines

**Issue: Out of memory errors in Pore-C**
- Solution: Reduce `--chunk_size` parameter (default: 20000)

**Issue: No monomers detected in plots**
- Solution: Check that BAM file contains `Xw:Z:` tags (contact walks)

## References

**Pipelines:**
- [nf-core/hic](https://github.com/nf-core/hic)
- [epi2me-labs/wf-pore-c](https://github.com/epi2me-labs/wf-pore-c)

**Publications:**
- Pore-C protocol: Ulahannan N, et al. Nat Biotechnol. 2022. PMID: [35637420](https://pubmed.ncbi.nlm.nih.gov/35637420/)
- Pore-C in Arabidopsis: Li Z, et al. Plant Biotechnol J. 2022;20(6):1009-11. PMID: [35313066](https://pubmed.ncbi.nlm.nih.gov/35313066/)

## License

MIT License

## Contact

Henderson Lab
Questions: Open an issue on GitHub
