# Reference Hi-C Datasets

This directory documents the published Hi-C datasets used as reference/comparison data for our Pore-C experiments.

## Overview

We process publicly available Hi-C data from published studies to:
1. Compare with our newly generated Pore-C data
2. Validate analysis pipelines
3. Benchmark centromeric contact patterns across technologies

## Datasets

### 1. Sakamoto et al. - Arabidopsis Hi-C

**Publication:** Sakamoto et al. (related to centromere structure in Arabidopsis thaliana)

**Dataset Details:**
- **Accession:** DRR327470, DRR327471
- **Organism:** *Arabidopsis thaliana*
- **Technology:** Illumina Hi-C (paired-end)
- **Enzyme:** DpnII
- **Data location:** `hic_papers/Sakamoto_et_al/`

**Processed using:**
```bash
nextflow run nf-core/hic \
   -profile singularity \
   --input samplesheet_Sakamoto.csv \
   --fasta reference.fa \
   --digestion dpnii \
   --outdir DRR327470_Sakamoto
```

**Sample sheet format:**
```csv
sample,fastq_1,fastq_2
DRR327470_Sakamoto,Sakamoto_et_al/DRR327470_1.fastq.gz,Sakamoto_et_al/DRR327470_2.fastq.gz
```

### 2. Teano et al. - Arabidopsis Hi-C

**Publication:** Teano et al. (centromere/chromatin organization studies)

**Dataset Details:**
- **Accession:** SRR13739209, SRR13739213 (biological replicates)
- **Organism:** *Arabidopsis thaliana*
- **Technology:** Illumina Hi-C (paired-end)
- **Enzyme:** DpnII (likely)
- **Data location:** `hic_papers/Teano_et_al/`

**Biological replicates:**
- Replicate 1: SRR13739209
- Replicate 2: SRR13739213
- Replicate 3: (additional replicate if applicable)

**Processed using:**
```bash
# Replicate 1
nextflow run nf-core/hic \
   -profile singularity \
   --input samplesheet_Teano.csv \
   --fasta reference.fa \
   --digestion dpnii \
   --outdir Teano_rep1_SRR13739209

# Replicate 2
nextflow run nf-core/hic \
   -profile singularity \
   --input samplesheet_Teano_rep2.csv \
   --fasta reference.fa \
   --digestion dpnii \
   --outdir Teano_rep2_SRR13739213
```

**Sample sheet format:**
```csv
sample,fastq_1,fastq_2
Teano_rep,Teano_et_al/SRR13739209_1.fastq.gz,Teano_et_al/SRR13739209_2.fastq.gz
```

## Data Download

### Downloading SRA Data

To download data from NCBI SRA:

```bash
# Install SRA Toolkit
conda install -c bioconda sra-tools

# Download Sakamoto et al. data
fastq-dump --split-files --gzip DRR327470
fastq-dump --split-files --gzip DRR327471

# Download Teano et al. data
fastq-dump --split-files --gzip SRR13739209
fastq-dump --split-files --gzip SRR13739213
```

### Using prefetch for faster downloads

```bash
# Prefetch then extract
prefetch DRR327470
fasterq-dump --split-files DRR327470
gzip DRR327470_*.fastq

prefetch SRR13739209
fasterq-dump --split-files SRR13739209
gzip SRR13739209_*.fastq
```

## Analysis Workflow

### Complete Pipeline

```bash
# 1. Download reference data
fastq-dump --split-files --gzip <accession>

# 2. Create sample sheet
echo "sample,fastq_1,fastq_2" > samplesheet.csv
echo "sample_name,path/to/R1.fastq.gz,path/to/R2.fastq.gz" >> samplesheet.csv

# 3. Run nf-core/hic pipeline
nextflow run nf-core/hic \
   -profile singularity \
   --input samplesheet.csv \
   --fasta reference.fa \
   --digestion dpnii \
   --outdir results_hic

# 4. Compare with Pore-C results
# - Contact matrices (.mcool, .hic)
# - TAD boundaries
# - Compartment analysis
# - Centromere interaction patterns
```

## Comparison Metrics

When comparing Hi-C vs Pore-C data:

### 1. Contact Resolution
- Hi-C: Limited by read length (~150bp paired-end)
- Pore-C: Multi-way contacts from ultra-long reads (>10kb)

### 2. Centromeric Coverage
- Hi-C with DpnII: May show bias due to overdigestion
- Pore-C with AlwI/NlaIII: More balanced coverage

### 3. Key Comparisons
- Contact frequency matrices
- Distance-decay relationships
- Trans/cis contact ratios
- Centromere-specific interaction patterns
- Loop detection and TAD boundaries

## Output Files

### nf-core/hic Pipeline Outputs

```
results/
├── hicpro/
│   ├── matrix/          # Contact matrices
│   ├── valid_pairs/     # Filtered valid pairs
│   └── stats/           # Mapping statistics
├── multiqc/
│   └── multiqc_report.html  # Quality control report
└── hicexplorer/
    ├── *.cool          # Cool format matrices
    └── *.hic           # Juicebox format matrices
```

## Expected Results

### Quality Metrics

**Good Hi-C library indicators:**
- Valid pairs ratio: >40%
- Trans/cis ratio: <0.3 for Arabidopsis
- Mapping rate: >70%
- PCR duplicates: <20%

### Centromere Patterns

In Arabidopsis Hi-C data, you should observe:
- Strong intra-centromeric interactions
- Reduced long-range contacts from centromeres
- Clear compartmentalization
- Potential centromere clustering

## Publications

### Sakamoto et al.

**Citation:** [Add full citation when available]
- DOI: [Add DOI]
- PubMed: [Add PMID]
- Focus: Centromere structure and organization

### Teano et al.

**Citation:** [Add full citation when available]
- DOI: [Add DOI]
- PubMed: [Add PMID]
- Focus: Chromatin organization and centromere dynamics

## Notes

### Data Processing Notes

1. **Reference genome:** All datasets aligned to Col-0.ragtag (chromosomes + organelles)
2. **Restriction enzyme:** DpnII used in both publications
3. **Replicates:** Multiple biological replicates available for validation
4. **Controls:** Published data serves as positive controls for pipeline validation

### Comparison Strategy

Our Pore-C data (generated with AlwI/NlaIII) provides:
- **Higher resolution** in centromeric regions (avoiding DpnII bias)
- **Multi-way contacts** (>2 loci per molecule)
- **Long-range structure** (ultra-long reads)

The published Hi-C data (DpnII) provides:
- **Established baseline** for Arabidopsis chromatin structure
- **Validation dataset** for pipeline correctness
- **Comparison point** for centromere digestion bias

## Reproducing the Analysis

### Step-by-step Guide

```bash
# 1. Navigate to hic_papers directory
cd hic_papers

# 2. Download data (if not already present)
bash scripts/download_data.sh

# 3. Run Hi-C pipeline for Sakamoto data
bash scripts/run_hic_nextflow_Sakamoto.sh

# 4. Run Hi-C pipeline for Teano data
bash scripts/run_hic_nextflow_Teano.sh

# 5. Compare outputs
# - View contact matrices in Juicebox
# - Compare with Pore-C matrices
# - Analyze centromere interaction differences
```

## Related Files

- Main Hi-C pipeline script: `../scripts/hic/run_hic_nextflow.sh`
- Sample sheets: `hic_papers/samplesheet_*.csv`
- Download scripts: `hic_papers/scripts/download_data.sh`

## Contact

For questions about the reference datasets or comparison analyses, please open an issue on GitHub.

## Data Availability

All datasets are publicly available from:
- **NCBI SRA:** https://www.ncbi.nlm.nih.gov/sra
- **DDBJ:** https://www.ddbj.nig.ac.jp/ (for DRR accessions)
- **ENA:** https://www.ebi.ac.uk/ena (European mirror)

## License

The reference datasets are subject to their original publication licenses. This documentation and analysis scripts are covered under the repository's MIT License.
