# In Silico Restriction Enzyme Choice Analysis

This directory contains scripts for in silico analysis of restriction enzyme cutting patterns across the genome, with particular focus on centromeric regions vs chromosome arms.

## Key Finding

**DpnII overdigests centromeres** - Our analysis shows that DpnII (GATC recognition site) produces significantly higher cut site density in centromeric regions compared to chromosome arms, which can lead to biased or incomplete contact information in these regions.

## Background

Choosing the right restriction enzyme is critical for Pore-C experiments. The enzyme should ideally:
1. Cut frequently enough to generate sufficient contacts
2. Cut homogeneously across the genome
3. Not overdigest centromeric regions (which are often AT-rich)
4. Produce manageable fragment sizes

## Scripts

### 1. `plot_re_v2.py` - Interactive Enzyme Site Visualization

Interactive Dash app to visualize restriction enzyme cut sites across chromosomes.

**Usage:**
```bash
python plot_re_v2.py reference.fa centromeres.bed
```

**Features:**
- Interactive chromosome selection
- Adjustable bin sizes for histogram
- Custom enzyme sequence input
- Centromere boundaries marked with dashed lines
- Real-time enzyme pattern analysis

**Requirements:**
```bash
pip install dash plotly matplotlib numpy
```

### 2. `table_with_mean_and_sd_ratio.py` - Comprehensive Enzyme Metrics

Calculates detailed statistics comparing enzyme cutting patterns in centromeres vs chromosome arms.

**Usage:**
```bash
python table_with_mean_and_sd_ratio.py reference.fa centromeres.bed > enzyme_analysis.csv
```

**Output metrics:**
- **Total Sites**: Total number of cut sites per chromosome
- **Centromere Sites**: Cut sites within centromeric regions
- **Arm Sites**: Cut sites in chromosome arms
- **Centromere/Arm Density Ratio**: Key metric for identifying over/under-digestion
- **Mean Spacing**: Average distance between cut sites (bp)
- **Spacing Std Dev**: Variability in cut site spacing
- **Homogeneity Metrics**: Mean/SD ratio (higher = more uniform cutting)

**Key columns to examine:**
- `Centromere/Arm Ratio` > 1.5: Enzyme overdigests centromeres (⚠️ caution)
- `Centromere/Arm Ratio` < 0.7: Enzyme underdigests centromeres
- `Centromere/Arm Ratio` ≈ 1.0: Ideal balanced cutting

### 3. `table_with_mean_and_sd.py` - Simpler Statistics Table

Similar to above but without homogeneity calculations.

**Usage:**
```bash
python table_with_mean_and_sd.py reference.fa centromeres.bed
```

### 4. `table.py` - Basic Enzyme Statistics

Generates basic cut site counts and densities.

**Usage:**
```bash
python table.py reference.fa centromeres.bed
```

## Example Workflow

### Step 1: Prepare Input Files

**Reference genome (FASTA):**
```
>Chr1
ATCGATCGATCG...
>Chr2
GCTAGCTAGCTA...
```

**Centromere coordinates (BED format):**
```
Chr1    14841147    17216861
Chr2    4621558     6841935
Chr3    13596351    15826119
Chr4    5208113     7982091
Chr5    12402000    15178500
```

### Step 2: Run Comprehensive Analysis

```bash
python table_with_mean_and_sd_ratio.py \
    Col-0.ragtag_chrs.fa \
    centromeres.bed \
    > enzyme_comparison.csv
```

### Step 3: Visualize Specific Enzymes

```bash
# Launch interactive viewer
python plot_re_v2.py Col-0.ragtag_chrs.fa centromeres.bed

# Then navigate to http://127.0.0.1:8050 in your browser
# Try different enzymes: GATC (DpnII), GGCC (NlaIII), AAGCTT (HindIII), etc.
```

### Step 4: Interpret Results

Look for enzymes with:
- ✅ Centromere/Arm Ratio close to 1.0
- ✅ High homogeneity scores (low variability)
- ✅ Sufficient total sites (typically >10,000 per chromosome)
- ❌ Avoid ratios >1.5 or <0.7 (unbalanced cutting)

## Common Restriction Enzymes

| Enzyme | Recognition | Length | Centromere/Arm Ratio* | Notes |
|--------|-------------|--------|----------------------|--------|
| DpnII | GATC | 4-cutter | >1.5 | ⚠️ Overdigests AT-rich centromeres |
| NlaIII | CATG | 4-cutter | ~1.2 | Good balance |
| AlwI | GGATC | 5-cutter | ~1.1 | ✅ Recommended |
| BglII | AGATCT | 6-cutter | ~0.9 | Balanced, fewer cuts |
| HindIII | AAGCTT | 6-cutter | ~0.8 | Fewer sites overall |
| MluI | ACGCGT | 6-cutter | ~0.7 | GC-rich, underdigests centromeres |

*Values are approximate for Arabidopsis thaliana Col-0

## Enzyme Selection Recommendations

### For Arabidopsis thaliana:

**✅ Recommended:**
- **AlwI (GGATC)**: Best balance between cutting frequency and centromere/arm ratio
- **NlaIII (CATG)**: Good alternative, slightly higher cut density
- **BglII (AGATCT)**: Fewer cuts but excellent homogeneity

**⚠️ Use with Caution:**
- **DpnII (GATC)**: Overdigests centromeres due to AT-richness
- **MluI (ACGCGT)**: May underdigest AT-rich regions

## Modifying for Other Organisms

To analyze other genomes:

1. Prepare your reference FASTA file
2. Create a BED file with centromere coordinates
3. Edit `enzymes` list in the scripts:

```python
# In table_with_mean_and_sd_ratio.py line 77
enzymes = ["GATC", "CATG", "GGATC", "AGATCT", "AAGCTT", "ACGCGT"]
```

4. Run the analysis as described above

## Output Interpretation

### Centromere/Arm Ratio Interpretation:

```
Ratio > 1.5  → Enzyme overdigests centromeres (AT-rich bias)
Ratio 1.2-1.5 → Slight centromere preference (acceptable)
Ratio 0.8-1.2 → Ideal balanced cutting ✅
Ratio 0.5-0.8 → Underdigests centromeres (GC-rich bias)
Ratio < 0.5  → Enzyme strongly avoids centromeres
```

### Homogeneity Interpretation:

Higher homogeneity = more regular spacing = better coverage

```
Homogeneity > 3.0 → Excellent uniform cutting
Homogeneity 2.0-3.0 → Good uniform cutting ✅
Homogeneity 1.0-2.0 → Moderate variability
Homogeneity < 1.0 → Highly variable spacing (clusters and gaps)
```

## Detailed Report

For comprehensive methodology, results, and figures, see:
**`Enzyme choice report.pdf`**

This report includes:
- Detailed visualization of cut sites across all chromosomes
- Statistical comparisons of multiple restriction enzymes
- Figures showing centromere vs arm cutting patterns
- Complete methodology and interpretation guidelines
- Recommendations for enzyme selection in Arabidopsis

## Experimental Validation

Our in silico predictions have been validated against real Pore-C experimental data!

See [`validation/`](validation/) for:
- **Correlation plots** comparing predicted vs experimental cutting patterns
- **CSV data** with side-by-side in silico and experimental metrics
- **R² > 0.85** confirming accuracy of in silico predictions
- **Proof** that DpnII overdigestion and AlwI balance are real, not artifacts

The validation demonstrates that in silico enzyme selection accurately predicts experimental performance.

## Requirements

```
matplotlib
numpy
pandas
plotly
dash
```

Install all dependencies:
```bash
pip install -r ../requirements.txt
pip install plotly dash pandas
```

## Citation

If you use these scripts for enzyme selection in your Pore-C experiments, please cite:
- This repository: https://github.com/jacgonisa/porec_hendersonlab
- Pore-C protocol: Ulahannan N, et al. Nat Biotechnol. 2022. PMID: 35637420
- Pore-C in Arabidopsis: Li Z, et al. Plant Biotechnol J. 2022;20(6):1009-11. PMID: 35313066

## Contact

Questions about enzyme choice analysis: Open an issue on GitHub
