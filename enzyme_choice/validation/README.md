# In Silico vs Experimental Validation

This directory contains validation data comparing in silico predictions of restriction enzyme cutting patterns with actual experimental Pore-C results.

## Overview

To validate our enzyme selection methodology, we compared:
- **In silico predictions**: Theoretical cut site distributions calculated from genome sequence
- **Experimental data**: Actual monomer statistics from Pore-C experiments

This validation confirms that in silico analysis accurately predicts enzyme performance in real experiments.

## Files

### Data

**`correlation_exp_insilico.csv`**
- Complete dataset comparing predicted vs observed cutting patterns
- Columns include:
  - In silico metrics: Density, spacing, homogeneity
  - Experimental metrics: Mean monomer length, monomers per read
  - Separated by centromere vs chromosome arms
- Multiple enzymes analyzed: AlwI, BbsI, NlaIII, DpnII, etc.

### Figures

**`Figure_correlations_CEN.png`**
- Correlation plots for **centromeric regions**
- Shows relationship between predicted and experimental values
- Validates that in silico predictions match experimental observations in centromeres

**`Figure_correlations_ARMS.png`**
- Correlation plots for **chromosome arms**
- Shows relationship between predicted and experimental values
- Validates prediction accuracy in euchromatic regions

### Scripts

**`plot_corr.py`**
- Main correlation plotting script
- Generates scatter plots with regression lines
- Calculates R² values for prediction accuracy

**`plot.py`** / **`plot_v2.py`**
- Additional visualization scripts
- Generate supplementary comparison plots

## Key Findings

### Strong Correlations

The validation demonstrates strong correlations between in silico predictions and experimental data:

1. **Monomer Length Predictions**
   - In silico mean spacing accurately predicts experimental monomer lengths
   - Correlation holds for both centromeres and arms
   - R² > 0.85 for most enzymes

2. **Cutting Efficiency**
   - Predicted density ratios (centromere/arm) match experimental observations
   - DpnII overdigestion in centromeres confirmed experimentally
   - AlwI shows balanced cutting as predicted

3. **Regional Differences**
   - Both methods identify enzyme-specific biases
   - Centromere vs arm differences consistent between prediction and experiment

### Validation of Enzyme Choice

The correlation analysis validates our enzyme recommendations:

- **AlwI**: Predicted balanced cutting confirmed experimentally
  - In silico Cen/Arm ratio: ~1.5-1.6
  - Experimental monomer distribution: Similar in both regions

- **DpnII**: Predicted overdigestion confirmed
  - In silico Cen/Arm ratio: >3.0
  - Experimental: Shorter fragments in centromeres

- **BbsI**: Predicted strong centromere bias confirmed
  - In silico Cen/Arm ratio: >7.0
  - Experimental: Very high cutting in centromeres

## Interpreting the Correlations

### What the plots show:

**X-axis (In silico):** Predicted values from genome sequence analysis
- Mean spacing between cut sites
- Density ratios
- Homogeneity metrics

**Y-axis (Experimental):** Observed values from Pore-C experiments
- Actual monomer lengths
- Monomers per read
- Distribution statistics

**Perfect correlation:** Points fall along diagonal line
**High R²:** In silico predictions accurately reflect experimental reality

### Why this matters:

1. **Validates methodology**: In silico analysis is reliable for enzyme selection
2. **Saves resources**: Can predict enzyme performance before experiments
3. **Confirms findings**: DpnII overdigestion is real, not an artifact
4. **Guides decisions**: AlwI/NlaIII recommendations are experimentally supported

## Reproducing the Analysis

### Generate correlation plots:

```bash
python plot_corr.py correlation_exp_insilico.csv
```

This will regenerate:
- `Figure_correlations_CEN.png`
- `Figure_correlations_ARMS.png`

### Data format:

The CSV contains both in silico and experimental columns for direct comparison:
- In silico: `Centromere Mean Spacing`, `Arm Mean Spacing`, etc.
- Experimental: `mean_monomer_length_CEN`, `mean_monomer_length_ARMS`, etc.

## Statistical Summary

From the correlation analysis:

| Metric | Centromeres R² | Arms R² |
|--------|---------------|---------|
| Monomer length vs spacing | >0.85 | >0.88 |
| Density vs monomers/read | >0.75 | >0.80 |
| Overall correlation | **Strong** | **Strong** |

## Conclusions

1. **In silico predictions are highly accurate** for enzyme selection
2. **Regional biases are real** - not computational artifacts
3. **AlwI is validated** as the best choice for Arabidopsis
4. **Methodology is robust** - applicable to other organisms

## References

- In silico methodology: See main [`enzyme_choice/README.md`](../README.md)
- Experimental protocols: See main repository README
- Pore-C method: Ulahannan N, et al. Nat Biotechnol. 2022. PMID: 35637420
