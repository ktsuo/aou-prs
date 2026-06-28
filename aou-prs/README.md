This directory contains scripts used to create plots for main and supplementary figures.

## Structure

- `setup/00_setup.R`: package loading, project-relative directories, and input/output helpers.
- `setup/01_load_main_data.R`: shared loading and harmonization of PRS accuracy inputs.
- `setup/02_plot_helpers.R`: shared plotting/statistical helper functions.
- `scripts/`: scripts for each plot.

## Running scripts

From the repository root, run individual scripts, for example:

```bash
Rscript scripts/01_scatterplots_individual_phenos.R
Rscript scripts/03_heatmaps_prs_accuracy.R
```

## Script map

| Script | Purpose |
|---|---|
| `01_scatterplots_individual_phenos.R` | Quantitative and binary PRS accuracy scatterplots per individual phenotype. |
| `02_export_best_prs_metrics.R` | Exports best-performing PRS/phenotype summary tables. |
| `03_heatmaps_prs_accuracy.R` | Quantitative and binary PRS accuracy heatmaps. |
| `04_boxplots_prs_accuracy.R` | Main PRS accuracy boxplots for quantitative and binary traits. |
| `05_downsampled_prs_results.R` | Downsampled PRS boxplot and scatterplot. |
| `06_pt_vs_prscs_relative_accuracy.R` | P+T versus PRS-CS relative accuracy violin plots. |
| `07_admixture_boxplot.R` | Admixture analysis boxplot. |
| `08_odds_ratios_risk_stratification.R` | Odds-ratio/risk-stratification line and slope plots. |
| `09_greml_ldsc_heritability.R` | GREML/LDSC heritability comparison plots. |
| `10_individual_accuracy_scatter.R` | Individual-level accuracy scatterplot. |
| `11_eur_meta_neglogp_scatterplots.R` | EUR versus meta-analysis `-log10(P)` scatterplots. |