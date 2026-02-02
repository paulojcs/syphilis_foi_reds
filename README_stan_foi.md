# Syphilis Force of Infection Analysis - Custom Stan Models

## Overview

This package provides custom Stan models for estimating the Force of Infection (FOI) of syphilis from age-stratified seroprevalence data in Brazilian blood donors.

## Files Included

### Stan Models (`stan_models/`)

1. **`foi_constant.stan`** - Constant FOI model (baseline)
   - Single FOI parameter for entire period
   - Simplest model, serves as null hypothesis

2. **`foi_piecewise.stan`** - Piecewise constant FOI
   - K periods with different FOI values
   - Flexible: define periods via `foi_index` input

3. **`foi_random_walk.stan`** - Smooth random walk FOI
   - Year-by-year FOI with smoothness constraint
   - Most flexible temporal structure

4. **`foi_random_walk_v2.stan`** - Optimized RW model
   - Non-centered parameterization for better MCMC
   - Optional seroreversion support

### R Scripts

1. **`syphilis_foi_stan_analysis.R`** - Full analysis script
   - Loads processed data from previous serofoi analysis
   - Fits all 4 models
   - Model comparison via LOO-CV
   - Publication-quality figures

2. **`syphilis_foi_stan_minimal.R`** - Standalone test version
   - Self-contained with inline Stan code
   - Works with simulated data for testing
   - Quick validation of model structure

### HPC Scripts

- **`foi_stan.slurm`** - SLURM batch submission script for USP aguia4 cluster

## Quick Start

### Option 1: HPC (Full Analysis)

```bash
# Copy files to HPC
scp stan_models/*.stan username@aguia4.hpc.usp.br:/scratch/7631403/reds_data/stan_models/
scp syphilis_foi_stan_analysis.R username@aguia4.hpc.usp.br:/scratch/7631403/reds_data/
scp foi_stan.slurm username@aguia4.hpc.usp.br:/scratch/7631403/reds_data/

# Submit job
ssh username@aguia4.hpc.usp.br
cd /scratch/7631403/reds_data
sbatch foi_stan.slurm
```

### Critical Notes

1. **Year Range**: FOI is estimated for years `min_birth_year` to `survey_year - 1`
   - For survey_year = 2023 and max_age = 65:
   - FOI years = 1958 to 2022 (65 years)

2. **Age Groups**: Use single-year ages (not grouped)

3. **First Donation Only**: Use first donation per donor for statistical independence

## Model Specifications

### Piecewise-3 Periods (REDS Eras)

- Period 1 (1958-2011): Historical + REDS-II (54 years)
- Period 2 (2012-2016): REDS-III (5 years)
- Period 3 (2017-2022): REDS-IV (6 years)

### Piecewise-5 Periods (Finer Resolution)

- Period 1 (1958-2010): Pre-observation baseline
- Period 2 (2011-2013): Early observation
- Period 3 (2014-2017): Mid observation
- Period 4 (2018-2020): Pre-COVID + COVID
- Period 5 (2021-2022): Recent

## Outputs

### Data Files

- `foi_estimates_constant.csv`
- `foi_estimates_piecewise_3.csv`
- `foi_estimates_piecewise_5.csv`
- `foi_estimates_random_walk.csv`
- `model_comparison_loo.csv`

### Figures

- `fig1_foi_model_comparison.png/pdf` - All models overlaid
- `fig2_posterior_predictive.png/pdf` - Observed vs predicted
- `fig3_foi_by_period.png/pdf` - Period-specific FOI
- `fig4_mcmc_diagnostics.png` - Trace plots

### R Workspace

- `syphilis_foi_stan_analysis.RData` - All fitted models and results

## Troubleshooting

### Stan Compilation Errors

```r
# On Windows, try:
rstan_options(auto_write = TRUE)
Sys.setenv(USE_CXX14 = 1)
```

### Convergence Issues

```r
# Increase adapt_delta
control = list(adapt_delta = 0.99, max_treedepth = 15)

# Increase iterations
iter = 10000, warmup = 5000
```

### Memory Issues

```r
# Use aggregated data (critical)
# Individual-level data (2M rows) is infeasible
# Age-year-aggregated (48x16 = 768 rows) works well
```

## Interpretation Guidelines

### Period 1 (Historical) Estimates

- Represents FOI over 54+ years
- Includes decades before direct observation

### Selection Bias

- Blood donors ≠ general population
- Permanent deferral of seropositives creates "survivor" effect
- FOI estimates may appear lower than population-level truth

### Model Comparison

- ΔELPD < 2 SE: Models equivalent, prefer simpler
- ΔELPD 2-5 SE: Weak preference
- ΔELPD > 5 SE: Strong preference for better model

## References

- Hozé et al. (2025). RSero: Reconstructing pathogen circulation from seroprevalence. _PLOS Comp Biol_.
- serofoi package: https://epiverse-trace.github.io/serofoi/
- Stan User's Guide: https://mc-stan.org/users/documentation/

## Author

Paulo - Research on syphilis epidemiology in Brazil
University of São Paulo (USP)
