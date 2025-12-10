# dist_lag_het

R package to estimate distributed lag regressions with heterogeneous treatment effects.

Based on the methodology described in Sections 4.2 and 4.3 of *"Treatment-Effect Estimation in Complex Designs under a Parallel-trends Assumption"* by Clément de Chaisemartin and Xavier D'Haultfoeuille.

## Installation

You can install the package directly from GitHub:

```r
# Install from GitHub using devtools
devtools::install_github("chaisemartinpackages/dist_lag_het")

# Or using remotes
remotes::install_github("chaisemartinpackages/dist_lag_het")


```

## Quick Start

### Estimate Models

#### Using the Unified Interface

```r
library(DistLagHet)
# Base model
result <- estim_RC_model_unified(
  K = 2,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  weights = weights,
  model = "base"
)

# View results
summary_RC_model(result)

# Estimate with bootstrap standard errors
result_boot <- estim_RC_model_unified(
  K = 2,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  model = "base",
  bootstrap = TRUE,  # Enable bootstrap
  B = 1000           # Number of bootstrap iterations
)

summary_RC_model(result_boot)
```

#### Direct Model Functions



## Data Format

Your data should be in long format with the following columns:

- **group identifier**: Column identifying groups (e.g., counties, individuals)
- **D**: Treatment level
- **delta_D**: First difference of treatment (ΔD)
- **Y**: Outcome variable
- **delta_Y**: First difference of outcome (ΔY)
- **X1, X2, ...**: Covariates

The first differences should exclude the first period (which would have NA values).

## Options

K: the number of treatment lags assumed to affect the current outcome. For instance, K=2 assumes that the current treatment and its first two lags affect the outcome.  

The package provides three model specifications:
1. **Base Model** (`model = "base"`): Standard distributed-lag model with K lags
2. **Interactions** (`model = "interactions"`): Includes interaction terms between current treatments and lags
3. **Full Dynamics** (`model = "full_dynamics"`): Allows all treatment lags up to period one to affect the outcome (then K does not need to be specified).

weights: to weight the estimation. 

## Functions

### Main Estimation Functions

- `estim_RC_model_unified()`: Unified interface for all models with optional bootstrap standard errors

### Utility Functions

- `generate_did_data()`: Generate synthetic panel data for testing
- `summary_RC_model()`: Print formatted summary of results

## Output Structure

All estimation functions return a list with:

- `gamma`: Coefficients for covariates
- `B_hat`: Treatment effect coefficients (betas)
- `Nobs`: Number of observations used for each coefficient
- `model`: Model type ("base", "full_dynamics", or "interactions")


## Dependencies

- R (>= 3.5.0)
- MASS

## License

MIT License - see LICENSE file for details

## Authors

- Clément de Chaisemartin
- Xavier D'Haultfoeuille
- Henri Fabre (package maintainer)

## Citation

If you use this package in your research, please cite:

```
de Chaisemartin, Clément and d'Haultfoeuille, Xavier, Treatment-Effect Estimation in Complex Designs under a Parallel-trends Assumption (September 04, 2025). Available at SSRN: https://ssrn.com/abstract=5440734 or http://dx.doi.org/10.2139/ssrn.5440734
```

## Support

For issues and questions:
- Open an issue on GitHub
- Contact the package maintainer : henri.fabre@ensae.fr

## References

de Chaisemartin, Clément and d'Haultfoeuille, Xavier, Treatment-Effect Estimation in Complex Designs under a Parallel-trends Assumption (September 04, 2025). Available at SSRN: https://ssrn.com/abstract=5440734 or http://dx.doi.org/10.2139/ssrn.5440734
