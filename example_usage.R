# ============================================================================
# EXAMPLE USAGE OF THE dist_lag_het PACKAGE
# ============================================================================
# This script demonstrates how to use the dist_lag_het package for treatment
# effect estimation in complex designs.

# Load the package functions
source("R/estim_RC_model.R")
source("R/data_utils.R")

# ============================================================================
# STEP 1: Generate or load your data
# ============================================================================

# Option A: Generate synthetic data for testing
data <- generate_did_data(
  n_groups = 50,        # Number of groups (e.g., counties)
  n_periods = 5,        # Number of time periods
  treatment_prob = 0.3, # Probability of treatment
  seed = 123            # For reproducibility
)

# Option B: Load your own data
# data <- read.csv("your_data.csv")

# Display data structure
cat("Data structure:\n")
str(data)
cat("\nFirst few rows:\n")
print(head(data))

# ============================================================================
# STEP 2: Estimate models
# ============================================================================

# --- BASE MODEL ---
cat("\n=== ESTIMATING BASE MODEL ===\n")
result_base <- estim_RC_model_unified(
  K = 2,                    # Number of lags
  data = data,
  group_col = "group",      # Column name for group ID
  deltaY_col = "delta_Y",   # Column name for ΔY
  deltaD_col = "delta_D",   # Column name for ΔD
  D_col = "D",             # Column name for D
  X_cols = c("X1", "X2"),  # Covariate columns
  model = "base"
)

summary_RC_model(result_base)

# --- FULL DYNAMICS MODEL ---
cat("\n=== ESTIMATING FULL DYNAMICS MODEL ===\n")
result_full <- estim_RC_model_unified(
  K = NULL,                 # K is ignored for full_dynamics
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  model = "full_dynamics"
)

summary_RC_model(result_full)

# --- INTERACTIONS MODEL ---
cat("\n=== ESTIMATING INTERACTIONS MODEL ===\n")
result_int <- estim_RC_model_unified(
  K = 2,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  model = "interactions"
)

summary_RC_model(result_int)

# ============================================================================
# STEP 3: Estimate model with bootstrap standard errors
# ============================================================================

cat("\n=== ESTIMATING WITH BOOTSTRAP STANDARD ERRORS ===\n")

# Estimate the model with bootstrap standard errors
# For demonstration, we use B=100 (use B=1000 or more for real analysis)
result_with_bootstrap <- estim_RC_model_unified(
  K = 2,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  model = "base",
  bootstrap = TRUE,     # Enable bootstrap
  B = 100,              # Number of bootstrap iterations
  conf_level = 0.95     # 95% confidence intervals
)

# Display results with bootstrap statistics
summary_RC_model(result_with_bootstrap)

# ============================================================================
# STEP 4: Save results (optional)
# ============================================================================

# Save to RData file
save(
  result_base,
  result_full,
  result_int,
  result_with_bootstrap,
  file = "estimation_results.RData"
)

# Save bootstrap results to CSV (if bootstrap was computed)
if (!is.null(result_with_bootstrap$se)) {
  results_table <- data.frame(
    Coefficient = names(result_with_bootstrap$B_hat),
    Estimate = as.numeric(result_with_bootstrap$B_hat),
    Std_Error = as.numeric(result_with_bootstrap$se),
    t_stat = as.numeric(result_with_bootstrap$t_stat),
    CI_lower = as.numeric(result_with_bootstrap$ci_lower),
    CI_upper = as.numeric(result_with_bootstrap$ci_upper),
    N_obs = result_with_bootstrap$Nobs
  )
  write.csv(results_table, "estimation_results.csv", row.names = FALSE)
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to estimation_results.RData and estimation_results.csv\n")
