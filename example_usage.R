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
# STEP 3: Bootstrap for standard errors (optional)
# ============================================================================

cat("\n=== COMPUTING BOOTSTRAP STANDARD ERRORS ===\n")

# Number of bootstrap iterations
B <- 100  # Use 1000 for real analysis

# Get unique groups
list_groups <- unique(data$group)
n_groups <- length(list_groups)

# Storage for bootstrap results
n_coef <- length(result_base$B_hat)
Bhat_boot <- matrix(0, B, n_coef)

# Progress indicator
pb <- txtProgressBar(min = 0, max = B, style = 3)

for (b in 1:B) {
  # Sample groups with replacement
  groups_boot <- sample(list_groups, n_groups, replace = TRUE)

  # Calculate weights (how many times each group appears)
  weights <- as.numeric(table(factor(groups_boot, levels = list_groups)))

  # Estimate model on bootstrap sample
  tryCatch({
    result_boot <- estim_RC_model_base(
      K = 2,
      data = data,
      group_col = "group",
      deltaY_col = "delta_Y",
      deltaD_col = "delta_D",
      D_col = "D",
      X_cols = c("X1", "X2"),
      weights = weights
    )
    Bhat_boot[b, ] <- result_boot$B_hat
  }, error = function(e) {
    Bhat_boot[b, ] <- NA
  })

  setTxtProgressBar(pb, b)
}
close(pb)

# Remove failed iterations
Bhat_boot <- Bhat_boot[complete.cases(Bhat_boot), , drop = FALSE]

# Compute standard errors
se_boot <- apply(Bhat_boot, 2, sd)

# Display results
cat("\n\n=== BOOTSTRAP RESULTS ===\n")
results_table <- data.frame(
  Coefficient = names(result_base$B_hat),
  Estimate = as.numeric(result_base$B_hat),
  Std_Error = se_boot,
  t_stat = as.numeric(result_base$B_hat) / se_boot,
  CI_lower = as.numeric(result_base$B_hat) - 1.96 * se_boot,
  CI_upper = as.numeric(result_base$B_hat) + 1.96 * se_boot,
  N_obs = result_base$Nobs
)
print(results_table, row.names = FALSE)

# ============================================================================
# STEP 4: Save results (optional)
# ============================================================================

# Save to RData file
save(
  result_base,
  result_full,
  result_int,
  results_table,
  file = "estimation_results.RData"
)

# Save to CSV
write.csv(results_table, "estimation_results.csv", row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to estimation_results.RData and estimation_results.csv\n")
