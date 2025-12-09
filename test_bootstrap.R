# Quick test of the integrated bootstrap functionality

# Load the package functions
source("R/estim_RC_model.R")
source("R/data_utils.R")

# Generate test data
set.seed(123)
data <- generate_did_data(
  n_groups = 30,
  n_periods = 5,
  treatment_prob = 0.3
)

cat("=== Testing Base Model (without bootstrap) ===\n")
# Estimate base model without bootstrap
result_base <- estim_RC_model_unified(
  K = 2,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  model = "base"
)

summary_RC_model(result_base)

cat("\n=== Testing Base Model WITH Bootstrap (10 iterations) ===\n")
# Test bootstrap with small number of iterations
result_boot <- estim_RC_model_unified(
  K = 2,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  model = "base",
  bootstrap = TRUE,
  B = 10,
  verbose = TRUE
)

summary_RC_model(result_boot)

cat("\n=== Testing Full Dynamics Model WITH Bootstrap ===\n")
result_full_boot <- estim_RC_model_unified(
  K = NULL,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  model = "full_dynamics",
  bootstrap = TRUE,
  B = 10,
  verbose = TRUE
)

summary_RC_model(result_full_boot)

cat("\n=== Testing Interactions Model WITH Bootstrap ===\n")
result_int_boot <- estim_RC_model_unified(
  K = 2,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  model = "interactions",
  bootstrap = TRUE,
  B = 10,
  verbose = TRUE
)

summary_RC_model(result_int_boot)

cat("\n=== All Tests Complete ===\n")
