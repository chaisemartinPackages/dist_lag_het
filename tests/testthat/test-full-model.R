test_that("estim_RC_model_full works with valid input", {
  # Generate test data with balanced panel
  data <- generate_did_data(n_groups = 30, n_periods = 5, seed = 303)

  # Estimate full dynamics model
  result <- estim_RC_model_full(
    K = NULL,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2")
  )

  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("gamma", "Nobs", "B_hat", "model"))

  # Check model type
  expect_equal(result$model, "full_dynamics")

  # Check dimensions
  expect_length(result$gamma, 2)  # 2 covariates

  # For full dynamics, number of coefficients = T - 1 where T is max periods
  # After removing first period for differences, we have 4 periods
  expect_length(result$B_hat, 4)
})

test_that("estim_RC_model_full handles never-switchers", {
  # Generate data where some groups never switch
  data <- generate_did_data(n_groups = 40, n_periods = 5, treatment_prob = 0.5, seed = 404)

  result <- estim_RC_model_full(
    K = NULL,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2")
  )

  expect_type(result, "list")
  expect_equal(result$model, "full_dynamics")
})
