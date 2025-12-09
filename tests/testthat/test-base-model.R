test_that("estim_RC_model_base works with valid input", {
  # Generate test data
  data <- generate_did_data(n_groups = 30, n_periods = 5, seed = 789)

  # Estimate base model
  result <- estim_RC_model_base(
    K = 2,
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

  # Check dimensions
  expect_length(result$gamma, 2)  # 2 covariates
  expect_length(result$B_hat, 3)  # K + 1 = 3
  expect_length(result$Nobs, 3)

  # Check model type
  expect_equal(result$model, "base")

  # Check that coefficients have proper names
  expect_equal(names(result$B_hat), c("beta_0", "beta_1", "beta_2"))
})

test_that("estim_RC_model_base handles different K values", {
  data <- generate_did_data(n_groups = 30, n_periods = 6, seed = 101)

  # Test with K = 1
  result1 <- estim_RC_model_base(
    K = 1,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2")
  )
  expect_length(result1$B_hat, 2)  # K + 1 = 2

  # Test with K = 3
  result3 <- estim_RC_model_base(
    K = 3,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2")
  )
  expect_length(result3$B_hat, 4)  # K + 1 = 4
})

test_that("estim_RC_model_base works with weights", {
  data <- generate_did_data(n_groups = 20, n_periods = 5, seed = 202)

  n_groups <- length(unique(data$group))
  weights <- rep(2, n_groups)

  result <- estim_RC_model_base(
    K = 2,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2"),
    weights = weights
  )

  expect_type(result, "list")
  expect_equal(result$model, "base")
})
