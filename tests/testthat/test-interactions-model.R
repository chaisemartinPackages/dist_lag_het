test_that("estim_RC_model_interactions works with valid input", {
  # Generate test data
  data <- generate_did_data(n_groups = 30, n_periods = 5, seed = 505)

  # Estimate interactions model
  result <- estim_RC_model_interactions(
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

  # Check model type
  expect_equal(result$model, "interactions")

  # Check dimensions
  expect_length(result$gamma, 2)  # 2 covariates

  # Number of parameters = (K+1)*(K+2)/2 = 3*4/2 = 6
  expect_length(result$B_hat, 6)

  # Check coefficient names include both main effects and interactions
  expect_true(any(grepl("beta_0,0", names(result$B_hat))))
  expect_true(any(grepl("beta_0,1", names(result$B_hat))))
})

test_that("estim_RC_model_interactions handles different K values", {
  data <- generate_did_data(n_groups = 30, n_periods = 6, seed = 606)

  # Test with K = 1: (1+1)*(1+2)/2 = 3 parameters
  result1 <- estim_RC_model_interactions(
    K = 1,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2")
  )
  expect_length(result1$B_hat, 3)

  # Test with K = 2: (2+1)*(2+2)/2 = 6 parameters
  result2 <- estim_RC_model_interactions(
    K = 2,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2")
  )
  expect_length(result2$B_hat, 6)
})
