test_that("estim_RC_model_base works with valid input", {
  data <- generate_did_data(n_groups = 30, n_periods = 5, seed = 789)

  result <- estim_RC_model_base(
    K = 2, data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2")
  )

  expect_type(result, "list")
  expect_named(result, c("gamma", "Nobs", "B_hat", "model"))
  expect_length(result$gamma, 2)
  expect_length(result$B_hat, 3)
  expect_length(result$Nobs, 3)
  expect_equal(result$model, "base")
  expect_equal(names(result$B_hat), c("beta_0", "beta_1", "beta_2"))
})

test_that("estim_RC_model_base handles different K values", {
  data <- generate_did_data(n_groups = 30, n_periods = 6, seed = 101)

  result1 <- estim_RC_model_base(
    K = 1, data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2")
  )
  expect_length(result1$B_hat, 2)

  result3 <- estim_RC_model_base(
    K = 3, data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2")
  )
  expect_length(result3$B_hat, 4)
})

test_that("estim_RC_model_base same_sample=FALSE matches default behavior", {
  data <- generate_did_data(n_groups = 40, n_periods = 5, seed = 42)

  result_default <- estim_RC_model_base(
    K = 2, data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2")
  )
  result_false <- estim_RC_model_base(
    K = 2, data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2"), same_sample = FALSE
  )

  expect_equal(result_default$B_hat, result_false$B_hat)
  expect_equal(result_default$Nobs, result_false$Nobs)
})

test_that("estim_RC_model_base same_sample=TRUE returns n_valid_groups and uniform Nobs", {
  data <- generate_did_data(n_groups = 50, n_periods = 5, seed = 123)

  result <- estim_RC_model_base(
    K = 2, data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2"), same_sample = TRUE
  )

  expect_true(!is.null(result$n_valid_groups))
  expect_true(is.numeric(result$n_valid_groups))
  expect_true(result$n_valid_groups >= 0)

  # All Nobs should be equal (same set of groups used for every coefficient)
  expect_true(length(unique(result$Nobs)) == 1)
})

test_that("estim_RC_model_base same_sample=TRUE n_valid_groups <= total groups", {
  data <- generate_did_data(n_groups = 40, n_periods = 5, seed = 55)
  n_groups <- length(unique(data$group))

  result <- estim_RC_model_base(
    K = 2, data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2"), same_sample = TRUE
  )

  expect_true(result$n_valid_groups <= n_groups)
})
