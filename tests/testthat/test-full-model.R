test_that("estim_RC_model_full works with valid input", {
  data <- generate_did_data(n_groups = 30, n_periods = 5, seed = 303)

  result <- estim_RC_model_full(
    data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2")
  )

  expect_type(result, "list")
  expect_named(result, c("gamma", "Nobs", "B_hat", "model"))
  expect_equal(result$model, "full_dynamics")
  expect_length(result$gamma, 2)
  # generate_did_data removes the t=1 NA row, so T_periods = n_periods - 1 = 4
  # dim = T_periods - 1 = 3 coefficients (code expects data WITH the t=1 NA row kept)
  expect_length(result$B_hat, 3)
})

test_that("estim_RC_model_full same_sample=FALSE matches default behavior", {
  data <- generate_did_data(n_groups = 30, n_periods = 5, seed = 42)

  result_default <- estim_RC_model_full(
    data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2")
  )
  result_false <- estim_RC_model_full(
    data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2"), same_sample = FALSE
  )

  expect_equal(result_default$B_hat, result_false$B_hat)
  expect_equal(result_default$Nobs, result_false$Nobs)
})

test_that("estim_RC_model_full same_sample=TRUE returns n_valid_groups and uniform Nobs", {
  data <- generate_did_data(n_groups = 50, n_periods = 5, seed = 77)

  result <- estim_RC_model_full(
    data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2"), same_sample = TRUE
  )

  expect_true(!is.null(result$n_valid_groups))
  expect_true(is.numeric(result$n_valid_groups))
  expect_true(length(unique(result$Nobs)) == 1)
})

test_that("estim_RC_model_full same_sample=TRUE n_valid_groups <= total groups", {
  data <- generate_did_data(n_groups = 40, n_periods = 5, seed = 88)
  n_groups <- length(unique(data$group))

  result <- estim_RC_model_full(
    data = data, group_col = "group",
    deltaY_col = "delta_Y", deltaD_col = "delta_D",
    D_col = "D", X_cols = c("X1", "X2"), same_sample = TRUE
  )

  expect_true(result$n_valid_groups <= n_groups)
})
