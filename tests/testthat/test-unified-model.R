test_that("estim_RC_model_unified dispatches to correct model", {
  data <- generate_did_data(n_groups = 30, n_periods = 5, seed = 707)

  # Test base model
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
  expect_equal(result_base$model, "base")

  # Test full_dynamics model
  result_full <- estim_RC_model_unified(
    K = NULL,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2"),
    model = "full_dynamics"
  )
  expect_equal(result_full$model, "full_dynamics")

  # Test interactions model
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
  expect_equal(result_int$model, "interactions")
})

test_that("estim_RC_model_unified validates K parameter", {
  data <- generate_did_data(n_groups = 20, n_periods = 5, seed = 808)

  # Should error when K is NULL for base model
  expect_error(
    estim_RC_model_unified(
      K = NULL,
      data = data,
      group_col = "group",
      deltaY_col = "delta_Y",
      deltaD_col = "delta_D",
      D_col = "D",
      X_cols = c("X1", "X2"),
      model = "base"
    ),
    "K must be specified"
  )

  # Should error when K is NULL for interactions model
  expect_error(
    estim_RC_model_unified(
      K = NULL,
      data = data,
      group_col = "group",
      deltaY_col = "delta_Y",
      deltaD_col = "delta_D",
      D_col = "D",
      X_cols = c("X1", "X2"),
      model = "interactions"
    ),
    "K must be specified"
  )

  # Should warn when K is provided for full_dynamics model
  expect_warning(
    estim_RC_model_unified(
      K = 2,
      data = data,
      group_col = "group",
      deltaY_col = "delta_Y",
      deltaD_col = "delta_D",
      D_col = "D",
      X_cols = c("X1", "X2"),
      model = "full_dynamics"
    ),
    "K is ignored"
  )
})

test_that("estim_RC_model_unified accepts default model argument", {
  data <- generate_did_data(n_groups = 20, n_periods = 5, seed = 909)

  # Default should be "base"
  result <- estim_RC_model_unified(
    K = 2,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2")
  )

  expect_equal(result$model, "base")
})
