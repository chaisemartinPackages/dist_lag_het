test_that("summary_RC_model prints without error", {
  data <- generate_did_data(n_groups = 20, n_periods = 5, seed = 1010)

  result <- estim_RC_model_base(
    K = 2,
    data = data,
    group_col = "group",
    deltaY_col = "delta_Y",
    deltaD_col = "delta_D",
    D_col = "D",
    X_cols = c("X1", "X2")
  )

  # Should not error
  expect_output(summary_RC_model(result))

  # Should return result invisibly
  return_val <- summary_RC_model(result)
  expect_equal(return_val, result)
})
