test_that("generate_did_data creates proper data structure", {
  data <- generate_did_data(n_groups = 10, n_periods = 4, seed = 123)

  # Check that data is a data frame
  expect_true(is.data.frame(data))

  # Check column names
  expected_cols <- c("group", "time", "D", "delta_D", "Y", "delta_Y", "X1", "X2")
  expect_equal(names(data), expected_cols)

  # Check dimensions (10 groups * 3 periods after removing first period)
  expect_equal(nrow(data), 10 * 3)

  # Check no NAs in delta_Y and delta_D (first period removed)
  expect_false(any(is.na(data$delta_Y)))
  expect_false(any(is.na(data$delta_D)))

  # Check D is binary
  expect_true(all(data$D %in% c(0, 1)))
})

test_that("generate_did_data respects seed", {
  data1 <- generate_did_data(n_groups = 10, n_periods = 4, seed = 456)
  data2 <- generate_did_data(n_groups = 10, n_periods = 4, seed = 456)

  # Should be identical with same seed
  expect_equal(data1, data2)
})
