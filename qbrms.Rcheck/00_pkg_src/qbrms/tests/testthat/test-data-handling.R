test_that("missing data handling works", {
  # Create data with missing values
  test_data <- data.frame(
    y = c(rnorm(40), rep(NA, 10)),
    x1 = rnorm(50),
    x2 = c(rnorm(45), rep(NA, 5))
  )
  
  # Test missing data removal
  clean_data <- handle_missing_data(y ~ x1 + x2, test_data)
  
  expect_true(nrow(clean_data) < nrow(test_data))
  expect_true(all(complete.cases(clean_data[c("y", "x1", "x2")])))
})

test_that("data handling preserves complete cases", {
  # Create complete data
  test_data <- data.frame(
    y = rnorm(30),
    x1 = rnorm(30),
    x2 = rnorm(30)
  )
  
  clean_data <- handle_missing_data(y ~ x1 + x2, test_data)
  expect_equal(nrow(clean_data), nrow(test_data))
})
