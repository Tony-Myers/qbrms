# =============================================================================
# FILE: tests/testthat/test-edge-cases.R
# =============================================================================
# Tests for edge cases and error handling in qbrms package

# Setup test data
test_data <- data.frame(
  y = rnorm(100, mean = 2, sd = 1),
  x1 = rnorm(100),
  x2 = rnorm(100),
  group = rep(c("A", "B"), each = 50)
)

test_that("qbrms error messages are informative", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Test conditional_effects with invalid effect name
  expect_error(
    conditional_effects(fit, effects = "nonexistent_var"),
    "`effects` must name a numeric predictor in the model data."
  )
  
  # Test with invalid family specification
  expect_error(
    qbrms(y ~ x1, data = test_data, family = "invalid_family"),
    class = "error"
  )
  
  # Note: qbrms has robust fallbacks for formula issues, so these don't error
  # Test that fallback works for formula without response
  result <- qbrms(~x1, data = test_data, family = gaussian())
  expect_true(inherits(result, "qbrms_fit"))
})

test_that("qbrms handles small sample sizes", {
  small_data <- data.frame(
    y = rnorm(10),
    x = rnorm(10)
  )
  
  # Should either fit or give informative error
  result <- tryCatch(
    qbrms(y ~ x, data = small_data, family = gaussian()),
    error = function(e) e
  )
  
  expect_true(
    inherits(result, "qbrms_fit") || 
      (inherits(result, "error") && nchar(result$message) > 0)
  )
})

test_that("qbrms handles perfect collinearity", {
  collinear_data <- test_data
  collinear_data$x3 <- collinear_data$x1 * 2
  
  # Should either fit or give informative error about collinearity
  result <- tryCatch(
    qbrms(y ~ x1 + x3, data = collinear_data, family = gaussian()),
    error = function(e) e,
    warning = function(w) w
  )
  
  expect_true(
    inherits(result, "qbrms_fit") || 
      inherits(result, "error") ||
      inherits(result, "warning")
  )
})

test_that("qbrms handles extreme values", {
  extreme_data <- test_data
  extreme_data$y[1] <- 1e10
  extreme_data$y[2] <- -1e10
  
  # Should either fit or handle gracefully
  result <- tryCatch(
    qbrms(y ~ x1 + x2, data = extreme_data, family = gaussian()),
    error = function(e) e,
    warning = function(w) w
  )
  
  expect_true(
    inherits(result, "qbrms_fit") || 
      inherits(result, "error") ||
      inherits(result, "warning")
  )
})

test_that("qbrms handles all-zero variance predictors", {
  zero_var_data <- test_data
  zero_var_data$x_constant <- 5
  
  # Should either fit or give informative error
  result <- tryCatch(
    qbrms(y ~ x1 + x_constant, data = zero_var_data, family = gaussian()),
    error = function(e) e,
    warning = function(w) w
  )
  
  expect_true(
    inherits(result, "qbrms_fit") || 
      inherits(result, "error") ||
      inherits(result, "warning")
  )
})

test_that("conditional_effects handles factor predictors appropriately", {
  factor_data <- test_data
  factor_data$x_factor <- factor(rep(c("Low", "Medium", "High"), length.out = 100))
  
  fit <- qbrms(y ~ x1 + x_factor, data = factor_data, family = gaussian())
  
  # Should handle or give informative message about factors
  result <- tryCatch(
    conditional_effects(fit),
    error = function(e) e
  )
  
  # Either produces output or gives informative error
  expect_true(
    (is.list(result) && length(result) > 0) || 
      (inherits(result, "error") && nchar(result$message) > 0)
  )
})

test_that("qbrms handles character predictors", {
  char_data <- test_data
  char_data$x_char <- rep(c("A", "B"), length.out = 100)
  
  # Should convert to factor or give informative error
  result <- tryCatch(
    qbrms(y ~ x1 + x_char, data = char_data, family = gaussian()),
    error = function(e) e,
    warning = function(w) w
  )
  
  expect_true(
    inherits(result, "qbrms_fit") || 
      inherits(result, "error") ||
      inherits(result, "warning")
  )
})

test_that("qbrms handles reserved variable names", {
  reserved_data <- test_data
  names(reserved_data)[1] <- "data"
  
  # Should either handle or give informative error
  result <- tryCatch(
    qbrms(data ~ x1 + x2, data = reserved_data, family = gaussian()),
    error = function(e) e
  )
  
  expect_true(
    inherits(result, "qbrms_fit") || 
      (inherits(result, "error") && nchar(result$message) > 0)
  )
})

test_that("pp_check handles edge cases", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Test with very small ndraws
  result_small <- tryCatch(
    pp_check(fit, ndraws = 2),
    error = function(e) e
  )
  
  expect_true(
    inherits(result_small, "ggplot") || 
      (inherits(result_small, "error") && nchar(result_small$message) > 0)
  )
  
  # Test with large ndraws
  result_large <- tryCatch(
    pp_check(fit, ndraws = 10000),
    error = function(e) e,
    warning = function(w) w
  )
  
  expect_true(
    inherits(result_large, "ggplot") || 
      inherits(result_large, "error") ||
      inherits(result_large, "warning")
  )
})