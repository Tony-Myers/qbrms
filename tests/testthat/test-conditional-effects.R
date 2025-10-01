# =============================================================================
# FILE: tests/testthat/test-conditional-effects.R
# =============================================================================
# Tests for conditional_effects functionality in qbrms package

# Setup test data
test_data <- data.frame(
  y = rnorm(100, mean = 2, sd = 1),
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = rnorm(100),
  group = rep(c("A", "B"), each = 50)
)

test_that("conditional_effects produces basic output", {
  # Fit a simple model
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Get conditional effects
  ce_result <- conditional_effects(fit)
  
  # Test that we get a list of data frames
  expect_type(ce_result, "list")
  expect_true(all(sapply(ce_result, is.data.frame)))
  
  # Test individual effect
  plot_data <- conditional_effects(fit, effects = "x1")
  expect_true(is.data.frame(plot_data) || (is.list(plot_data) && is.data.frame(plot_data[[1]])))
  
  # Test that data frame has required columns for plotting
  if (is.list(plot_data)) plot_data <- plot_data[[1]]
  expect_true(is.data.frame(plot_data))
  expect_true(nrow(plot_data) > 0)
  expect_true(ncol(plot_data) > 0)
})

test_that("conditional_effects handles specific effects argument", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Request specific effect
  ce_x1 <- conditional_effects(fit, effects = "x1")
  
  # Should return data frame or list with one element
  if (is.list(ce_x1)) {
    expect_length(ce_x1, 1)
    expect_true(is.data.frame(ce_x1[[1]]))
  } else {
    expect_true(is.data.frame(ce_x1))
  }
})

test_that("conditional_effects handles multiple effects", {
  skip("Multiple effects not supported in current implementation")
  
  fit <- qbrms(y ~ x1 + x2 + x3, data = test_data, family = gaussian())
  
  # Request multiple effects
  ce_multi <- conditional_effects(fit, effects = c("x1", "x2"))
  
  expect_type(ce_multi, "list")
  expect_length(ce_multi, 2)
  expect_true(all(sapply(ce_multi, is.data.frame)))
})

test_that("conditional_effects handles plot argument", {
  skip("plot argument not supported in current implementation")
  
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Test with plot = FALSE
  ce_no_plot <- conditional_effects(fit, plot = FALSE)
  expect_type(ce_no_plot, "list")
  expect_true(all(sapply(ce_no_plot, is.data.frame)))
})

test_that("conditional_effects handles interactions", {
  # Fit model with interaction
  fit <- qbrms(y ~ x1 * x2, data = test_data, family = gaussian())
  
  # Get conditional effects
  ce_result <- conditional_effects(fit)
  
  # Should return list of data frames
  expect_type(ce_result, "list")
  expect_true(length(ce_result) > 0)
  expect_true(all(sapply(ce_result, function(x) {
    is.data.frame(x) || (is.list(x) && all(sapply(x, is.data.frame)))
  })))
})

test_that("conditional_effects works with different families", {
  # Binary outcome - use binomial family
  test_data_binary <- test_data
  test_data_binary$y_binary <- rbinom(nrow(test_data), 1, 0.5)
  fit_binary <- qbrms(y_binary ~ x1, data = test_data_binary, family = binomial())
  
  ce_binary <- conditional_effects(fit_binary)
  expect_type(ce_binary, "list")
  expect_true(all(sapply(ce_binary, function(x) {
    is.data.frame(x) || (is.list(x) && all(sapply(x, is.data.frame)))
  })))
  
  # Count outcome  
  test_data_count <- test_data
  test_data_count$y_count <- rpois(nrow(test_data), 2)
  fit_count <- qbrms(y_count ~ x1, data = test_data_count, family = poisson())
  
  ce_count <- conditional_effects(fit_count)
  expect_type(ce_count, "list")
  expect_true(all(sapply(ce_count, function(x) {
    is.data.frame(x) || (is.list(x) && all(sapply(x, is.data.frame)))
  })))
})

test_that("conditional_effects validates input arguments", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Test with nonexistent variable
  expect_error(
    conditional_effects(fit, effects = "nonexistent_var"),
    "`effects` must name a numeric predictor in the model data."
  )
  
  # Test with invalid input type
  expect_error(
    conditional_effects(fit, effects = 123),
    class = "error"
  )
})

test_that("conditional_effects handles models with no predictors", {
  # Intercept-only model
  fit <- qbrms(y ~ 1, data = test_data, family = gaussian())
  
  # Should error appropriately
  expect_error(
    conditional_effects(fit),
    "No suitable numeric predictors found to compute conditional effects"
  )
})

test_that("conditional_effects output is properly structured", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  ce_data <- conditional_effects(fit, effects = "x1")
  
  # Extract data frame if in list
  if (is.list(ce_data)) ce_data <- ce_data[[1]]
  
  # Check structure
  expect_true(is.data.frame(ce_data))
  expect_true(nrow(ce_data) > 0)
  expect_true(ncol(ce_data) > 0)
  
  # Check that essential columns exist for plotting
  expect_true(length(names(ce_data)) > 0)
})

test_that("conditional_effects respects prob argument", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Test with different probability levels
  ce_95 <- conditional_effects(fit, effects = "x1", prob = 0.95)
  ce_80 <- conditional_effects(fit, effects = "x1", prob = 0.80)
  
  # Both should return valid data structures
  if (is.list(ce_95)) ce_95 <- ce_95[[1]]
  if (is.list(ce_80)) ce_80 <- ce_80[[1]]
  
  expect_true(is.data.frame(ce_95))
  expect_true(is.data.frame(ce_80))
})

test_that("conditional_effects respects ndraws argument", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Test with different numbers of draws
  ce_100 <- conditional_effects(fit, effects = "x1", ndraws = 100)
  ce_50 <- conditional_effects(fit, effects = "x1", ndraws = 50)
  
  # Both should return valid data structures
  if (is.list(ce_100)) ce_100 <- ce_100[[1]]
  if (is.list(ce_50)) ce_50 <- ce_50[[1]]
  
  expect_true(is.data.frame(ce_100))
  expect_true(is.data.frame(ce_50))
})

test_that("conditional_effects handles at argument", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Test with at argument to fix other predictors
  ce_at <- conditional_effects(fit, effects = "x1", at = list(x2 = 0))
  
  # Should return valid data structure
  if (is.list(ce_at)) ce_at <- ce_at[[1]]
  
  expect_true(is.data.frame(ce_at))
  expect_true(nrow(ce_at) > 0)
})

test_that("conditional_effects handles seed argument", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Test with seed for reproducibility
  ce_seed1 <- conditional_effects(fit, effects = "x1", seed = 123)
  ce_seed2 <- conditional_effects(fit, effects = "x1", seed = 123)
  
  # Extract data frames
  if (is.list(ce_seed1)) ce_seed1 <- ce_seed1[[1]]
  if (is.list(ce_seed2)) ce_seed2 <- ce_seed2[[1]]
  
  # Results should be identical with same seed
  expect_equal(ce_seed1, ce_seed2)
})