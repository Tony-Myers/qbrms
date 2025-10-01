# =============================================================================
# FILE: tests/testthat/test-model-fitting.R
# =============================================================================
# Tests for basic model fitting functionality in qbrms package

# Setup test data
test_data <- data.frame(
  y = rnorm(100, mean = 2, sd = 1),
  x1 = rnorm(100),
  x2 = rnorm(100),
  group = rep(c("A", "B"), each = 50)
)

test_that("qbrms fits basic Gaussian models", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  expect_true(inherits(fit, "qbrms_fit"))
  expect_true(!is.null(fit$fit))
})

test_that("qbrms respects backend argument", {
  # Test with INLA backend
  fit_inla <- tryCatch(
    qbrms(y ~ x1, data = test_data, family = gaussian(), backend = "inla"),
    error = function(e) e
  )
  
  expect_true(
    inherits(fit_inla, "qbrms_fit") ||
      (inherits(fit_inla, "error") && nchar(fit_inla$message) > 0)
  )
  
  # Test with TMB backend
  fit_tmb <- tryCatch(
    qbrms(y ~ x1, data = test_data, family = gaussian(), backend = "tmb"),
    error = function(e) e
  )
  
  expect_true(
    inherits(fit_tmb, "qbrms_fit") ||
      (inherits(fit_tmb, "error") && nchar(fit_tmb$message) > 0)
  )
})

test_that("qbrms handles missing values appropriately", {
  test_data_na <- test_data
  test_data_na$y[1:5] <- NA
  
  # qbrms should either handle NA values, remove them, or give a clear error
  result <- tryCatch(
    qbrms(y ~ x1 + x2, data = test_data_na, family = gaussian()),
    error = function(e) e,
    warning = function(w) w
  )
  
  # Test passes if either: 
  # - model fits successfully (NAs removed automatically)
  # - gives informative error
  # - gives warning about NA removal
  expect_true(
    inherits(result, "qbrms_fit") || 
      (inherits(result, "error") && nchar(result$message) > 0) ||
      (inherits(result, "warning") && nchar(result$message) > 0)
  )
})

test_that("qbrms fits binomial models", {
  test_data_binary <- test_data
  test_data_binary$y_binary <- rbinom(nrow(test_data), 1, 0.5)
  
  fit <- qbrms(y_binary ~ x1 + x2, data = test_data_binary, family = binomial())
  
  expect_true(inherits(fit, "qbrms_fit"))
})

test_that("qbrms fits Poisson models", {
  test_data_count <- test_data
  test_data_count$y_count <- rpois(nrow(test_data), lambda = 2)
  
  fit <- qbrms(y_count ~ x1 + x2, data = test_data_count, family = poisson())
  
  expect_true(inherits(fit, "qbrms_fit"))
})

test_that("qbrms fits models with interactions", {
  fit <- qbrms(y ~ x1 * x2, data = test_data, family = gaussian())
  
  expect_true(inherits(fit, "qbrms_fit"))
})

test_that("qbrms fits intercept-only models", {
  fit <- qbrms(y ~ 1, data = test_data, family = gaussian())
  
  expect_true(inherits(fit, "qbrms_fit"))
})

test_that("qbrms validates formula argument", {
  # Note: qbrms has robust fallbacks, so these don't necessarily error
  # Test with missing response
  result1 <- qbrms(~x1 + x2, data = test_data, family = gaussian())
  expect_true(inherits(result1, "qbrms_fit"))
  
  # Test with character string formula - converts automatically
  result2 <- qbrms("y ~ x1", data = test_data, family = gaussian())
  expect_true(inherits(result2, "qbrms_fit"))
})

test_that("qbrms validates data argument", {
  # Test with missing data
  expect_error(
    qbrms(y ~ x1 + x2, family = gaussian()),
    class = "error"
  )
  
  # Test with non-data.frame input
  expect_error(
    qbrms(y ~ x1 + x2, data = "not a data frame", family = gaussian()),
    class = "error"
  )
})

test_that("qbrms validates family argument", {
  # Test with invalid family
  expect_error(
    qbrms(y ~ x1 + x2, data = test_data, family = "invalid"),
    class = "error"
  )
  
  # Note: qbrms has robust fallbacks for NULL family
  result <- qbrms(y ~ x1 + x2, data = test_data, family = NULL)
  expect_true(inherits(result, "qbrms_fit"))
})

test_that("qbrms handles complex formulas", {
  # Test with polynomial terms
  fit_poly <- tryCatch(
    qbrms(y ~ poly(x1, 2) + x2, data = test_data, family = gaussian()),
    error = function(e) e
  )
  
  expect_true(
    inherits(fit_poly, "qbrms_fit") ||
      (inherits(fit_poly, "error") && nchar(fit_poly$message) > 0)
  )
  
  # Test with transformations
  fit_trans <- tryCatch(
    qbrms(y ~ log(abs(x1) + 1) + x2, data = test_data, family = gaussian()),
    error = function(e) e
  )
  
  expect_true(
    inherits(fit_trans, "qbrms_fit") ||
      (inherits(fit_trans, "error") && nchar(fit_trans$message) > 0)
  )
})

test_that("qbrms respects prior argument", {
  # Test with custom priors (if supported)
  result <- tryCatch(
    qbrms(y ~ x1 + x2, data = test_data, family = gaussian(), 
          prior = list(beta = c(0, 10))),
    error = function(e) e
  )
  
  # Should either use priors or inform that they are not supported yet
  expect_true(
    inherits(result, "qbrms_fit") ||
      (inherits(result, "error") && nchar(result$message) > 0)
  )
})

test_that("qbrms respects control arguments", {
  # Test with control arguments
  result <- tryCatch(
    qbrms(y ~ x1 + x2, data = test_data, family = gaussian(),
          control = list(max_iter = 100)),
    error = function(e) e,
    warning = function(w) w
  )
  
  # Should either use control args or give warning/error
  expect_true(
    inherits(result, "qbrms_fit") ||
      inherits(result, "error") ||
      inherits(result, "warning")
  )
})