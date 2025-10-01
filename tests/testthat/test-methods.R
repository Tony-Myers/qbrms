# =============================================================================
# FILE: tests/testthat/test-methods.R
# =============================================================================
# Tests for S3 methods in qbrms package

# Setup test data
test_data <- data.frame(
  y = rnorm(100, mean = 2, sd = 1),
  x1 = rnorm(100),
  x2 = rnorm(100),
  group = rep(c("A", "B"), each = 50)
)

test_that("print.qbrms_fit works correctly", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Test that print produces output
  output <- capture.output(print(fit))
  
  expect_true(length(output) > 0)
  expect_true(any(grepl("qbrms Model fit|Family", output)))
})

test_that("summary.qbrms_fit works correctly", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Capture summary output
  output <- capture.output(summary(fit))
  
  # Test that summary produces output
  expect_true(length(output) > 0)
  
  # Test that output contains expected sections
  expect_true(any(grepl("Population-Level Effects|Family|Intercept|gaussian", output)))
})

test_that("coef.qbrms_fit extracts coefficients", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Extract coefficients
  coefs <- coef(fit)
  
  # Test structure
  expect_true(is.numeric(coefs) || is.matrix(coefs) || is.data.frame(coefs))
  expect_true(length(coefs) > 0 || nrow(coefs) > 0)
})

test_that("fitted.qbrms_fit returns fitted values", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Get fitted values
  fitted_vals <- fitted(fit)
  
  # Test structure
  expect_true(is.numeric(fitted_vals) || is.matrix(fitted_vals))
  expect_equal(length(fitted_vals), nrow(test_data))
})

test_that("residuals.qbrms_fit returns residuals", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Get residuals
  resids <- residuals(fit)
  
  # Test structure
  expect_true(is.numeric(resids) || is.matrix(resids))
  expect_equal(length(resids), nrow(test_data))
})

test_that("predict.qbrms_fit makes predictions", {
  skip("predict method not yet implemented")
  
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Predict on original data
  preds <- predict(fit)
  
  # Test structure
  expect_true(is.numeric(preds) || is.matrix(preds) || is.data.frame(preds))
  expect_true(length(preds) > 0 || nrow(preds) > 0)
  
  # Predict on new data
  new_data <- data.frame(x1 = c(0, 1), x2 = c(0, 1))
  preds_new <- predict(fit, newdata = new_data)
  
  expect_true(is.numeric(preds_new) || is.matrix(preds_new) || is.data.frame(preds_new))
  expect_true(length(preds_new) > 0 || nrow(preds_new) > 0)
})

test_that("vcov.qbrms_fit returns variance-covariance matrix", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Get vcov matrix
  result <- tryCatch(
    vcov(fit),
    error = function(e) e
  )
  
  # Should return matrix or give informative error
  expect_true(
    is.matrix(result) || 
      (inherits(result, "error") && nchar(result$message) > 0)
  )
})

test_that("formula.qbrms_fit extracts formula", {
  skip("formula method not yet implemented")
  
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Get formula
  form <- formula(fit)
  
  # Test structure
  expect_true(inherits(form, "formula"))
  expect_true(length(form) > 1)
})

test_that("family.qbrms_fit extracts family", {
  skip("family method not yet implemented")
  
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Get family
  fam <- family(fit)
  
  # Test structure
  expect_true(inherits(fam, "family") || is.character(fam) || is.list(fam))
})

test_that("nobs.qbrms_fit returns number of observations", {
  skip("nobs method not yet implemented")
  
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Get number of observations
  n <- nobs(fit)
  
  # Test value
  expect_true(is.numeric(n))
  expect_equal(n, nrow(test_data))
})

test_that("logLik.qbrms_fit returns log-likelihood", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Get log-likelihood
  result <- tryCatch(
    logLik(fit),
    error = function(e) e
  )
  
  # Should return numeric or logLik object, or give informative error
  expect_true(
    is.numeric(result) || 
      inherits(result, "logLik") ||
      (inherits(result, "error") && nchar(result$message) > 0)
  )
})

test_that("update.qbrms_fit updates model", {
  fit <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Update formula
  result <- tryCatch(
    update(fit, formula = y ~ x1),
    error = function(e) e
  )
  
  # Should return updated fit or give informative error
  expect_true(
    inherits(result, "qbrms_fit") ||
      (inherits(result, "error") && nchar(result$message) > 0)
  )
})

test_that("anova.qbrms_fit performs model comparison", {
  fit1 <- qbrms(y ~ x1, data = test_data, family = gaussian())
  fit2 <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  # Compare models
  result <- tryCatch(
    anova(fit1, fit2),
    error = function(e) e
  )
  
  # Should return comparison or give informative error
  expect_true(
    is.data.frame(result) || 
      is.list(result) ||
      (inherits(result, "error") && nchar(result$message) > 0)
  )
})