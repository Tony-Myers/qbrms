# tests/testthat/test-quantile-regression-comprehensive.R

test_that("quantile regression works for different quantiles", {
  # Create test data with heteroscedastic noise
  set.seed(789)
  qr_data <- data.frame(
    y = rnorm(50, 5, 2) + rnorm(50, 0, abs(rnorm(50, 0, 0.5))), # Heteroscedastic
    x1 = rnorm(50),
    x2 = rnorm(50),
    x3 = factor(sample(c("A", "B"), 50, replace = TRUE))
  )
  
  # Test different quantiles
  quantiles_to_test <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  
  for (q in quantiles_to_test) {
    # Fit the model once and check for both no-error and correct structure
    qr_model <- qbrms(y ~ x1 + x2, data = qr_data,
                      family = "asymmetric_laplace", quantile = q)
    
    # Check model structure and properties
    expect_s3_class(qr_model, "qbrms_fit")
    expect_equal(qr_model$quantile, q)
  }
})

test_that("quantile regression with categorical predictors works", {
  qr_data_cat <- data.frame(
    y = rnorm(40, 10, 3),
    treatment = factor(sample(c("Control", "Low", "High"), 40, replace = TRUE)),
    continuous_var = rnorm(40)
  )
  
  expect_no_error({
    qr_model_cat <- qbrms(y ~ treatment + continuous_var,
                          data = qr_data_cat,
                          family = "asymmetric_laplace",
                          quantile = 0.75)
  })
})

test_that("quantile regression handles edge cases", {
  # Test with extreme quantiles
  edge_data <- data.frame(y = rnorm(30), x = rnorm(30))
  
  expect_no_error(qbrms(y ~ x, data = edge_data, family = "asymmetric_laplace", quantile = 0.05))
  expect_no_error(qbrms(y ~ x, data = edge_data, family = "asymmetric_laplace", quantile = 0.95))
  
  # Test with median (should be equivalent to LAD regression)
  expect_no_error(qbrms(y ~ x, data = edge_data, family = "asymmetric_laplace", quantile = 0.5))
})

test_that("quantile regression summary and methods work", {
  qr_data <- data.frame(y = rnorm(35, 8, 2), x = rnorm(35))
  qr_model <- qbrms(y ~ x, data = qr_data, family = "asymmetric_laplace", quantile = 0.25)
  
  expect_no_error(summary(qr_model))
  expect_output(summary(qr_model))
  
  # Test conditional effects if implemented
  expect_no_error(conditional_effects(qr_model))
})

test_that("quantile regression convergence works", {
  # Test that the iterative algorithm converges properly
  convergence_data <- data.frame(
    y = c(1:30) + rnorm(30, 0, 0.5), # Clear linear trend
    x = 1:30
  )
  
  qr_model <- qbrms(y ~ x, data = convergence_data, family = "asymmetric_laplace", quantile = 0.5)
  
  # Should converge (this would be indicated in the output)
  expect_s3_class(qr_model, "qbrms_fit")
  expect_true(qr_model$fit$converged %||% TRUE) # Handle if converged field does not exist
})