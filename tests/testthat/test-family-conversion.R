# =============================================================================
# FILE: tests/testthat/test-family-conversion.R
# =============================================================================

test_that("convert_family_to_inla handles edge cases gracefully", {
  # NULL input - should default
  result <- convert_family_to_inla(NULL)
  expect_equal(result, "gaussian")
  
  # Empty character - should error with the actual error message format
  expect_error(
    convert_family_to_inla(""),
    "Family '' is not supported in qbrms"
  )
  
  # Invalid function
  bad_function <- function() stop("Error!")
  expect_error(
    convert_family_to_inla(bad_function),
    "Family.*is not supported in qbrms"
  )
})

test_that("validate_family_quantile works correctly", {
  # Should pass for valid combinations
  expect_true(validate_family_quantile("asymmetric_laplace", 0.5))
  expect_true(validate_family_quantile("asymmetric_laplace", 0.9))
  
  # Should pass for gaussian with NULL quantile  
  expect_silent(validate_family_quantile("gaussian", NULL))
  
  # Should error for invalid quantile values
  expect_error(
    validate_family_quantile("asymmetric_laplace", 0),
    "quantile must be a single numeric value between 0 and 1"
  )
  
  expect_error(
    validate_family_quantile("asymmetric_laplace", 1.5),
    "quantile must be a single numeric value between 0 and 1"
  )
  
  expect_error(
    validate_family_quantile("asymmetric_laplace", c(0.25, 0.75)),
    "quantile must be a single numeric value between 0 and 1"
  )
  
  # Should error for families that don't support quantiles - updated error message
  expect_error(
    validate_family_quantile("gaussian", 0.5),
    "Family 'gaussian' does not support quantile regression"
  )
  
  expect_error(
    validate_family_quantile("poisson", 0.9),
    "Family 'poisson' does not support quantile regression"
  )
})
