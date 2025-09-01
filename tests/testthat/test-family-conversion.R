test_that("family conversion works correctly", {
  # Test basic families
  expect_equal(convert_family_to_inla(gaussian()), "gaussian")
  expect_equal(convert_family_to_inla("binomial"), "binomial")
  expect_equal(convert_family_to_inla("poisson"), "poisson")
  
  # Test asymmetric laplace
  al_result <- convert_family_to_inla("asymmetric_laplace", quantile = 0.3)
  expect_equal(al_result$family, "asymmetric_laplace")
  expect_equal(al_result$quantile, 0.3)
  
  # Test function input
  expect_equal(convert_family_to_inla(binomial()), "binomial")
  
  # Test unsupported family
  expect_error(convert_family_to_inla("unsupported_family"))
})

test_that("family helper functions exist", {
  expect_s3_class(negbinomial(), "family")
  expect_s3_class(weibull(), "family")
  expect_s3_class(zero_inflated_poisson(), "family")
  expect_s3_class(skew_normal(), "family")
  expect_s3_class(student_t(), "family")
})
