# test-prior-predictive.R (snippet focusing on plotting fixes)

library(testthat)

test_that("prior predictive plots render correctly", {
  skip_if_not_installed("ggplot2")
  
  set.seed(123)
  n <- 30
  
  data <- data.frame(
    y = rnorm(n),
    x = rnorm(n)
  )
  
  prior_fit <- qbrms(y ~ x, data = data, family = gaussian(), sample_prior = "only", verbose = FALSE)
  
  p1 <- pp_check(prior_fit, type = "dens_overlay")
  expect_s3_class(p1, "ggplot")
  expect_output(print(p1))
  
  p2 <- pp_check(prior_fit, type = "hist")
  expect_s3_class(p2, "ggplot")
  expect_output(print(p2))
  
  p3 <- pp_check(prior_fit, type = "scatter_avg")
  expect_s3_class(p3, "ggplot")
  expect_output(print(p3))
})

