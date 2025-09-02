# Complete replacement for your test-posterior-methods.R file

test_that("posterior predictive checks produce plots and handle all types", {
  data <- data.frame(y = rnorm(50), x = rnorm(50))
  
  # Fit model with verbose=FALSE and suppress all messages
  fit <- suppressMessages(qbrms(y ~ x, data = data, verbose = FALSE))
  
  # Test all plot types with complete message suppression
  expect_no_error(suppressMessages(suppressWarnings({
    p1 <- pp_check(fit, type = "dens_overlay")
    p2 <- pp_check(fit, type = "hist") 
    p3 <- pp_check(fit, type = "scatter")
    p4 <- pp_check(fit, type = "scatter_avg")
  })))
  
  # Test that plots are ggplot objects (if ggplot2 is available)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    expect_s3_class(p3, "ggplot")
    expect_s3_class(p4, "ggplot")
  }
})

test_that("pp_check handles different model types", {
  data <- data.frame(y = rnorm(30), x = rnorm(30))
  
  # Test with different families, all with verbose=FALSE
  suppressMessages({
    fit_gaussian <- qbrms(y ~ x, data = data, family = gaussian(), verbose = FALSE)
    fit_quantile <- qbrms(y ~ x, data = data, family = asymmetric_laplace(), quantile = 0.5, verbose = FALSE)
  })
  
  expect_no_error(suppressMessages({
    pp_check(fit_gaussian)
    pp_check(fit_quantile)
  }))
})

test_that("pp_check handles prior predictive correctly", {
  data <- data.frame(y = rnorm(20), x = rnorm(20))
  
  # Prior predictive check
  prior_fit <- suppressMessages(qbrms(y ~ x, data = data, sample_prior = "only", verbose = FALSE))
  
  expect_no_error(suppressMessages({
    p_prior <- pp_check(prior_fit, type = "dens_overlay")
  }))
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(p_prior, c("qbrms_prior_ggplot", "ggplot"))
  }
})