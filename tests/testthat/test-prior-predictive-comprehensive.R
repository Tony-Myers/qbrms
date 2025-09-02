# Complete replacement for your test-prior-predictive-comprehensive.R file

test_that("prior predictive plots render correctly", {
  data <- data.frame(y = rnorm(50), x = rnorm(50))
  
  # Generate prior predictive samples with verbose=FALSE
  prior_fit <- suppressMessages(qbrms(y ~ x, data = data, sample_prior = "only", verbose = FALSE))
  
  # Create plots with message suppression
  suppressMessages({
    p1 <- pp_check(prior_fit, type = "dens_overlay")
    p2 <- pp_check(prior_fit, type = "hist") 
    p3 <- pp_check(prior_fit, type = "scatter")
  })
  
  # Test that plots produce output when printed
  # Use expect_output instead of testing print directly
  expect_output(print(p1), "Prior Predictive|\\[|ggplot", perl = TRUE)
  expect_output(print(p2), "Prior Predictive|\\[|ggplot", perl = TRUE)  
  expect_output(print(p3), "Prior Predictive|\\[|ggplot", perl = TRUE)
  
  # Test plot classes
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(p1, c("qbrms_prior_ggplot", "ggplot"))
    expect_s3_class(p2, "ggplot")
    expect_s3_class(p3, "ggplot")
  }
})

test_that("prior predictive sampling works with different families", {
  data <- data.frame(y = rnorm(30), x = rnorm(30))
  
  suppressMessages({
    # Test different families
    prior_gaussian <- qbrms(y ~ x, data = data, family = gaussian(), 
                            sample_prior = "only", verbose = FALSE)
    prior_laplace <- qbrms(y ~ x, data = data, family = asymmetric_laplace(), 
                           sample_prior = "only", verbose = FALSE)
  })
  
  # Check that prior samples are generated
  expect_true(!is.null(prior_gaussian$prior_samples))
  expect_true(!is.null(prior_laplace$prior_samples))
  
  # Check that samples are matrices
  expect_true(is.matrix(prior_gaussian$prior_samples))
  expect_true(is.matrix(prior_laplace$prior_samples))
})

test_that("prior predictive checks handle various plot types", {
  data <- data.frame(y = rnorm(25), x = rnorm(25))
  
  prior_fit <- suppressMessages(qbrms(y ~ x, data = data, sample_prior = "only", verbose = FALSE))
  
  # Test all plot types work without error
  expect_no_error(suppressMessages({
    p_overlay <- pp_check(prior_fit, type = "dens_overlay")
    p_hist <- pp_check(prior_fit, type = "hist")
    p_scatter <- pp_check(prior_fit, type = "scatter")
    p_scatter_avg <- pp_check(prior_fit, type = "scatter_avg")
  }))
  
  # Ensure they're all plottable objects
  plots <- list(p_overlay, p_hist, p_scatter, p_scatter_avg)
  expect_true(all(sapply(plots, function(p) {
    inherits(p, "ggplot") || inherits(p, c("qbrms_plot", "qbrms_prior_plot"))
  })))
})