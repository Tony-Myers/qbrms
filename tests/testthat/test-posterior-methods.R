# test-posterior-methods.R

library(testthat)

test_that("posterior predictive checks produce plots and handle all types", {
  skip_if_not_installed("ggplot2")
  
  set.seed(123)
  n <- 50
  test_data <- data.frame(
    y = rnorm(n, 5, 2),
    x = rnorm(n)
  )
  
  model <- qbrms::qbrms(y ~ x, data = test_data, family = gaussian(), verbose = FALSE)
  
  expect_silent({
    p1 <- pp_check(model, type = "dens_overlay", ndraws = 10)
    print(p1)
    p2 <- pp_check(model, type = "scatter_avg", ndraws = 10)
    print(p2)
    p3 <- pp_check(model, type = "hist", ndraws = 10)
    print(p3)
  })
  
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("pp_check errors on unsupported plot types", {
  skip_if_not_installed("ggplot2")
  
  set.seed(123)
  test_data <- data.frame(
    y = rnorm(50),
    x = rnorm(50)
  )
  
  model <- qbrms::qbrms(y ~ x, data = test_data, family = gaussian(), verbose = FALSE)
  
  expect_error(
    pp_check(model, type = "invalid_type"),
    "Unsupported pp_check type"
  )
})

test_that("pp_check works with quantile regression models", {
  skip_if_not_installed("ggplot2")
  set.seed(456)
  
  n <- 30
  test_data <- data.frame(
    y = rnorm(n, 3, 1),
    x = rnorm(n)
  )
  
  model <- qbrms::qbrms(y ~ x, data = test_data, family = "asymmetric_laplace",
                        quantile = 0.5, verbose = FALSE)
  
  expect_silent({
    p_hist <- pp_check(model, type = "hist", ndraws = 5)
    print(p_hist)
  })
  
  expect_s3_class(p_hist, "ggplot")
})

test_that("pp_check handles data with extreme values gracefully", {
  skip_if_not_installed("ggplot2")
  
  set.seed(789)
  n <- 25
  test_data <- data.frame(
    y = c(rnorm(20, 0, 1), rep(100, 5)),  # 5 extreme outliers at 100
    x = rnorm(n)
  )
  
  model <- qbrms::qbrms(y ~ x, data = test_data, family = gaussian(), verbose = FALSE)
  
  expect_silent({
    p_dens <- pp_check(model, type = "dens_overlay", ndraws = 5)
    print(p_dens)
  })
  
  expect_s3_class(p_dens, "ggplot")
})
