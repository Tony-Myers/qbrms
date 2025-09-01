# test-model-fitting.R

library(testthat)

test_that("basic fixed and mixed effects model fitting without factor arithmetic errors", {
  set.seed(123)
  n <- 50
  
  # Create data with proper factor grouping; no arithmetic on factors
  data <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    group = factor(rep(1:10, each = 5))
  )
  
  # Introduce missing values to check missing data handling
  data$y[1:10] <- NA
  
  # Fit fixed effects model
  expect_no_error({
    fixed_mod <- qbrms(y ~ x1, data = data, family = gaussian(), verbose = FALSE)
  })
  expect_s3_class(fixed_mod, "qbrms_fit")
  
  # Remove missing rows for mixed effects fit
  complete_data <- data[complete.cases(data), ]
  
  expect_true(is.factor(complete_data$group))
  
  # Fit mixed effects model ensuring no arithmetic on factor
  expect_no_error({
    mixed_mod <- qbrms(y ~ x1 + (1 | group), data = complete_data, family = gaussian(), verbose = FALSE)
  })
  expect_s3_class(mixed_mod, "qbrms_fit")
  expect_equal(mixed_mod$group_var, "group")
})

test_that("binomial model fits with successes/failures matrix", {
  set.seed(123)
  n <- 30
  trials <- sample(5:10, n, replace = TRUE)
  successes <- rbinom(n, trials, 0.3)
  failures <- trials - successes
  
  successes <- as.integer(successes)
  failures <- as.integer(failures)
  
  data_binom <- data.frame(
    x = rnorm(n),
    successes = successes,
    failures = failures
  )
  
  expect_true(all(successes >= 0 & successes <= trials))
  expect_true(all(failures >= 0))
  
  expect_no_error({
    model_binom <- qbrms(
      cbind(successes, failures) ~ x,
      data = data_binom,
      family = binomial(),
      verbose = TRUE  # enable verbose to see INLA diagnostics
    )
  })
  expect_s3_class(model_binom, "qbrms_fit")
})

test_that("asymmetric laplace quantile regression fits", {
  set.seed(123)
  n <- 35
  
  data <- data.frame(
    y = rnorm(n, 7, 2),
    x = rnorm(n)
  )
  
  expect_no_error({
    qr_mod <- qbrms(y ~ x, data = data, family = "asymmetric_laplace", quantile = 0.25, verbose = FALSE)
  })
  expect_s3_class(qr_mod, "qbrms_fit")
  expect_equal(qr_mod$model_type, "quantile_regression")
})

test_that("ordinal regression requires at least 3 levels", {
  set.seed(123)
  n <- 40
  
  data <- data.frame(
    x = rnorm(n),
    y = factor(sample(1:2, n, replace = TRUE), ordered = TRUE)
  )
  
  expect_error(
    qbrms(y ~ x, data = data, family = cumulative()),
    "Need at least 3 ordinal levels"
  )
})

test_that("ordinal regression fits with 3 or more ordered levels", {
  set.seed(123)
  n <- 60
  
  data <- data.frame(
    treatment = factor(rep(c("A", "B"), each = n/2)),
    age = rnorm(n),
    baseline = rnorm(n)
  )
  
  linpred <- 0.5 * (data$treatment == "B") + 0.02 * data$age + 0.3 * data$baseline + rnorm(n)
  
  data$response <- cut(linpred, 
                       breaks = c(-Inf, -0.5, 0, 0.5, 1, Inf),
                       labels = c("Very Poor", "Poor", "Average", "Good", "Excellent"),
                       ordered_result = TRUE)
  
  expect_true(is.ordered(data$response))
  
  expect_no_error({
    mod_ord <- qbrms(response ~ treatment + age + baseline, data = data, family = cumulative(), verbose = FALSE)
  })
  expect_s3_class(mod_ord, c("ordinal_augmented_qbrms_fit", "qbrms_fit"))
  expect_equal(length(mod_ord$ordinal_levels), 5)
})

