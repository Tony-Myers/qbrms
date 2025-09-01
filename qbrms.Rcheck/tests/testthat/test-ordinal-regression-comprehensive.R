# tests/testthat/test-ordinal-regression-comprehensive.R

library(testthat)

test_that("ordinal regression errors on less than 3 levels", {
  set.seed(123)
  
  data_2cat <- data.frame(
    x = rnorm(30),
    y = factor(sample(1:2, 30, replace = TRUE), ordered = TRUE)
  )
  
  # Expect error for <3 levels
  expect_error(qbrms(y ~ x, data = data_2cat, family = cumulative()),
               "Need at least 3 ordinal levels")
})

test_that("ordinal data augmentation approach works properly", {
  set.seed(456)
  n <- 60
  
  data <- data.frame(
    treatment = factor(rep(c("A", "B"), each = n/2)),
    age = rnorm(n, 50, 10),
    baseline_score = rnorm(n)
  )
  
  linear_pred <- 0.5 * (data$treatment == "B") + 0.02 * data$age + 0.3 * data$baseline_score + rnorm(n, 0, 0.5)
  
  data$response <- cut(linear_pred,
                       breaks = c(-Inf, -0.5, 0, 0.5, 1, Inf),
                       labels = c("Very Poor", "Poor", "Average", "Good", "Excellent"),
                       ordered_result = TRUE)
  
  expect_no_error({
    model_ord <- qbrms(response ~ treatment + age + baseline_score, data = data,
                       family = cumulative(), verbose = FALSE)
  })
  
  expect_s3_class(model_ord, c("ordinal_augmented_qbrms_fit", "qbrms_fit"))
  expect_equal(length(model_ord$ordinal_levels), 5)
})

test_that("ordinal binary decomposition approach works", {
  set.seed(789)
  n <- 40
  
  data <- data.frame(
    x = rnorm(n),
    y = factor(sample(1:4, n, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2)),
               ordered = TRUE)
  )
  
  expect_no_error({
    model_bin <- qbrms_ordinal_binary(y ~ x, data = data, verbose = FALSE)
  })
  
  expect_s3_class(model_bin, c("ordinal_binary_qbrms_fit", "qbrms_fit"))
  expect_equal(length(model_bin$binary_models), 3)
})
