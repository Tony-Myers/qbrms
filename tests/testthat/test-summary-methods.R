# tests/testthat/test-summary-methods.R

test_that("summary and print methods produce correct output", {
  set.seed(123)
  test_data <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    group = factor(rep(1:5, each = 10))
  )
  
  # --- Test a fixed effects model ---
  model1 <- qbrms(y ~ x, data = test_data, family = gaussian(), verbose = FALSE)
  
  # Test print method
  print_output1 <- capture.output(print(model1))
  expect_gt(length(print_output1), 3)
  expect_true(any(grepl("qbrms Model fit", print_output1)))
  
  # Test summary method
  summary_output1 <- capture.output(summary(model1))
  expect_true(any(grepl("Population-Level Effects:", summary_output1)))
  expect_true(any(grepl("Family Specific Parameters:", summary_output1)))
  
  # --- Test a mixed effects model ---
  model2 <- qbrms(y ~ x + (1|group), data = test_data, family = gaussian(), verbose = FALSE)
  
  # Test print method
  print_output2 <- capture.output(print(model2))
  expect_true(any(grepl("Random Effects:  group", print_output2)))
  
  # Test summary method
  summary_output2 <- capture.output(summary(model2))
  expect_true(any(grepl("Group-Level Effects:", summary_output2)))
})

test_that("timing information is captured", {
  set.seed(456)
  test_data <- data.frame(y = rnorm(30), x = rnorm(30))
  model <- qbrms(y ~ x, data = test_data, family = gaussian(), verbose = FALSE)
  
  expect_true(!is.null(model$timing))
  expect_true(is.numeric(model$timing$total_seconds))
  expect_true(is.character(model$timing$formatted_duration))
})