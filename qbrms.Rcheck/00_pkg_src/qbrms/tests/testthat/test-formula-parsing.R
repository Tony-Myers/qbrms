test_that("formula parsing works correctly", {
  # Test basic formula
  test_data <- data.frame(
    y = rnorm(50),
    x1 = rnorm(50),
    x2 = rnorm(50),
    group = factor(rep(1:5, each = 10))
  )
  
  # Test fixed effects formula
  components1 <- parse_formula_components(y ~ x1 + x2, test_data)
  expect_false(components1$has_random_effects)
  expect_false(components1$is_binomial_trials)
  
  # Test mixed effects formula
  components2 <- parse_formula_components(y ~ x1 + (1|group), test_data)
  expect_true(components2$has_random_effects)
  expect_false(components2$is_binomial_trials)
  
  # Test binomial trials
  binom_data <- data.frame(
    success = rbinom(30, 10, 0.3),
    trials = rep(10, 30),
    x = rnorm(30)
  )
  
  components3 <- parse_formula_components(success | trials(trials) ~ x, binom_data)
  expect_true(components3$is_binomial_trials)
})

test_that("predictor variable extraction works", {
  test_data <- data.frame(
    y = rnorm(50),
    x1 = rnorm(50),
    x2 = factor(rep(letters[1:5], each = 10)),
    x3 = rnorm(50)
  )
  
  predictors <- get_predictor_variables(y ~ x1 + x2 + x3, test_data)
  
  expect_true("x1" %in% predictors$numeric_vars)
  expect_true("x3" %in% predictors$numeric_vars)
  expect_true("x2" %in% predictors$categorical_vars)
  expect_equal(length(predictors$all_vars), 3)
})
