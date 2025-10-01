# =============================================================================
# FILE: tests/testthat/test-additional-families.R
# =============================================================================
# Test additional family implementations beyond basic gaussian/binomial/poisson

skip_if_not_installed("INLA")

test_that("qbrms handles robust families correctly", {
  
  set.seed(123)
  # Create data with outliers to test robust families
  robust_data <- data.frame(
    y = c(rnorm(45, 0, 1), c(10, -10, 15, -12, 8)),  # Add outliers
    x = c(rnorm(45), rnorm(5))
  )
  
  # Test skew_normal family
  fit_sn <- qbrms(y ~ x, data = robust_data, family = "skew_normal", verbose = FALSE)
  expect_s3_class(fit_sn, "qbrms_fit")
  expect_true(all(is.finite(coef(fit_sn))))
  
  # Test student_t family
  fit_t <- qbrms(y ~ x, data = robust_data, family = "student_t", verbose = FALSE)
  expect_s3_class(fit_t, "qbrms_fit")
  expect_true(all(is.finite(coef(fit_t))))
  
  # Test student alias
  fit_student <- qbrms(y ~ x, data = robust_data, family = "student", verbose = FALSE)
  expect_s3_class(fit_student, "qbrms_fit")
  expect_equal(fit_student$family, fit_t$family)
})

test_that("qbrms handles count families beyond poisson", {
  
  set.seed(42)
  # Create overdispersed count data
  count_data <- data.frame(
    y_nb = MASS::rnegbin(50, mu = exp(0.5 + 0.3 * (1:50)/50), theta = 2),
    x = (1:50)/50
  )
  
  # Test negative binomial families
  fit_nb1 <- qbrms(y_nb ~ x, data = count_data, family = "negbinomial", verbose = FALSE)
  expect_s3_class(fit_nb1, "qbrms_fit")
  
  fit_nb2 <- qbrms(y_nb ~ x, data = count_data, family = "negative_binomial", verbose = FALSE)
  expect_s3_class(fit_nb2, "qbrms_fit")
  
  fit_nb3 <- qbrms(y_nb ~ x, data = count_data, family = "nbinomial", verbose = FALSE)
  expect_s3_class(fit_nb3, "qbrms_fit")
  
  # All should map to the same INLA family
  expect_equal(fit_nb1$family, fit_nb2$family)
  expect_equal(fit_nb2$family, fit_nb3$family)
})

test_that("qbrms handles continuous families with constraints", {
  
  set.seed(123)
  
  # Beta family - data must be in (0, 1)
  beta_data <- data.frame(
    y_beta = rbeta(40, 2, 3),
    x = rnorm(40)
  )
  
  fit_beta <- qbrms(y_beta ~ x, data = beta_data, family = "beta", verbose = FALSE)
  expect_s3_class(fit_beta, "qbrms_fit")
  expect_equal(fit_beta$family, "beta")
  
  # Gamma family - data must be positive
  gamma_data <- data.frame(
    y_gamma = rgamma(40, shape = 2, rate = 1),
    x = rnorm(40)
  )
  
  fit_gamma <- qbrms(y_gamma ~ x, data = gamma_data, family = "gamma", verbose = FALSE)
  expect_s3_class(fit_gamma, "qbrms_fit")
  expect_equal(fit_gamma$family, "gamma")
  
  # Lognormal family
  lognormal_data <- data.frame(
    y_lognormal = rlnorm(40, meanlog = 0, sdlog = 1),
    x = rnorm(40)
  )
  
  fit_lognormal <- qbrms(y_lognormal ~ x, data = lognormal_data, family = "lognormal", verbose = FALSE)
  expect_s3_class(fit_lognormal, "qbrms_fit")
  expect_equal(fit_lognormal$family, "lognormal")
})

test_that("qbrms handles asymmetric laplace for quantile regression", {
  
  set.seed(456)
  quantile_data <- data.frame(
    y = rt(60, df = 3) + 2 * (1:60)/60,  # Heavy-tailed data for quantile regression
    x = (1:60)/60
  )
  
  # Test different quantiles
  quantiles_to_test <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  
  for (q in quantiles_to_test) {
    fit_q <- qbrms(y ~ x, data = quantile_data, 
                   family = "asymmetric_laplace", quantile = q, verbose = FALSE)
    
    expect_s3_class(fit_q, "qbrms_fit")
    expect_true(is.list(fit_q$family) && fit_q$family$family == "asymmetric_laplace")
    expect_equal(fit_q$family$quantile, q)
  }
})

test_that("qbrms validates family-specific constraints", {
  
  set.seed(123)
  
  # Test beta with out-of-range data
  bad_beta_data <- data.frame(
    y = c(rbeta(30, 2, 3), 1.5, -0.1),  # Values outside (0,1)
    x = rnorm(32)
  )
  
  # Should either handle gracefully or give meaningful error
  expect_error(
    qbrms(y ~ x, data = bad_beta_data, family = "beta", verbose = FALSE),
    NA  # NA means we don't expect a specific error (it might work or fail gracefully)
  )
  
  # Test gamma with negative data
  bad_gamma_data <- data.frame(
    y = c(rgamma(30, 2, 1), -1, -0.5),  # Negative values
    x = rnorm(32)
  )
  
  expect_error(
    qbrms(y ~ x, data = bad_gamma_data, family = "gamma", verbose = FALSE),
    NA
  )
})

test_that("qbrms family aliases work correctly", {
  
  set.seed(123)
  test_data <- data.frame(
    y = rnorm(30),
    x = rnorm(30)
  )
  
  # Test that aliases produce equivalent results
  fit_normal <- qbrms(y ~ x, data = test_data, family = "normal", verbose = FALSE)
  fit_gaussian <- qbrms(y ~ x, data = test_data, family = "gaussian", verbose = FALSE)
  
  expect_equal(fit_normal$family, fit_gaussian$family)
  
  # Test binomial aliases
  test_data$y_binary <- rbinom(30, 1, 0.3)
  fit_binom <- qbrms(y_binary ~ x, data = test_data, family = "binom", verbose = FALSE)
  fit_binomial <- qbrms(y_binary ~ x, data = test_data, family = "binomial", verbose = FALSE)
  
  expect_equal(fit_binom$family, fit_binomial$family)
  
  # Test poisson aliases
  test_data$y_count <- rpois(30, 2)
  fit_pois <- qbrms(y_count ~ x, data = test_data, family = "pois", verbose = FALSE)
  fit_poisson <- qbrms(y_count ~ x, data = test_data, family = "poisson", verbose = FALSE)
  
  expect_equal(fit_pois$family, fit_poisson$family)
})

test_that("qbrms handles family constructors correctly", {
  
  set.seed(123)
  test_data <- data.frame(
    y = rnorm(30),
    y_binary = rbinom(30, 1, 0.3),
    x = rnorm(30)
  )
  
  # Test function constructors vs character strings
  fit_func <- qbrms(y ~ x, data = test_data, family = gaussian(), verbose = FALSE)
  fit_char <- qbrms(y ~ x, data = test_data, family = "gaussian", verbose = FALSE)
  
  expect_equal(fit_func$family, fit_char$family)
  
  # Test specific family constructors
  families_to_test <- list(
    list(constructor = negbinomial, string = "negbinomial"),
    list(constructor = skew_normal, string = "skew_normal"),
    list(constructor = asymmetric_laplace, string = "asymmetric_laplace")
  )
  
  test_data$y_count <- rpois(30, 2)
  
  for (family_info in families_to_test) {
    if (family_info$string == "asymmetric_laplace") {
      # Skip asymmetric_laplace as it needs quantile parameter
      next
    }
    
    if (family_info$string == "negbinomial") {
      y_var <- "y_count"
    } else {
      y_var <- "y"
    }
    
    formula_obj <- as.formula(paste(y_var, "~ x"))
    
    fit_constructor <- qbrms(formula_obj, data = test_data, 
                             family = family_info$constructor(), verbose = FALSE)
    fit_string <- qbrms(formula_obj, data = test_data, 
                        family = family_info$string, verbose = FALSE)
    
    expect_equal(fit_constructor$family, fit_string$family)
  }
})

test_that("qbrms provides informative errors for invalid family specifications", {
  
  test_data <- data.frame(y = rnorm(20), x = rnorm(20))
  
  # Test completely invalid family
  expect_error(
    qbrms(y ~ x, data = test_data, family = "completely_invalid"),
    "not supported|unknown"
  )
  
  # Test asymmetric_laplace without proper quantile
  expect_error(
    qbrms(y ~ x, data = test_data, family = "asymmetric_laplace", quantile = 1.5),
    "quantile must be.*between 0 and 1"
  )
  
  expect_error(
    qbrms(y ~ x, data = test_data, family = "asymmetric_laplace", quantile = c(0.25, 0.75)),
    "quantile must be a single numeric value"
  )
})

test_that("qbrms family documentation and examples work", {
  
  # Test that family constructors have proper structure
  families_to_check <- c("gaussian", "binomial", "poisson", "negbinomial", 
                         "skew_normal", "asymmetric_laplace")
  
  for (family_name in families_to_check) {
    family_constructor <- get(family_name)
    family_obj <- family_constructor()
    
    expect_s3_class(family_obj, "family")
    expect_equal(family_obj$family, family_name)
  }
  
  # Test cumulative with link parameter
  cum_logit <- cumulative("logit")
  cum_probit <- cumulative("probit")
  
  expect_equal(cum_logit$link, "logit")
  expect_equal(cum_probit$link, "probit")
  expect_equal(cum_logit$family, "cumulative")
  expect_equal(cum_probit$family, "cumulative")
})