# =============================================================================
# Additional features: model comparison, diagnostics, visualisation, export
# =============================================================================

test_that("compare_models works correctly", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  
  set.seed(123)
  test_data <- data.frame(
    y  = rnorm(80),
    x1 = rnorm(80),
    x2 = rnorm(80)
  )
  
  fit1 <- qbrms(y ~ x1,      data = test_data, family = gaussian())
  fit2 <- qbrms(y ~ x1 + x2, data = test_data, family = gaussian())
  
  comparison <- compare_models(
    Model1    = fit1,
    Model2    = fit2,
    criterion = "waic"
  )
  
  # Basic structure checks
  expect_s3_class(comparison, "qbrms_comparison")
  expect_true("comparison_table" %in% names(comparison))
  expect_true(is.data.frame(comparison$comparison_table))
  expect_equal(nrow(comparison$comparison_table), 2L)
  
  # best_model is optional in the current implementation:
  # if present, it must be one of the provided model names.
  if (!is.null(comparison$best_model)) {
    expect_true(comparison$best_model %in% c("Model1", "Model2"))
  }
})

test_that("visualise_prior handles different distributions", {
  skip_on_cran()
  
  expect_no_error(visualise_prior("normal(0, 10)", add_reference = FALSE))
  expect_no_error(visualise_prior("student_t(3, 0, 5)", add_reference = FALSE))
})

test_that("diagnose_model runs without error for a simple Gaussian model", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  
  set.seed(456)
  dat <- data.frame(
    y = rnorm(60),
    x = rnorm(60)
  )
  
  fit <- qbrms(y ~ x, data = dat, family = gaussian())
  
  # Core requirement: it should run without error
  expect_no_error({
    diag_res <- diagnose_model(fit, verbose = FALSE)
  })
  
  # Light structural check, without assuming specific names
  diag_res <- diagnose_model(fit, verbose = FALSE)
  expect_true(is.list(diag_res))
  expect_true(length(diag_res) >= 1)
})

test_that("export_model works for different formats without error", {
  skip_on_cran()
  
  set.seed(789)
  sim_dat <- data.frame(
    y  = rnorm(20),
    x1 = rnorm(20)
  )
  
  # If export_model expects a model specification, use qbrm();
  # if it expects a fitted model, change this to qbrms(...) as needed.
  spec <- qbrm(
    formula = y ~ x1,
    data    = sim_dat,
    family  = gaussian()
  )
  
  temp_r   <- tempfile(fileext = ".R")
  temp_txt <- tempfile(fileext = ".txt")
  
  on.exit(unlink(c(temp_r, temp_txt)), add = TRUE)
  
  expect_no_error(export_model(spec, temp_r,  format = "R"))
  expect_true(file.exists(temp_r))
  
  expect_no_error(export_model(spec, temp_txt, format = "text"))
  expect_true(file.exists(temp_txt))
})