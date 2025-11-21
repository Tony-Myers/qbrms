# =============================================================================
# FILE: tests/testthat/helper-test-utilities.R
# =============================================================================

#' Skip if INLA is not available or broken
#' @keywords internal
skip_if_no_inla <- function() {
  # 1. Check if installed
  if (!requireNamespace("INLA", quietly = TRUE)) {
    testthat::skip("INLA package not installed")
  }
  
  # 2. Check if functional (binary check)
  # We try a lightweight operation to see if the internal binary responds
  inla_works <- tryCatch({
    # This call usually triggers an internal check
    INLA::inla.setOption(num.threads = 1) 
    TRUE
  }, error = function(e) FALSE)
  
  if (!inla_works) {
    testthat::skip("INLA installed but binary executable is missing/broken")
  }
}

#' Create test fit object for use in tests
#' @keywords internal
create_test_fit <- function() {
  # CRITICAL: Must skip if INLA is broken to avoid crash
  skip_if_no_inla() 
  
  set.seed(123)
  test_data <- data.frame(
    y = rnorm(50, 2 + 0.5 * (1:50)/50, 1),
    x1 = (1:50)/50,
    x2 = rnorm(50),
    group = factor(rep(1:5, length.out = 50))
  )
  
  suppressMessages(suppressWarnings(
    qbrms(y ~ x1 + x2, data = test_data, family = gaussian(), verbose = FALSE)
  ))
}

#' Create test data for conditional effects testing
#' @keywords internal
create_conditional_test_data <- function(n = 60, seed = 123) {
  set.seed(seed)
  data.frame(
    y = rnorm(n, 2 + 0.5 * (1:n)/n, 1),
    x_cont = (1:n)/n,
    x_cat = factor(rep(c("A", "B", "C"), length.out = n)),
    x_int = rnorm(n),
    stringsAsFactors = FALSE
  )
}

#' Skip test if package not available
#' @keywords internal
skip_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    testthat::skip(paste("Package", pkg, "not available"))
  }
}