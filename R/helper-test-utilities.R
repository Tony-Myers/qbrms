# =============================================================================
# FILE: tests/testthat/helper-test-utilities.R
# =============================================================================
# Test helper functions that are shared across multiple test files

#' Create test fit object for use in tests
#' @keywords internal
create_test_fit <- function() {
  set.seed(123)
  test_data <- data.frame(
    y = rnorm(50, 2 + 0.5 * (1:50)/50, 1),
    x1 = (1:50)/50,
    x2 = rnorm(50),
    group = factor(rep(1:5, length.out = 50))
  )
  
  # Suppress messages/warnings during testing
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