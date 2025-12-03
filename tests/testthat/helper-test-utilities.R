# =============================================================================
# FILE: tests/testthat/helper-test-utilities.R
# =============================================================================

#' Skip if INLA is not available or broken
#' @keywords internal
skip_if_no_inla <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    testthat::skip("INLA package not installed")
  }
  inla_works <- tryCatch({
    INLA::inla.setOption(num.threads = 1) 
    TRUE
  }, error = function(e) FALSE)
  
  if (!inla_works) {
    testthat::skip("INLA installed but binary executable is missing/broken")
  }
}

#' Create test fit object for use in tests
#' @return A qbrms_fit object
#' @keywords internal
create_test_fit <- function() {
  skip_if_no_inla()
  test_data <- data.frame(
    y = rnorm(50, 2, 1),
    x1 = (1:50)/50,
    x2 = rnorm(50),
    group = factor(rep(1:5, length.out = 50))
  )
  suppressMessages(suppressWarnings(
    qbrms(y ~ x1 + x2, data = test_data, family = gaussian(), verbose = FALSE)
  ))
}