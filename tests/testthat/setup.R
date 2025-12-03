# =============================================================================
# FILE: tests/testthat/setup.R
# =============================================================================

library(testthat)

# Conditionally load packages for testing with SAFETY CATCH
if (requireNamespace("INLA", quietly = TRUE)) {
  tryCatch({
    INLA::inla.setOption(inla.mode = "experimental")
    INLA::inla.setOption(num.threads = "1:1") 
  }, error = function(e) {
    message("INLA found but not functional (binary missing). Skipping configuration.")
  })
}

# Set up testing environment options
options(
  qbrms.verbose = FALSE, 
  width = 80 
)