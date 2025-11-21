# =============================================================================
# FILE: tests/testthat/setup.R
# =============================================================================
# Setup file for testthat tests - runs before all tests

# Load required libraries
library(testthat)

# Conditionally load packages for testing with SAFETY CATCH
if (requireNamespace("INLA", quietly = TRUE)) {
  # INLA-specific setup if available
  # We wrap this in tryCatch because on some CRAN/CI runners (especially macOS),
  # the package is installed but the binary is missing, which causes a crash here.
  tryCatch({
    INLA::inla.setOption(inla.mode = "experimental")
    INLA::inla.setOption(num.threads = "1:1") # Single thread for reproducibility
  }, error = function(e) {
    message("INLA found but not functional (binary missing). Skipping configuration.")
  })
}

# Set up testing environment
options(warn = 1) # Print warnings as they occur

# Helper functions available to all tests
create_simple_test_data <- function(n = 50, seed = 123) {
  set.seed(seed)
  data.frame(
    y_continuous = rnorm(n, 2 + 0.5 * (1:n)/n, 1),
    y_binary = rbinom(n, 1, plogis(-1 + 2 * (1:n)/n)),
    y_count = rpois(n, exp(0.5 + 0.3 * (1:n)/n)),
    x_continuous = (1:n)/n,
    x_factor = factor(rep(LETTERS[1:5], length.out = n)),
    x_numeric = rnorm(n),
    group = factor(rep(1:ceiling(n/10), length.out = n)),
    stringsAsFactors = FALSE
  )
}

# Set global test options
options(
  qbrms.verbose = FALSE, # Suppress verbose output during testing
  width = 80 # Consistent output width
)

# Ensure reproducibility
set.seed(42)

# Clean up function to run after tests
.test_cleanup <- function() {
  # Remove any temporary objects
  if (exists(".test_temp_objects")) {
    rm(list = .test_temp_objects, envir = .GlobalEnv)
  }
  
  # Reset options
  options(warn = 0)
  
  # Garbage collection
  gc()
}

# Register cleanup
reg.finalizer(.GlobalEnv, .test_cleanup, onexit = TRUE)