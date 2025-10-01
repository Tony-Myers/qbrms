#' Diagnose Binomial Mixed Effects Models
#'
#' @description
#' Diagnose potential issues in binomial mixed effects models before fitting
#'
#' @param formula Model formula with mixed effects
#' @param data Data frame containing variables
#' @param verbose Logical; print diagnostic information (default: TRUE)
#'
#' @return List with diagnostic information
#' @export
diagnose_binomial_mixed <- function(formula, data, verbose = TRUE) {
  
  if (verbose) cat("=== Binomial Mixed Effects Diagnostics ===\n")
  
  # Use internal diagnostic function that already exists
  result <- .diagnose_binomial_issues(formula, data, verbose)
  
  if (verbose) {
    if (result$has_issues) {
      cat("\nISSUES DETECTED:\n")
      for (issue in result$issues) {
        cat("  * ", issue, "\n")
      }
      cat("\nRECOMMENDATIONS:\n")
      cat("  * Consider using qbrms_binomial_regularised()\n")
      cat("  * Check group sizes and balance\n")
      cat("  * Consider data augmentation if needed\n")
    } else {
      cat("OK: No major issues detected\n")
    }
  }
  
  invisible(result)
}
