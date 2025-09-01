# =============================================================================
# R/data_handling.R
# =============================================================================

#' Handle Missing Data
#'
#' @description
#' Process missing data in model variables, similar to brms handling.
#'
#' @param formula Model formula.
#' @param data Data frame.
#' @param verbose Logical, whether to print messages.
#' @return Data frame with complete cases for model variables.
#'
#' @keywords internal
handle_missing_data <- function(formula, data, verbose = TRUE) {
  if (verbose) cat("Processing missing data...\n")
  
  all_vars <- all.vars(formula)
  
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  
  original_n <- nrow(data)
  complete_data <- data[complete.cases(data[, all_vars, drop = FALSE]), ]
  final_n <- nrow(complete_data)
  
  if (final_n < original_n) {
    if (verbose) {
      cat("Missing data handling:\n")
      cat("  Original:", original_n, "observations\n")
      cat("  Final:", final_n, "observations\n")
      cat("  Removed:", original_n - final_n, "rows with missing values\n")
    }
  } else {
    if (verbose) cat("No missing values detected in model variables\n")
  }
  
  return(complete_data)
}