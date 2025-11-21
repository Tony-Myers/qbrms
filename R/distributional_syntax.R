# =============================================================================
# R/distributional_syntax.R
# =============================================================================

#' Create a Bayesian Formula
#'
#' @description
#' Function to set up a model formula for use in \code{qbrms}, allowing
#' specification of distributional parameters (e.g., sigma) in addition to the
#' mean structure.
#'
#' @param formula Main model formula (for the mean/location parameter).
#' @param ... Additional formulas for distributional parameters (e.g., \code{sigma ~ x}).
#' @param flist Optional list of formulas (for internal use).
#' @param family Same as in \code{qbrms()} (optional here).
#' @param nl Logical; indicating if the model is non-linear (not yet fully supported).
#'
#' @return An object of class \code{brmsformula} (and \code{qbrmsformula})
#'   containing the parsed formulas.
#'
#' @details
#' This function mimics the \code{brms::bf()} syntax to allow users familiar
#' with brms to define distributional models.
#'
#' Supported distributional parameters depend on the family:
#' \itemize{
#'   \item \code{gaussian}: \code{sigma} (residual standard deviation)
#'   \item \code{student_t}: \code{sigma}, \code{nu} (degrees of freedom)
#'   \item \code{lognormal}: \code{sigma} (shape parameter)
#'   \item \code{beta}: \code{phi} (precision)
#'   \item \code{simplex}: \code{phi} (precision)
#' }
#'
#' @examples
#' \dontrun{
#' # Standard model
#' f1 <- bf(y ~ x)
#'
#' # Distributional model (heteroscedasticity)
#' # Sigma varies by group
#' f2 <- bf(y ~ x, sigma ~ group)
#' }
#'
#' @export
bf <- function(formula, ..., flist = NULL, family = NULL, nl = FALSE) {
  
  # Start with the main formula
  out <- list(formula = as.formula(formula))
  
  # Capture distributional formulas from ...
  dpars <- list(...)
  
  # Validate dpars are formulas
  if (length(dpars) > 0) {
    is_formula <- vapply(dpars, function(x) inherits(x, "formula"), logical(1))
    if (!all(is_formula)) {
      stop("All additional arguments to bf() must be formulas (e.g., sigma ~ x)", call. = FALSE)
    }
    
    # Extract parameter names from LHS of formulas (e.g., "sigma" from sigma ~ x)
    dpar_names <- vapply(dpars, function(x) {
      as.character(x[[2]])
    }, character(1))
    
    names(dpars) <- dpar_names
    
    # Clean formulas (remove LHS for storage)
    for (nm in names(dpars)) {
      # We keep the full formula for now to allow parsing later
      dpars[[nm]] <- list(formula = dpars[[nm]])
    }
  }
  
  # Add to output object
  out$dpars <- dpars
  out$nl <- nl
  out$family <- family
  
  # Assign classes for compatibility with your existing parser
  class(out) <- c("qbrmsformula", "brmsformula", "list")
  
  return(out)
}

#' Print method for qbrms formulas
#' @param x A qbrmsformula object
#' @param ... Unused
#' @export
print.qbrmsformula <- function(x, ...) {
  cat("qbrms Formula:\n")
  print(x$formula)
  
  if (length(x$dpars) > 0) {
    cat("\nDistributional Parameters:\n")
    for (nm in names(x$dpars)) {
      cat(paste0("  ", nm, ": "))
      # Extract RHS for display
      f <- x$dpars[[nm]]$formula
      cat(paste(deparse(f[[3]]), collapse = ""), "\n")
    }
  }
  invisible(x)
}