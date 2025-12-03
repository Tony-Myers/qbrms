# =============================================================================
#s3_methods.R - S3 Methods for qbrms Objects
# =============================================================================

# =============================================================================
# SECTION 1: PRINT AND SUMMARY METHODS FOR qbrms_fit
# =============================================================================

#' Print Method for qbrms_fit Objects
#'
#' @description
#' Prints a summary of a fitted qbrms model object.
#'
#' @param x A \code{qbrms_fit} object.
#' @param digits Number of decimal places for output (default: 2).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   fit <- qbrms(mpg ~ hp, data = mtcars, family = gaussian(), verbose = FALSE)
#'   print(fit)
#' }
#' }
#'
#' @export
#' @method print qbrms_fit
print.qbrms_fit <- function(x, digits = 2, ...) {
  cat("qbrms Model Fit\n\n")
  
  # Print formula
  if (!is.null(x$original_formula)) {
    cat("Formula:", deparse(x$original_formula), "\n")
  }
  
  # Print data info
  if (!is.null(x$data)) {
    cat("Data:   ", nrow(x$data), "observations\n")
  }
  
  # Print family
  fam <- x$family
  if (is.character(fam)) {
    cat("Family: ", fam, "\n")
  } else if (is.list(fam) && !is.null(fam$family)) {
    cat("Family: ", fam$family, "\n")
  }
  
  cat("\n")
  
  # Print coefficients if available
  if (!is.null(x$fit) && !is.null(x$fit$summary.fixed)) {
    cat("Population-Level Effects:\n")
    summ <- x$fit$summary.fixed
    
    # Format numeric columns
    numeric_cols <- sapply(summ, is.numeric)
    summ[numeric_cols] <- lapply(summ[numeric_cols], 
                                 function(col) round(col, digits))
    print(summ, quote = FALSE, right = TRUE)
  }
  
  invisible(x)
}

#' Summary Method for qbrms_fit Objects
#'
#' @description
#' Provides a detailed summary of a fitted qbrms model.
#'
#' @param object A \code{qbrms_fit} object.
#' @param digits Number of decimal places for output (default: 2).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.qbrms_fit"} containing:
#' \itemize{
#'   \item \code{formula}: The model formula.
#'   \item \code{family}: The distribution family.
#'   \item \code{nobs}: Number of observations.
#'   \item \code{fixed}: Data frame of fixed effects estimates.
#'   \item \code{random}: Random effects summary (if applicable).
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   fit <- qbrms(mpg ~ hp, data = mtcars, family = gaussian(), verbose = FALSE)
#'   summary(fit)
#' }
#' }
#'
#' @export
#' @method summary qbrms_fit
summary.qbrms_fit <- function(object, digits = 2, ...) {
  
  cat("Model: qbrms Bayesian Regression (INLA backend)\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  # Formula
  if (!is.null(object$original_formula)) {
    cat("Formula:", deparse(object$original_formula), "\n")
  }
  
  # Family
  fam <- object$family
  if (is.character(fam)) {
    cat("Family: ", fam, "\n")
  } else if (is.list(fam) && !is.null(fam$family)) {
    cat("Family: ", fam$family, "\n")
  }
  
  # Data
  if (!is.null(object$data)) {
    cat("Observations:", nrow(object$data), "\n")
  }
  
  cat("\n")
  
  # Population-Level Effects
  cat("Population-Level Effects:\n")
  if (!is.null(object$fit) && !is.null(object$fit$summary.fixed)) {
    summ <- object$fit$summary.fixed
    numeric_cols <- sapply(summ, is.numeric)
    summ[numeric_cols] <- lapply(summ[numeric_cols], 
                                 function(col) round(col, digits))
    print(summ, quote = FALSE, right = TRUE)
  } else {
    cat("  (No fixed effects summary available)\
")
  }
  
  cat("\n")
  
  # Group-Level Effects
  cat("Group-Level Effects:\n")
  if (isTRUE(object$has_random_effects) && !is.null(object$fit$summary.random)) {
    for (nm in names(object$fit$summary.random)) {
      cat("  ~", nm, "\n")
      re_summ <- object$fit$summary.random[[nm]]
      if (is.data.frame(re_summ) && nrow(re_summ) > 0) {
        cat("    Groups:", nrow(re_summ), "\n")
      }
    }
  } else {
    cat("  (none)\n")
  }
  
  cat("\n")
  
  # Fitting time
  if (!is.null(object$fitting_time)) {
    cat("Fitting time:", round(object$fitting_time, 2), "seconds\n")
  }
  
  # Create summary object
  result <- list(
    formula = object$original_formula,
    family = object$family,
    nobs = if (!is.null(object$data)) nrow(object$data) else NA,
    fixed = if (!is.null(object$fit)) object$fit$summary.fixed else NULL,
    random = if (!is.null(object$fit)) object$fit$summary.random else NULL
  )
  class(result) <- "summary.qbrms_fit"
  
  invisible(result)
}

#' Print Method for summary.qbrms_fit Objects
#'
#' @param x A \code{summary.qbrms_fit} object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
#' @method print summary.qbrms_fit
print.summary.qbrms_fit <- function(x, ...) {
  cat("qbrms Model Summary\n")
  if (!is.null(x$formula)) {
    cat("Formula:", deparse(x$formula), "\n")
  }
  if (!is.null(x$fixed)) {
    cat("\nFixed Effects:\n")
    print(x$fixed)
  }
  invisible(x)
}

# =============================================================================
# SECTION 2: PRINT AND SUMMARY METHODS FOR qbrmb_fit (Binomial regularised)
# =============================================================================

#' Print Method for qbrmb_fit Objects
#'
#' @description
#' Prints a summary of a regularised binomial qbrms model.
#'
#' @param x A \code{qbrmb_fit} object.
#' @param digits Number of decimal places for output (default: 2).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
#' @method print qbrmb_fit
print.qbrmb_fit <- function(x, digits = 2, ...) {
  cat("qbrmb Regularised Binomial Model Fit\n\n")
  
  if (!is.null(x$original_formula)) {
    cat("Formula:", deparse(x$original_formula), "\n")
  }
  
  if (!is.null(x$data)) {
    cat("Data:   ", nrow(x$data), "observations\n")
  }
  
  cat("Family:  binomial (regularised)\n\n")
  
  if (!is.null(x$fit) && !is.null(x$fit$summary.fixed)) {
    cat("Population-Level Effects:\n")
    summ <- x$fit$summary.fixed
    numeric_cols <- sapply(summ, is.numeric)
    summ[numeric_cols] <- lapply(summ[numeric_cols], 
                                 function(col) round(col, digits))
    print(summ, quote = FALSE, right = TRUE)
  }
  
  invisible(x)
}

#' Summary Method for qbrmb_fit Objects
#'
#' @description
#' Provides a detailed summary of a regularised binomial qbrms model.
#'
#' @param object A \code{qbrmb_fit} object.
#' @param digits Number of decimal places for output (default: 2).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.qbrmb_fit"} containing model summary information.
#'
#' @export
#' @method summary qbrmb_fit
summary.qbrmb_fit <- function(object, digits = 2, ...) {
  
  cat("Model: qbrmb Regularised Binomial Regression\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  if (!is.null(object$original_formula)) {
    cat("Formula:", deparse(object$original_formula), "\n")
  }
  
  cat("Family:  binomial\n")
  
  if (!is.null(object$data)) {
    cat("Observations:", nrow(object$data), "\n")
  }
  
  cat("\n")
  
  cat("Population-Level Effects:\n")
  if (!is.null(object$fit) && !is.null(object$fit$summary.fixed)) {
    summ <- object$fit$summary.fixed
    numeric_cols <- sapply(summ, is.numeric)
    summ[numeric_cols] <- lapply(summ[numeric_cols], 
                                 function(col) round(col, digits))
    print(summ, quote = FALSE, right = TRUE)
  } else {
    cat("  (No fixed effects summary available)\n")
  }
  
  result <- list(
    formula = object$original_formula,
    family = "binomial",
    nobs = if (!is.null(object$data)) nrow(object$data) else NA,
    fixed = if (!is.null(object$fit)) object$fit$summary.fixed else NULL
  )
  class(result) <- "summary.qbrmb_fit"
  
  invisible(result)
}

