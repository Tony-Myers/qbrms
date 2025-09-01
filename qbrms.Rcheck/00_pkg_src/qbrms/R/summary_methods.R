# =============================================================================
# R/summary_methods.R
# =============================================================================

#' Summary Method for qbrms Models
#'
#' @description
#' Provides brms-style summary output for qbrms model fits.
#'
#' @param object A `qbrms_fit` model object.
#' @param ... Additional arguments (currently ignored).
#' @return The object invisibly, after printing summary.
#'
#' @export
#' @method summary qbrms_fit
summary.qbrms_fit <- function(object, ...) {
  
  # Model information
  cat("Formula:", deparse(object$original_formula), "\n")
  cat("Data: ", nrow(object$data), " observations\n")
  cat("Family:", extract_family_name(object$family), "\n")
  
  # Timing information
  if (!is.null(object$timing)) {
    cat("Runtime:", object$timing$formatted_duration, "\n")
  } else if (!is.null(object$fitting_time)) {
    cat("Runtime:", format_duration(object$fitting_time), "\n")
  }
  cat("\n")
  
  # Handle different model types
  if (object$model_type == "quantile_regression") {
    cat("Quantile Regression (tau =", object$quantile, ")\n\n")
    print(object$fit$summary.fixed)
    
  } else if (inherits(object, "ordinal_qbrms_fit")) {
    cat("Ordinal Regression Model\n")
    if (!is.null(object$ordinal_levels)) {
      cat("Response levels:", paste(object$ordinal_levels, collapse = " < "), "\n")
    }
    cat("\n")
    
    # For binary decomposition models
    if (inherits(object, "ordinal_binary_qbrms_fit")) {
      cat("Binary Decomposition Approach\n")
      cat("Number of threshold models:", length(object$binary_models), "\n\n")
      for (i in 1:length(object$binary_models)) {
        cat("Threshold", i, ":\n")
        print(object$binary_models[[i]]$summary.fixed)
        cat("\n")
      }
    } else {
      # Augmented data approach
      print(object$fit$summary.fixed)
    }
    
  } else {
    # Standard INLA models
    cat("Population-Level Effects:\n")
    print(object$fit$summary.fixed)
    cat("\n")
    
    # Random effects
    if (object$has_random_effects && !is.null(object$group_var)) {
      cat("Group-Level Effects:\n")
      cat("~", object$group_var, "\n")
      re_summary <- get_random_effects_sd_summary(object$fit, object$group_var)
      if (!is.null(re_summary)) {
        cat("     Estimate Est.Error l-95% CI u-95% CI\n")
        cat("sd   ", sprintf("%.4f", re_summary$mean), 
            "  ", sprintf("%.4f", re_summary$sd),
            "  ", sprintf("%.4f", re_summary$q0.025),
            "  ", sprintf("%.4f", re_summary$q0.975), "\n")
      }
      cat("\n")
    }
    
    # Family parameters
    if (!is.null(object$fit$summary.hyperpar) && nrow(object$fit$summary.hyperpar) > 0) {
      cat("Family Specific Parameters:\n")
      print(object$fit$summary.hyperpar)
      cat("\n")
    }
    
    # Model fit statistics
    metrics <- extract_model_metrics(object$fit)
    if (!is.null(metrics$dic)) {
      cat("Model fit statistics:\n")
      cat("DIC: ", round(metrics$dic, 2), "\n")
    }
    if (!is.null(metrics$waic)) {
      cat("WAIC: ", round(metrics$waic, 2), "\n")
    }
  }
  
  invisible(object)
}

#' Print Method for qbrms Models
#'
#' @description
#' Print method for qbrms model objects.
#'
#' @param x A `qbrms_fit` model object.
#' @param ... Additional arguments (currently ignored).
#' @return The object invisibly, after printing.
#'
#' @export
#' @method print qbrms_fit
print.qbrms_fit <- function(x, ...) {
  cat("qbrms Model fit\n")
  cat("Formula:", deparse(x$original_formula), "\n")
  cat("Data: ", nrow(x$data), " observations\n")
  cat("Family:", extract_family_name(x$family), "\n")
  
  if (x$has_random_effects && !is.null(x$group_var)) {
    cat("Random Effects: ", x$group_var, "\n")
  }
  
  if (!is.null(x$timing)) {
    cat("Runtime:", x$timing$formatted_duration, "\n")
  } else if (!is.null(x$fitting_time)) {
    cat("Runtime:", format_duration(x$fitting_time), "\n")
  }
  
  invisible(x)
}

#' Summary Method for Augmented Ordinal Models
#'
#' @description
#' Summary method for ordinal models fitted using data augmentation.
#'
#' @param object An `ordinal_augmented_qbrms_fit` object.
#' @param ... Additional arguments.
#' @return The object invisibly, after printing the summary.
#'
#' @export
#' @method summary ordinal_augmented_qbrms_fit
summary.ordinal_augmented_qbrms_fit <- function(object, ...) {
  cat("\nAugmented Ordinal qbrms Model Summary\n")
  cat("====================================\n\n")
  
  cat("Formula:", deparse(object$original_formula), "\n")
  cat("Response levels:", paste(object$ordinal_levels, collapse = " < "), "\n")
  cat("Data:", nrow(object$data), "observations\n")
  
  if (!is.null(object$timing)) {
    cat("Runtime:", object$timing$formatted_duration, "\n")
  }
  cat("\n")
  
  cat("Population-Level Effects:\n")
  print(object$fit$summary.fixed)
  
  invisible(object)
}

#' Summary Method for Binary Decomposition Ordinal Models
#'
#' @description
#' Summary method for ordinal models fitted using binary decomposition.
#'
#' @param object An `ordinal_binary_qbrms_fit` object.
#' @param ... Additional arguments.
#' @return The object invisibly, after printing summary.
#'
#' @export
#' @method summary ordinal_binary_qbrms_fit
summary.ordinal_binary_qbrms_fit <- function(object, ...) {
  cat("\nBinary Decomposition Ordinal qbrms Model Summary\n")
  cat("===============================================\n\n")
  
  cat("Formula:", deparse(object$original_formula), "\n")
  cat("Response levels:", paste(object$ordinal_levels, collapse = " < "), "\n")
  cat("Number of binary models:", length(object$binary_models), "\n")
  cat("Data:", nrow(object$data), "observations\n")
  
  if (!is.null(object$timing)) {
    cat("Runtime:", object$timing$formatted_duration, "\n")
  }
  cat("\n")
  
  for (i in 1:length(object$binary_models)) {
    cat("Binary Model", i, "(threshold at level", i, "):\n")
    print(object$binary_models[[i]]$summary.fixed)
    cat("\n")
  }
  
  invisible(object)
}

#' Print Method for Augmented Ordinal Models
#'
#' @description
#' Print method for augmented ordinal model objects.
#'
#' @param x An `ordinal_augmented_qbrms_fit` object.
#' @param ... Additional arguments.
#' @return The object invisibly, after printing.
#'
#' @export
#' @method print ordinal_augmented_qbrms_fit
print.ordinal_augmented_qbrms_fit <- function(x, ...) {
  cat("Augmented Ordinal qbrms Model fit\n")
  cat("Formula:", deparse(x$original_formula), "\n")
  cat("Response levels:", paste(x$ordinal_levels, collapse = " < "), "\n")
  cat("Data:", nrow(x$data), "observations\n")
  
  if (!is.null(x$timing)) {
    cat("Runtime:", x$timing$formatted_duration, "\n")
  }
  
  invisible(x)
}

#' Print Method for Binary Decomposition Ordinal Models
#'
#' @description
#' Print method for binary decomposition ordinal model objects.
#'
#' @param x An `ordinal_binary_qbrms_fit` object.
#' @param ... Additional arguments.
#' @return The object invisibly, after printing.
#'
#' @export
#' @method print ordinal_binary_qbrms_fit
print.ordinal_binary_qbrms_fit <- function(x, ...) {
  cat("Binary Decomposition Ordinal qbrms Model fit\n")
  cat("Formula:", deparse(x$original_formula), "\n")
  cat("Response levels:", paste(x$ordinal_levels, collapse = " < "), "\n")
  cat("Binary models:", length(x$binary_models), "\n")
  cat("Data:", nrow(x$data), "observations\n")
  
  if (!is.null(x$timing)) {
    cat("Runtime:", x$timing$formatted_duration, "\n")
  }
  
  invisible(x)
}

#' Summary for General Ordinal qbrms Fits
#'
#' @description
#' Generic summary method for ordinal qbrms fits.
#'
#' @param object An `ordinal_qbrms_fit` object.
#' @param ... Additional arguments.
#'
#' @export
summary.ordinal_qbrms_fit <- function(object, ...) {
  # Dispatch to more specific methods
  if (inherits(object, "ordinal_binary_qbrms_fit")) {
    summary.ordinal_binary_qbrms_fit(object, ...)
  } else if (inherits(object, "ordinal_augmented_qbrms_fit")) {
    summary.ordinal_augmented_qbrms_fit(object, ...)
  } else {
    # Fallback
    summary.qbrms_fit(object, ...)
  }
}