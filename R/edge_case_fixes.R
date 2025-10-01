# =============================================================================
# R/edge_case_fixes.R - Edge Case Handling and Multinomial Support
# =============================================================================

#' Validate data before model fitting
#'
#' @description
#' Perform lightweight diagnostics on the response and predictors to catch
#' common issues that derail model fitting (missingness, zero variance,
#' impossible values for specific families, simple multicollinearity flags,
#' and sample-size sanity checks).
#'
#' @param formula A model formula.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param family A family object such as \code{gaussian()}, \code{binomial()},
#'   or \code{poisson()}. Used for basic, family-specific checks.
#' @param verbose Logical; print a summary of issues found.
#'
#' @return A list with elements \code{valid} (logical), \code{errors}
#'   (character vector), \code{warnings} (character vector), and
#'   \code{n_complete} (integer count of complete cases across the variables
#'   in the model).
#' @keywords internal
validate_model_data <- function(formula, data, family, verbose = TRUE) {
  
  issues <- list()
  warnings <- character()
  errors <- character()
  
  # Extract variables
  response_var <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  
  # Check 1: Response variable exists and has variation
  if (!(response_var %in% names(data))) {
    errors <- c(errors, paste("Response variable", response_var, "not found in data"))
  } else {
    y <- data[[response_var]]
    
    # Check for all missing
    if (all(is.na(y))) {
      errors <- c(errors, "All values of response variable are missing")
    }
    
    # Check for zero variance
    if (is.numeric(y) && stats::var(y, na.rm = TRUE) < 1e-10) {
      errors <- c(errors, "Response variable has zero or near-zero variance")
    }
    
    # Family-specific checks
    fam_name <- extract_family_name(family)
    
    if (fam_name == "binomial") {
      if (is.numeric(y)) {
        if (!all(y %in% c(0, 1, NA))) {
          errors <- c(errors, "Binomial response must be 0/1 or logical")
        }
        if (all(y == 0, na.rm = TRUE) || all(y == 1, na.rm = TRUE)) {
          errors <- c(errors, "Perfect separation in binomial response (all 0s or all 1s)")
        }
      }
    }
    
    if (fam_name == "poisson") {
      if (is.numeric(y)) {
        if (any(y < 0, na.rm = TRUE)) {
          errors <- c(errors, "Poisson response cannot have negative values")
        }
        if (any(y != round(y), na.rm = TRUE)) {
          warnings <- c(warnings, "Poisson response should be integer counts")
        }
      }
    }
    
    if (fam_name == "multinomial") {
      if (!is.factor(y) && !is.character(y)) {
        warnings <- c(warnings, "Multinomial response should be a factor or character")
      }
      if (is.factor(y) && nlevels(y) < 2) {
        errors <- c(errors, "Multinomial response needs at least 2 categories")
      }
    }
  }
  
  # Check 2: Predictors exist and have variation
  for (pred in predictors) {
    if (!(pred %in% names(data))) {
      errors <- c(errors, paste("Predictor", pred, "not found in data"))
    } else {
      x <- data[[pred]]
      
      if (all(is.na(x))) {
        warnings <- c(warnings, paste("Predictor", pred, "has all missing values"))
      }
      
      if (is.numeric(x) && !all(is.na(x))) {
        if (stats::var(x, na.rm = TRUE) < 1e-10) {
          warnings <- c(warnings, paste("Predictor", pred, "has zero or near-zero variance"))
        }
        
        # Check for extreme values
        x_range <- range(x, na.rm = TRUE)
        if (abs(x_range[1]) > 1e6 || abs(x_range[2]) > 1e6) {
          warnings <- c(warnings, paste("Predictor", pred, "has extreme values"))
        }
      }
      
      if (is.factor(x)) {
        if (nlevels(x) == 1) {
          warnings <- c(warnings, paste("Factor", pred, "has only one level"))
        }
        if (nlevels(x) > 50) {
          warnings <- c(warnings, paste("Factor", pred, "has many levels (", nlevels(x), ")"))
        }
      }
    }
  }
  
  # Check 3: Sample size issues
  n_complete <- sum(stats::complete.cases(data[, c(response_var, predictors), drop = FALSE]))
  n_params <- length(predictors) + 1  # Including intercept
  
  if (n_complete < n_params) {
    errors <- c(errors, paste0("Too few observations (", n_complete, ") for number of parameters (", n_params, ")"))
  } else if (n_complete < 10 * n_params) {
    warnings <- c(warnings, paste("Low observation to parameter ratio:", n_complete / n_params))
  }
  
  # Check 4: Multicollinearity (simple check)
  if (length(predictors) > 1) {
    numeric_preds <- predictors[sapply(predictors, function(p) is.numeric(data[[p]]))]
    if (length(numeric_preds) > 1) {
      cor_matrix <- stats::cor(data[, numeric_preds, drop = FALSE], use = "complete.obs")
      high_cor <- which(abs(cor_matrix) > 0.95 & cor_matrix != 1, arr.ind = TRUE)
      if (nrow(high_cor) > 0) {
        for (i in 1:nrow(high_cor)) {
          if (high_cor[i, 1] < high_cor[i, 2]) {  # Avoid duplicates
            warnings <- c(
              warnings,
              paste(
                "High correlation between",
                numeric_preds[high_cor[i, 1]],
                "and",
                numeric_preds[high_cor[i, 2]]
              )
            )
          }
        }
      }
    }
  }
  
  # Print results
  if (verbose) {
    if (length(errors) > 0) {
      cat("ERRORS found in data validation:\n")
      for (e in errors) cat("  - ", e, "\n", sep = "")
    }
    
    if (length(warnings) > 0) {
      cat("Warnings from data validation:\n")
      for (w in warnings) cat("  - ", w, "\n", sep = "")
    }
    
    if (length(errors) == 0 && length(warnings) == 0) {
      cat("Data validation passed\n")
    }
  }
  
  # Return validation results
  list(
    valid = length(errors) == 0,
    errors = errors,
    warnings = warnings,
    n_complete = n_complete
  )
}

#' Coefficients for multinomial qbrms fits
#'
#' @description
#' Extract a concatenated vector of coefficients from a
#' \code{qbrms_multinomial_fit}, combining the per-category binary submodels
#' if present.
#'
#' @param object A \code{qbrms_multinomial_fit}.
#' @param ... Unused.
#'
#' @return A named numeric vector of coefficients. If coefficient
#'   information is not available, a minimal intercept-only vector is returned.
#' @method coef qbrms_multinomial_fit
#' @export
coef.qbrms_multinomial_fit <- function(object, ...) {
  if (!is.null(object$fit$models)) {
    # Extract coefficients from each binary model
    all_coefs <- list()
    
    for (cat_name in names(object$fit$models)) {
      model <- object$fit$models[[cat_name]]
      if (!is.null(model$summary.fixed)) {
        coefs <- model$summary.fixed[, "mean"]
        names(coefs) <- paste0(cat_name, ":", rownames(model$summary.fixed))
        all_coefs[[cat_name]] <- coefs
      }
    }
    
    # Combine into single vector
    return(unlist(all_coefs))
  }
  
  # Fallback
  return(c("(Intercept)" = 0))
}

#' Summary method for multinomial qbrms fits
#'
#' @description
#' Print a readable summary of a \code{qbrms_multinomial_fit}, including its
#' reference category, the list of categories, and per-category fixed-effect
#' summaries when available.
#'
#' @param object A \code{qbrms_multinomial_fit}.
#' @param digits Integer; number of decimal places to display.
#' @param ... Unused.
#'
#' @return The input \code{object}, returned invisibly.
#' @method summary qbrms_multinomial_fit
#' @export
summary.qbrms_multinomial_fit <- function(object, digits = 2, ...) {
  cat("Multinomial qbrms Model\n\n")
  
  cat("Reference category: ", object$fit$reference_category, "\n", sep = "")
  cat("Categories: ", paste(object$fit$categories, collapse = ", "), "\n\n", sep = "")
  
  if (!is.null(object$fit$models)) {
    for (cat_name in names(object$fit$models)) {
      cat("\n--- Model for category: ", cat_name, " ---\n", sep = "")
      model <- object$fit$models[[cat_name]]
      
      if (!is.null(model$summary.fixed)) {
        formatted_summary <- format_numeric_df(model$summary.fixed, digits = digits)
        print(formatted_summary, quote = FALSE, right = TRUE)
      } else {
        cat("  (No coefficient information available)\n")
      }
    }
  }
  
  cat("\nNumber of observations:", nrow(object$data), "\n")
  
  invisible(object)
}

#' Print method for multinomial qbrms fits
#'
#' @description
#' Shorthand print method that delegates to
#' \code{summary.qbrms_multinomial_fit()}.
#'
#' @param x A \code{qbrms_multinomial_fit}.
#' @param digits Integer; number of decimal places to display.
#' @param ... Unused.
#'
#' @return The input \code{x}, returned invisibly.
#' @method print qbrms_multinomial_fit
#' @export
print.qbrms_multinomial_fit <- function(x, digits = 2, ...) {
  summary.qbrms_multinomial_fit(x, digits = digits, ...)
}

#' Enhanced prior predictive sampling with edge-case handling
#'
#' @description
#' A defensive variant of the prior-predictive sampler that tolerates
#' problematic design matrices and extreme prior draws. It attempts to build a
#' design matrix, samples coefficient vectors from user-supplied prior
#' specifications (if available), and then simulates prior-predictive outcomes
#' according to the chosen family.
#'
#' @param formula A model formula.
#' @param data A \code{data.frame} with the variables used in \code{formula}.
#' @param family A family object, for example \code{gaussian()}, \code{binomial()},
#'   or \code{poisson()}.
#' @param prior A prior specification understood by
#'   \code{extract_prior_specifications()} and \code{sample_from_prior_spec()}.
#'   If \code{NULL}, weak normal priors are used internally.
#' @param ndraws Number of prior draws to simulate.
#' @param verbose Logical; print progress information.
#'
#' @return A numeric matrix of simulated outcomes with shape \code{ndraws x nrow(data)}.
#' @keywords internal
#' @noRd
generate_prior_predictive_samples_robust <- function(formula, data, family = gaussian(),
                                                     prior = NULL, ndraws = 1000,
                                                     verbose = TRUE) {
  
  # Validate data first
  validation <- validate_model_data(formula, data, family, verbose = FALSE)
  if (!validation$valid) {
    stop("Data validation failed:\n", paste(validation$errors, collapse = "\n"))
  }
  
  # Warn about issues
  if (length(validation$warnings) > 0 && verbose) {
    cat("Warnings during prior sampling:\n")
    for (w in validation$warnings) cat("  - ", w, "\n", sep = "")
  }
  
  if (verbose) cat("Generating ", ndraws, " prior predictive samples...\n", sep = "")
  
  response_var <- all.vars(formula)[1]
  n_obs <- nrow(data)
  
  # Generate design matrix with error handling
  X <- tryCatch({
    stats::model.matrix(formula, data)
  }, error = function(e) {
    if (verbose) cat("Warning: Could not create full design matrix. Using simplified version.\n")
    # Create minimal design matrix
    predictors <- all.vars(formula)[-1]
    X_simple <- matrix(1, nrow = n_obs, ncol = 1)  # Intercept
    colnames(X_simple) <- "(Intercept)"
    
    for (pred in predictors) {
      if (pred %in% names(data) && is.numeric(data[[pred]])) {
        X_simple <- cbind(X_simple, data[[pred]])
        colnames(X_simple)[ncol(X_simple)] <- pred
      }
    }
    X_simple
  })
  
  n_coefs <- ncol(X)
  coef_names <- colnames(X)
  
  # Extract prior specifications
  prior_specs <- extract_prior_specifications(prior, coef_names, verbose)
  
  # Generate prior samples for coefficients with bounds checking
  prior_coef_samples <- matrix(NA_real_, nrow = ndraws, ncol = n_coefs)
  colnames(prior_coef_samples) <- coef_names
  
  for (i in seq_len(ndraws)) {
    coef_sample <- numeric(n_coefs)
    names(coef_sample) <- coef_names
    
    for (j in seq_len(n_coefs)) {
      coef_name <- coef_names[j]
      prior_spec <- prior_specs[[coef_name]]
      
      # Sample with bounds
      sample_val <- sample_from_prior_spec(prior_spec)
      
      # Bound extreme values to prevent numerical issues
      sample_val <- pmax(pmin(sample_val, 100), -100)
      coef_sample[j] <- sample_val
    }
    
    prior_coef_samples[i, ] <- coef_sample
  }
  
  # Generate predictions from prior coefficient samples
  yrep <- matrix(NA_real_, nrow = ndraws, ncol = n_obs)
  family_name <- extract_family_name(family)
  
  for (i in seq_len(ndraws)) {
    linear_pred <- as.numeric(X %*% prior_coef_samples[i, ])
    
    # Bound linear predictor to prevent numerical issues
    linear_pred <- pmax(pmin(linear_pred, 20), -20)
    
    # Generate observations based on family with error handling
    tryCatch({
      if (family_name %in% c("gaussian", "normal")) {
        error_sd <- abs(stats::rnorm(1, 0, 1))
        if (error_sd < 0.01) error_sd <- 0.01  # Minimum SD
        if (error_sd > 100) error_sd <- 100    # Maximum SD
        yrep[i, ] <- stats::rnorm(n_obs, linear_pred, error_sd)
        
      } else if (family_name == "binomial") {
        probs <- stats::plogis(linear_pred)
        yrep[i, ] <- stats::rbinom(n_obs, 1, probs)
        
      } else if (family_name == "poisson") {
        lambda <- exp(linear_pred)
        lambda <- pmin(lambda, 1e6)  # Cap at reasonable maximum
        yrep[i, ] <- stats::rpois(n_obs, lambda)
        
      } else if (family_name == "multinomial") {
        # For multinomial, generate categorical outcomes
        n_cats <- length(unique(data[[response_var]]))
        if (n_cats < 2) n_cats <- 3  # Default to 3 categories
        yrep[i, ] <- sample(seq_len(n_cats), n_obs, replace = TRUE)
        
      } else {
        # Default fallback
        yrep[i, ] <- stats::rnorm(n_obs, linear_pred, 1)
      }
    }, error = function(e) {
      if (verbose) cat("Warning: Prior sample ", i, " failed. Using fallback.\n", sep = "")
      yrep[i, ] <<- stats::rnorm(n_obs, 0, 1)
    })
  }
  
  if (verbose) {
    cat("Prior predictive samples generated:\n")
    rng <- range(yrep, na.rm = TRUE)
    cat("  Sample range: ", round(rng[1], 2), " to ", round(rng[2], 2), "\n", sep = "")
    cat("  Sample mean: ", round(mean(yrep, na.rm = TRUE), 2), "\n", sep = "")
  }
  
  return(yrep)
}

#' Fallback model fitting for edge cases
#'
#' @description
#' Very simple fitting strategies used as a last resort when primary fitting
#' fails. These provide coefficient means and approximate standard deviations
#' sufficient for downstream summaries and plotting.
#'
#' @param formula A model formula.
#' @param data A \code{data.frame}.
#' @param family A family object.
#' @param verbose Logical; print progress information.
#'
#' @return A list with a \code{summary.fixed} data frame (columns \code{mean}
#'   and \code{sd}) and small metadata flags.
#' @keywords internal
fit_fallback_model <- function(formula, data, family, verbose = TRUE) {
  if (verbose) cat("Using simplified fallback model fitting...\n")
  
  response_var <- all.vars(formula)[1]
  y <- data[[response_var]]
  
  # Create very simple model depending on family
  family_name <- extract_family_name(family)
  
  if (family_name == "gaussian") {
    # Use lm as fallback
    lm_fit <- stats::lm(formula, data = data)
    
    # Convert to INLA-like structure
    summary_fixed <- data.frame(
      mean = stats::coef(lm_fit),
      sd = sqrt(diag(stats::vcov(lm_fit))),
      row.names = names(stats::coef(lm_fit))
    )
    
    fit <- list(
      summary.fixed = summary_fixed,
      converged = TRUE,
      fallback = TRUE
    )
    
  } else if (family_name == "binomial") {
    # Use glm as fallback
    glm_fit <- stats::glm(formula, data = data, family = stats::binomial())
    
    summary_fixed <- data.frame(
      mean = stats::coef(glm_fit),
      sd = sqrt(diag(stats::vcov(glm_fit))),
      row.names = names(stats::coef(glm_fit))
    )
    
    fit <- list(
      summary.fixed = summary_fixed,
      converged = isTRUE(glm_fit$converged),
      fallback = TRUE
    )
    
  } else {
    # Very basic fallback - just return means
    fit <- list(
      summary.fixed = data.frame(
        mean = mean(y, na.rm = TRUE),
        sd = stats::sd(y, na.rm = TRUE),
        row.names = "(Intercept)"
      ),
      converged = FALSE,
      fallback = TRUE
    )
  }
  
  return(fit)
}

#' Safe construction of model matrices
#'
#' @description
#' Attempt to build \code{model.matrix()} with additional checks for
#' size and rank, and fall back to an intercept-only matrix on failure.
#'
#' @param formula A model formula.
#' @param data A \code{data.frame}.
#' @param max_cols Integer; warn if more than this many columns are produced.
#'
#' @return A numeric matrix suitable for downstream computations.
#' @keywords internal
safe_model_matrix <- function(formula, data, max_cols = 100) {
  X <- tryCatch({
    mm <- stats::model.matrix(formula, data)
    
    # Check for excessive columns (likely from high-cardinality factors)
    if (ncol(mm) > max_cols) {
      warning("Model matrix has ", ncol(mm), " columns. Consider reducing factor levels.")
    }
    
    # Check for rank deficiency
    qr_X <- qr(mm)
    if (qr_X$rank < ncol(mm)) {
      warning("Model matrix is rank deficient. Some coefficients may not be estimable.")
    }
    
    mm
  }, error = function(e) {
    warning("Could not create model matrix: ", e$message)
    # Return minimal matrix
    matrix(1, nrow = nrow(data), ncol = 1,
           dimnames = list(NULL, "(Intercept)"))
  })
  
  return(X)
}
