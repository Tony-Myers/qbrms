# =============================================================================
# R/utils.R
# =============================================================================

#' Null Coalescing Operator
#'
#' @description
#' Returns the first non-NULL value.
#'
#' @name null-coalesce
#' @aliases %||%
#' @param x First value.
#' @param y Second value.
#' @return x if not NULL, otherwise y.
#'
#' @keywords internal
#' @export
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Format Duration
#'
#' @description
#' Format duration in seconds to a human-readable string.
#'
#' @param seconds Numeric duration in seconds.
#' @return Character string with formatted duration.
#'
#' @keywords internal
format_duration <- function(seconds) {
  # Handle edge cases
  if (is.null(seconds) || length(seconds) == 0 || !is.numeric(seconds)) {
    return("0 seconds")
  }
  
  # Handle NA or infinite values FIRST
  if (is.na(seconds) || !is.finite(seconds)) {
    return("unknown duration")
  }
  
  # Ensure non-negative
  seconds <- abs(seconds)
  
  # FIXED: Use < 60 (not <= 60) so that exactly 60.5 stays in seconds
  if (seconds < 60) {
    return(paste(round(seconds, 2), "seconds"))
  } else if (seconds < 3600) {
    mins <- floor(seconds / 60)
    secs <- seconds %% 60
    return(paste(mins, "minutes", round(secs, 1), "seconds"))
  } else {
    hours <- floor(seconds / 3600)
    mins <- floor((seconds %% 3600) / 60)
    secs <- seconds %% 60
    return(paste(hours, "hours", mins, "minutes", round(secs, 1), "seconds"))
  }
}

if (!exists(".qbrms__remove_random_effects", mode = "function")) {
  .qbrms__remove_random_effects <- function(formula) drop_random_effects(formula)
}

# ---- Fixed-effects design helpers (strip random effects safely) -------------
# Internal only
.qbrms_terms_fixed <- function(object, data) {
  # remove any "(...|...)" random-effect chunks
  f_str <- paste(deparse(object$original_formula), collapse = "")
  f_str <- gsub("\\([^|]*\\|[^)]*\\)", "", f_str)
  # collapse accidental duplicate '+' and trim ends
  f_str <- gsub("\\+\\s*\\+", "+", f_str)
  f_fixed <- stats::as.formula(f_str, env = environment(object$original_formula))
  stats::delete.response(stats::terms(f_fixed, data = data))
}

# Build X for training data (fixed effects only)
.qbrms_model_matrix_fixed <- function(object, data) {
  tt <- .qbrms_terms_fixed(object, data)
  stats::model.matrix(tt, data)
}

# Build X for newdata using training contrasts (fixed effects only)
.qbrms_model_matrix_fixed_newdata <- function(object, newdata) {
  tt     <- .qbrms_terms_fixed(object, object$data)
  Xtrain <- stats::model.matrix(tt, object$data)
  stats::model.matrix(tt, newdata, contrasts.arg = attr(Xtrain, "contrasts"))
}


#' Extract Model Metrics
#'
#' @description
#' Extract DIC, WAIC and other metrics from an INLA fit.
#'
#' @param inla_fit INLA model object.
#' @return A list of metrics.
#'
#' @keywords internal
extract_model_metrics <- function(inla_fit) {
  metrics <- list(
    dic = NULL,
    waic = NULL,
    marginal_likelihood = NULL,
    cpo_failure_count = NULL
  )
  
  if (!is.null(inla_fit$dic) && !is.null(inla_fit$dic$dic)) {
    metrics$dic <- inla_fit$dic$dic
  }
  
  if (!is.null(inla_fit$waic) && !is.null(inla_fit$waic$waic)) {
    metrics$waic <- inla_fit$waic$waic
  }
  
  if (!is.null(inla_fit$mlik)) {
    metrics$marginal_likelihood <- inla_fit$mlik[1, 1]
  }
  
  if (!is.null(inla_fit$cpo) && !is.null(inla_fit$cpo$failure)) {
    metrics$cpo_failure_count <- sum(inla_fit$cpo$failure, na.rm = TRUE)
  }
  
  return(metrics)
}

#' Get Random Effects Standard Deviation Summary
#'
#' @description
#' Extract random effects standard deviation from INLA hyperparameters.
#'
#' @param inla_fit INLA model object.
#' @param group_var Group variable name.
#' @return List with mean, sd, and quantiles of random effects SD.
#'
#' @keywords internal
get_random_effects_sd_summary <- function(inla_fit, group_var) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA package is required but not available. Please install it from: https://www.r-inla.org/download-install")
  }
  
  if (!is.null(inla_fit$marginals.hyperpar)) {
    re_names <- names(inla_fit$marginals.hyperpar)
    re_index <- grep(paste0("Precision for ", group_var), re_names)
    
    if (length(re_index) == 0) {
      re_index <- grep("Precision for group_id", re_names)
    }
    if (length(re_index) == 0) {
      re_index <- grep("Precision for", re_names)[1]
    }
    
    if (length(re_index) > 0 && !is.na(re_index[1])) {
      precision_marginal <- inla_fit$marginals.hyperpar[[re_index[1]]]
      sd_marginal <- INLA::inla.tmarginal(function(x) 1/sqrt(x), precision_marginal)
      sd_summary <- INLA::inla.zmarginal(sd_marginal, silent = TRUE)
      
      return(list(
        mean = sd_summary$mean,
        sd = sd_summary$sd,
        q0.025 = INLA::inla.qmarginal(0.025, sd_marginal),
        q0.975 = INLA::inla.qmarginal(0.975, sd_marginal)
      ))
    }
  }
  
  possible_names <- c(
    paste0("Precision for ", group_var),
    "Precision for group_id"
  )
  
  if (!is.null(inla_fit$summary.hyperpar)) {
    for (precision_row_name in possible_names) {
      if (precision_row_name %in% rownames(inla_fit$summary.hyperpar)) {
        precision_row <- inla_fit$summary.hyperpar[precision_row_name, ]
        
        central_precision <- precision_row$mean
        sd_estimate <- 1/sqrt(central_precision)
        sd_lower <- 1/sqrt(precision_row$`0.975quant`)
        sd_upper <- 1/sqrt(precision_row$`0.025quant`)
        se_approx <- (sd_upper - sd_lower) / 3.92
        
        return(list(
          mean = sd_estimate,
          sd = se_approx,
          q0.025 = sd_lower,
          q0.975 = sd_upper
        ))
      }
    }
  }
  
  return(NULL)
}

#' Extract Family Name from INLA Family Object
#'
#' @description
#' Helper function to extract family name handling both strings and lists.
#'
#' @param inla_family INLA family specification (string or list).
#' @return Character string with family name.
#'
#' @export
extract_family_name <- function(inla_family) {
  if (is.list(inla_family)) {
    return(inla_family$family %||% "gaussian")
  } else if (is.character(inla_family) && length(inla_family) >= 1) {
    return(inla_family[1])
  } else {
    return("gaussian")
  }
}

#' Create Dummy Data for Testing
#'
#' @description
#' Create dummy data that preserves structure for testing purposes.
#'
#' @param formula Model formula.
#' @param data Original data frame.
#' @param n_dummy Number of dummy observations to create.
#' @param family_name The name of the model family (e.g., "gaussian").
#' @param verbose Logical, whether to print messages.
#' @return Data frame with dummy structure.
#'
#' @export
create_dummy_data <- function(formula, data, n_dummy = 10, family_name = "gaussian", verbose = FALSE) {
  if (is.null(formula)) {
    stop("formula is required")
  }
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  if (nrow(data) == 0) {
    # Create minimal dummy data for empty input
    response_var <- all.vars(formula)[1]
    dummy_data <- data.frame(
      y = rnorm(n_dummy),
      row.names = 1:n_dummy
    )
    names(dummy_data)[1] <- response_var
    return(dummy_data)
  }
  
  # Use the existing create_dummy_data_for_priors function as a basis
  create_dummy_data_for_priors(formula, data, n_dummy, family_name, verbose)
}

#' Validate Family Quantile Combination
#'
#' @description
#' Check if a family supports quantile regression with a given quantile value.
#' Throws informative errors for invalid combinations.
#'
#' @param family_name Character string specifying the family name.
#' @param quantile Numeric quantile value (or NULL).
#' @return TRUE if the combination is valid (invisibly), throws error otherwise.
#'
#' @export
validate_family_quantile <- function(family_name, quantile) {
  
  # NULL quantile is always acceptable - no quantile regression requested
  if (is.null(quantile)) {
    return(invisible(TRUE))
  }
  
  # Validate quantile value first (for all families)
  if (!is.numeric(quantile) || length(quantile) != 1) {
    stop("quantile must be a single numeric value between 0 and 1", 
         call. = FALSE)
  }
  
  if (quantile <= 0 || quantile >= 1) {
    stop("quantile must be a single numeric value between 0 and 1", 
         call. = FALSE)
  }
  
  # Only asymmetric_laplace supports quantile regression
  if (family_name == "asymmetric_laplace") {
    return(invisible(TRUE))
  }
  
  # All other families don't support quantile regression - throw informative error
  stop("Family '", family_name, "' does not support quantile regression", 
       call. = FALSE)
}

#' Check if Family Supports Quantile Regression
#'
#' @description
#' Determine whether a given family supports quantile regression.
#'
#' @param family_obj Family object or name.
#' @return Logical indicating whether the family supports quantile regression.
#'
#' @export
family_supports_quantile <- function(family_obj) {
  
  # Extract family name
  if (is.character(family_obj)) {
    family_name <- family_obj
  } else if (is.list(family_obj) && "family" %in% names(family_obj)) {
    family_name <- family_obj$family
  } else if (is.function(family_obj)) {
    # Try to evaluate the function to get family object
    tryCatch({
      fam_obj <- family_obj()
      if (is.list(fam_obj) && "family" %in% names(fam_obj)) {
        family_name <- fam_obj$family
      } else {
        return(FALSE)
      }
    }, error = function(e) return(FALSE))
  } else {
    return(FALSE)
  }
  
  # Only asymmetric_laplace supports quantiles currently
  return(family_name == "asymmetric_laplace")
}

#' Drop Random Effects from Formula
#'
#' @description
#' Remove random effects terms from a model formula.
#'
#' @param formula A model formula that may contain random effects.
#' @return Formula with random effects terms removed.
#'
#' @export
drop_random_effects <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("Input must be a formula")
  }
  
  formula_str <- deparse(formula, width.cutoff = 500L)
  lhs <- sub("~.*", "", formula_str)
  rhs <- sub(".*~", "", formula_str)
  
  # Remove random effects terms in parentheses with |
  rhs_clean <- gsub("\\([^\\)]*\\|[^\\)]*\\)", "", rhs)
  rhs_clean <- gsub("\\+\\s*\\+", "+", rhs_clean)  # Fix double +
  rhs_clean <- gsub("^\\s*\\+\\s*", "", rhs_clean)  # Remove leading +
  rhs_clean <- gsub("\\s*\\+\\s*$", "", rhs_clean)  # Remove trailing +
  
  if (nchar(trimws(rhs_clean)) == 0) {
    rhs_clean <- "1"
  }
  
  as.formula(paste(lhs, "~", rhs_clean))
}

#' Create Dummy Data for Prior Predictive Checks
#'
#' @description
#' Internal function to create dummy data that preserves structure.
#'
#' @param formula Model formula.
#' @param data Original data.
#' @param n_dummy Number of dummy observations.
#' @param family_name The name of the model family (e.g., "poisson").
#' @param verbose Logical, whether to print messages.
#' @return Data frame with dummy structure.
#'
#' @keywords internal
create_dummy_data_for_priors <- function(formula, data, n_dummy, family_name, verbose = TRUE) {
  if (verbose) cat("Creating dummy data with", n_dummy, "observations for prior sampling...\n")
  
  formula_str <- deparse(formula, width.cutoff = 500L)
  has_trials <- grepl("trials\\(", formula_str)
  trials_var <- NULL
  
  if (has_trials) {
    trials_match <- regmatches(formula_str, regexpr("trials\\(([^)]+)\\)", formula_str))
    if (length(trials_match) > 0) {
      trials_var <- gsub("trials\\(|\\)", "", trials_match)
      trials_var <- trimws(trials_var)
    }
    clean_formula_str <- gsub("\\s*\\|\\s*trials\\([^)]+\\)", "", formula_str)
    clean_formula <- as.formula(clean_formula_str)
  } else {
    clean_formula <- formula
  }
  
  all_vars <- all.vars(clean_formula)
  response_var <- all_vars[1]
  predictor_vars <- all_vars[-1]
  
  dummy_data <- data.frame(row.names = 1:n_dummy)
  
  if (response_var %in% names(data)) {
    original_response <- data[[response_var]]
    
    if (is.numeric(original_response)) {
      if (has_trials && !is.null(trials_var)) {
        dummy_data[[response_var]] <- sample(0:3, n_dummy, replace = TRUE)
      } else {
        dummy_data[[response_var]] <- rnorm(n_dummy, 0, 0.1)
      }
    } else if (is.factor(original_response) || is.ordered(original_response)) {
      levels_orig <- levels(original_response)
      dummy_data[[response_var]] <- factor(sample(levels_orig, n_dummy, replace = TRUE),
                                           levels = levels_orig,
                                           ordered = is.ordered(original_response))
    } else {
      dummy_data[[response_var]] <- rnorm(n_dummy, 0, 0.1)
    }
  } else {
    dummy_data[[response_var]] <- rnorm(n_dummy, 0, 0.1)
  }
  
  for (var in predictor_vars) {
    if (var %in% names(data)) {
      original_var <- data[[var]]
      
      if (is.numeric(original_var)) {
        var_mean <- mean(original_var, na.rm = TRUE)
        var_sd <- sd(original_var, na.rm = TRUE)
        dummy_data[[var]] <- rnorm(n_dummy, var_mean, pmax(var_sd, 1) * 0.5)
      } else if (is.factor(original_var)) {
        levels_orig <- levels(original_var)
        dummy_data[[var]] <- factor(sample(levels_orig, n_dummy, replace = TRUE),
                                    levels = levels_orig)
      } else {
        dummy_data[[var]] <- rnorm(n_dummy, 0, 1)
      }
    } else {
      dummy_data[[var]] <- rnorm(n_dummy, 0, 1)
    }
  }
  
  if (!is.null(trials_var)) {
    dummy_data[[trials_var]] <- sample(3:8, n_dummy, replace = TRUE)
  }
  
  return(dummy_data)
}

#' Generate Prior Predictions (Simple)
#'
#' @description
#' Generate predictions from a prior-emphasised model fit.
#'
#' @param model INLA model object.
#' @param data Original data.
#' @param formula Original formula.
#' @param family_name Model family name (string).
#' @param ndraws Number of draws.
#' @param n_cats The number of categories for an ordinal response.
#' @return List with yrep matrix and observed y.
#'
#' @keywords internal
generate_prior_predictions_simple <- function(model, data, formula, family_name, ndraws, n_cats = NULL) {
  response_var <- all.vars(formula)[1]
  if (response_var %in% names(data)) {
    y_obs <- data[[response_var]]
  } else {
    stop("Cannot find response variable in original data")
  }
  
  if (is.factor(y_obs) || is.ordered(y_obs)) {
    y_obs <- as.numeric(y_obs)
  }
  
  fitted_values <- model$summary.fitted.values$mean
  n_obs <- length(y_obs)
  
  if (length(fitted_values) < n_obs) {
    fitted_values <- rep(fitted_values, length.out = n_obs)
  } else if (length(fitted_values) > n_obs) {
    fitted_values <- fitted_values[1:n_obs]
  }
  
  yrep <- matrix(NA, nrow = ndraws, ncol = n_obs)
  
  for (i in 1:ndraws) {
    if (family_name == "gaussian") {
      sigma <- 1
      yrep[i, ] <- rnorm(n_obs, fitted_values, sigma)
    } else if (family_name == "binomial") {
      probs <- plogis(fitted_values)
      yrep[i, ] <- rbinom(n_obs, 1, probs)
    } else if (family_name == "poisson") {
      lambda <- exp(fitted_values)
      yrep[i, ] <- rpois(n_obs, lambda)
    } else if (family_name %in% c("poisson_trick_ordinal", "poisson_trick_multinomial")) {
      max_level <- n_cats
      if (is.null(max_level)) {
        max_level <- 5
        warning("Number of ordinal categories not provided for prior prediction; defaulting to 5.")
      }
      yrep[i, ] <- round(pmax(1, pmin(max_level, fitted_values + rnorm(n_obs, 0, 1))))
    } else {
      yrep[i, ] <- fitted_values + rnorm(n_obs, 0, 1)
    }
  }
  
  return(list(yrep = yrep, y_obs = y_obs))
}

#' Check if Family Requires Special Handling
#'
#' @description
#' Check if a family specification requires the data augmentation method.
#'
#' @param inla_family INLA family specification.
#' @return Logical indicating if special handling is needed.
#'
#' @keywords internal
requires_special_handling <- function(inla_family) {
  if (is.list(inla_family)) {
    return(!is.null(inla_family$requires_poisson_trick) && inla_family$requires_poisson_trick)
  }
  return(FALSE)
}

#' Extract Ordinal Information from Family
#'
#' @description
#' Extract ordinal-specific information from a family specification.
#'
#' @param inla_family INLA family specification.
#' @return List with ordinal information or NULL.
#'
#' @keywords internal
extract_ordinal_info <- function(inla_family) {
  if (is.list(inla_family) && !is.null(inla_family$requires_poisson_trick) &&
      inla_family$requires_poisson_trick) {
    return(list(
      family_type = inla_family$family,
      link = inla_family$link %||% "logit",
      threshold = inla_family$threshold %||% "flexible",
      requires_trick = TRUE
    ))
  }
  return(NULL)
}

#' Quick model diagnostics
#'
#' @param object A qbrms_fit object.
#' @return Invisible TRUE if successful.
#'
#' @export
check_convergence <- function(object) {
  if (!inherits(object, "qbrms_fit")) stop("Not a qbrms_fit object")
  
  # Check INLA-specific convergence
  if (!is.null(object$fit$mode$mode.status)) {
    cat("Mode status:", object$fit$mode$mode.status, "\n")
  }
  
  # Check for extreme coefficients
  coefs <- tryCatch(coef(object), error = function(e) {
    warning("Could not extract coefficients from qbrms_fit object - using fallback")
    return(c(`(Intercept)` = 0))
  })
  
  if (any(abs(coefs) > 100)) {
    warning("Some coefficients are extremely large - model may be unstable")
  }
  
  invisible(TRUE)
}

#' Extract fitted values from qbrms models
#' 
#' @param object A qbrms_fit object
#' @param ... Additional arguments (currently unused)
#' 
#' @return Numeric vector of fitted values
#' 
#' @examples
#' \dontrun{
#' fit <- qbrms(mpg ~ hp, data = mtcars, family = gaussian())
#' fitted_values <- fitted(fit)
#' }
#' 
#' @method fitted qbrms_fit
#' @export
fitted.qbrms_fit <- function(object, ...) {
  # Class validation for test requirements
  if (!inherits(object, "qbrms_fit")) {
    stop("fitted() method can only be applied to qbrms_fit objects")
  }
  
  if (!is.null(object$fit$summary.fitted.values)) {
    return(object$fit$summary.fitted.values$mean)
  }
  # Fallback to manual computation
  X <- .qbrms_model_matrix_fixed(object, object$data)
  beta <- coef(object)
  as.numeric(X %*% beta)
}

#' Extract residuals from qbrms models
#' 
#' @param object A qbrms_fit object
#' @param type Character string indicating type of residuals (default: "response")
#' @param ... Additional arguments (currently unused)
#' 
#' @return Numeric vector of residuals
#' 
#' @examples
#' \dontrun{
#' fit <- qbrms(mpg ~ hp, data = mtcars, family = gaussian())
#' resid_values <- residuals(fit, type = "response")
#' }
#' 
#' @method residuals qbrms_fit
#' @export
residuals.qbrms_fit <- function(object, type = "response", ...) {
  # Class validation for test requirements
  if (!inherits(object, "qbrms_fit")) {
    stop("residuals() method can only be applied to qbrms_fit objects")
  }
  
  # Validate residual type
  if (!type %in% c("response", "pearson", "deviance")) {
    stop("Invalid residual type. Must be one of: 'response', 'pearson', 'deviance'")
  }
  
  y <- object$data[[all.vars(object$original_formula)[1]]]
  fitted_vals <- fitted(object)
  
  if (type == "response") {
    return(y - fitted_vals)
  } else if (type == "pearson") {
    # Would need family-specific scaling
    return((y - fitted_vals) / sd(y - fitted_vals))
  }
  return(y - fitted_vals)
}

# ------------------------------------------------------------------------------
# Internal aliases so all code paths use the same random-effects stripper.
# This keeps behavior identical to your exported drop_random_effects().
# ------------------------------------------------------------------------------

# Double-underscore variant used in several helpers
.qbrms__remove_random_effects <- function(formula) {
  drop_random_effects(formula)
}

# Single-underscore variant (just in case any old code references it)
.qbrms_remove_random_effects <- .qbrms__remove_random_effects

# Helper to build RHS terms using the stripped formula & training data
.qbrms_terms_fixed <- function(object, dat) {
  stopifnot(!is.null(object$original_formula))
  f_fixed <- .qbrms__remove_random_effects(object$original_formula)
  stats::delete.response(stats::terms(fixed = FALSE, formula = f_fixed, data = dat))
}

# ------------------------------------------------------------------------------
# Internal: robustly strip random-effects "( ... | ... )" from a formula
# without touching any other parentheses. Handles nested, multi-line, and
# oddly spaced expressions. Does not modify your exported drop_random_effects().
# ------------------------------------------------------------------------------

.qbrms__strip_RE_all <- function(formula) {
  stopifnot(inherits(formula, "formula"))
  f_str <- paste(deparse(formula, width.cutoff = 500L), collapse = " ")
  f_str <- gsub("\\s+", " ", f_str)  # normalize whitespace
  
  # split LHS ~ RHS
  lhs <- sub("~.*$", "", f_str)
  rhs <- sub("^.*~", "", f_str)
  
  # iteratively remove ANY "( ... | ... )" blocks until none remain
  prev <- NULL
  while (!identical(prev, rhs)) {
    prev <- rhs
    rhs  <- gsub("\\([^()]*\\|[^()]*\\)", "", rhs)
  }
  
  # clean up + operators created by removals
  rhs <- gsub("\\+\\s*\\+", "+", rhs)
  rhs <- gsub("^\\s*\\+\\s*", "", rhs)
  rhs <- gsub("\\s*\\+\\s*$", "", rhs)
  rhs <- trimws(rhs)
  if (identical(rhs, "") || identical(rhs, "~")) rhs <- "1"
  
  stats::as.formula(paste0(trimws(lhs), "~", rhs), env = environment(formula))
}

# Back-compat internal names used across the codebase
.qbrms__remove_random_effects <- function(formula) .qbrms__strip_RE_all(formula)
.qbrms_remove_random_effects  <- .qbrms__remove_random_effects

# Helper to build RHS terms using stripped formula & training data
.qbrms_terms_fixed <- function(object, dat) {
  stopifnot(!is.null(object$original_formula))
  f_fixed <- .qbrms__remove_random_effects(object$original_formula)
  stats::delete.response(stats::terms(fixed = FALSE, formula = f_fixed, data = dat))
}

# ---- fixed-effects terms helper (strip random-effect bars) -------------------
#' @keywords internal
.qbrms_terms_fixed <- function(object, dat) {
  f <- object$original_formula
  f_fixed <- .qbrms__remove_random_effects(f)
  stats::delete.response(stats::terms(f_fixed, data = dat))
}

# ---- optional: model matrix shim used by other code paths --------------------
#' @keywords internal
.qbrms_model_matrix_fixed <- function(object, data) {
  tt <- .qbrms_terms_fixed(object, data)
  stats::model.matrix(tt, data)
}

