# =============================================================================
# R/model_fitting.R
# =============================================================================

# =============================================================================
# SECTION 1: INTERNAL HELPER FUNCTIONS (DEFINED FIRST)
# =============================================================================

#' Generate Prior Predictive Samples
#' @description Internal helper to generate samples from the prior predictive distribution.
#' @param formula Model formula.
#' @param data Data frame.
#' @param family Model family.
#' @param prior Prior specifications.
#' @param ndraws Number of draws.
#' @param verbose Logical; print progress.
#' @param ... Additional args.
#' @return Matrix of prior predictive samples.
#' @keywords internal
generate_prior_samples <- function(formula, data, family = gaussian(),
                                   prior = NULL, ndraws = 100,
                                   verbose = TRUE, ...) {
  n_dummy <- 50
  inla_family <- convert_family_to_inla(family)
  family_name <- extract_family_name(inla_family)
  
  dummy_data <- create_dummy_data_for_priors(formula, data, n_dummy, family_name, verbose = verbose)
  
  n_cats <- NULL
  if (family_name %in% c("poisson_trick_ordinal", "poisson_trick_multinomial")) {
    if (verbose) cat("Ordinal family detected, using Gaussian approximation for prior samples.\n")
    response_var <- all.vars(formula)[1]
    if (response_var %in% names(data) && is.factor(data[[response_var]])) {
      n_cats <- nlevels(data[[response_var]])
      if (verbose) cat("Detected", n_cats, "levels for the ordinal response.\n")
    }
    family_name <- "gaussian"
  }
  
  clean_formula <- formula
  trials_data <- NULL
  formula_str <- deparse(formula, width.cutoff = 500L)
  if (grepl("trials\\(", formula_str)) {
    term_obj <- terms(formula)
    vars <- attr(term_obj, "variables")
    response_var <- as.character(vars[[2]])
    predictor_labels <- attr(term_obj, "term.labels")
    predictor_labels <- predictor_labels[!grepl("^trials\\(", predictor_labels)]
    clean_formula <- reformulate(predictor_labels, response = response_var)
    
    trials_match <- regmatches(formula_str, regexpr("trials\\(([^)]+)\\)", formula_str))
    trials_var <- gsub("trials\\(|\\)", "", trials_match)
    trials_data <- dummy_data[[trials_var]]
  }
  
  control_compute <- list(dic = FALSE, waic = FALSE, cpo = FALSE, config = TRUE)
  control_family <- if (family_name == "gaussian") 
    list(hyper = list(prec = list(param = c(0.01, 0.01), prior = "loggamma")))
  else list()
  control_fixed <- list(mean.intercept = 0, prec.intercept = 0.16, mean = 0, prec = 0.16)
  
  model <- NULL
  tryCatch({
    if (family_name == "binomial" && !is.null(trials_data)) {
      model <- INLA::inla(clean_formula, data = dummy_data, family = "binomial",
                          control.compute = control_compute,
                          control.fixed = control_fixed,
                          control.family = control_family,
                          verbose = FALSE, ...)
    } else {
      model <- INLA::inla(clean_formula, data = dummy_data, family = family_name,
                          control.compute = control_compute,
                          control.fixed = control_fixed,
                          control.family = control_family,
                          verbose = FALSE, ...)
    }
  }, error = function(e) {
    if (verbose) cat("Prior sampling failed: falling back to Gaussian. Error:", e$message, "\n")
    model <<- INLA::inla(clean_formula, data = dummy_data, family = "gaussian",
                         control.compute = control_compute,
                         control.fixed = control_fixed,
                         control.family = control_family,
                         verbose = FALSE, ...)
    family_name <<- "gaussian"
  })
  
  prior_preds <- generate_prior_predictions_simple(model, data, formula, family_name, ndraws, n_cats = n_cats)
  return(prior_preds$yrep)
}

#' Convert Binomial Formulas Using 'trials()' Syntax
#' @keywords internal
convert_binomial_formula <- function(formula, data, formula_components, verbose = TRUE) {
  if (!formula_components$is_binomial_trials) {
    return(list(formula = formula, data = data, Ntrials = NULL))
  }
  
  if (verbose) cat("Processing binomial trials syntax.\n")
  
  formula_str <- deparse(formula)
  lhs <- strsplit(formula_str, "~")[[1]][1]
  rhs <- strsplit(formula_str, "~")[[1]][2]
  
  # Extract trial variable name
  trials_match <- regmatches(formula_str, regexpr("trials\\(([^)]+)\\)", formula_str))
  if (length(trials_match) == 0) stop("Cannot parse trials variable from formula")
  trials_var <- gsub("trials\\(|\\)", "", trials_match)
  
  success_var <- trimws(strsplit(lhs, "\\|")[[1]][1])
  
  # Create new formula with cbind(success, failures)
  new_formula <- as.formula(
    paste0("cbind(", success_var, ", ", trials_var, " - ", success_var, ") ~ ", rhs)
  )
  
  if (verbose) cat("Converted formula to:", deparse(new_formula), "\n")
  
  clean_data <- data
  if (!(success_var %in% names(clean_data))) stop("Success variable not in data")
  if (!(trials_var %in% names(clean_data))) stop("Trials variable not in data")
  
  clean_data[[trials_var]] <- pmax(as.integer(clean_data[[trials_var]]), 1L)
  clean_data[[success_var]] <- pmin(pmax(as.integer(clean_data[[success_var]]), 0L), clean_data[[trials_var]])
  
  # Remove incomplete rows
  complete_rows <- complete.cases(clean_data[[success_var]], clean_data[[trials_var]])
  clean_data <- clean_data[complete_rows, , drop = FALSE]
  
  list(formula = new_formula, data = clean_data)
}

#' Fit Mixed Effects Model Using INLA
#' @keywords internal
fit_mixed_effects_model_improved <- function(formula, data, inla_family, control.compute, verbose = TRUE, ...) {
  fit_start <- Sys.time()
  if (verbose) cat("Fitting mixed effects model...\n")
  
  data <- handle_missing_data(formula, data, verbose = verbose)
  family_name <- extract_family_name(inla_family)
  model <- NULL
  
  tryCatch({
    formula_components <- parse_formula_components(formula, data)
    
    if (family_name == "binomial" && formula_components$is_binomial_trials) {
      binom_setup <- convert_binomial_formula(formula, data, formula_components, verbose)
      clean_data <- binom_setup$data
      clean_formula <- binom_setup$formula
      
      if (verbose) {
        cat("Debug: Response range:", range(clean_data[[all.vars(clean_formula)[1]]], na.rm = TRUE), "\n")
      }
      
      model <- INLA::inla(
        formula = clean_formula,
        data = clean_data,
        family = family_name,
        control.compute = control.compute,
        verbose = FALSE,
        ...
      )
    } else {
      model <- INLA::inla(
        formula = formula,
        data = data,
        family = family_name,
        control.compute = control.compute,
        verbose = FALSE,
        ...
      )
    }
  }, error = function(e) {
    stop("INLA model fitting failed: ", e$message)
  })
  
  fit_end <- Sys.time()
  fitting_time <- as.numeric(difftime(fit_end, fit_start, units = "secs"))
  
  if (verbose) cat("Model fitting completed in", round(fitting_time, 2), "seconds.\n")
  
  list(fit = model, fitting_time = fitting_time)
}

#' Fit Fixed Effects Model Using INLA
#' @keywords internal
fit_fixed_effects_model_improved <- function(formula, data, inla_family, control.compute, verbose = TRUE, ...) {
  fit_start <- Sys.time()
  if (verbose) cat("Fitting fixed effects model...\n")
  
  data <- handle_missing_data(formula, data, verbose = verbose)
  family_name <- extract_family_name(inla_family)
  model <- NULL
  
  tryCatch({
    formula_components <- parse_formula_components(formula, data)
    
    if (family_name == "binomial" && formula_components$is_binomial_trials) {
      binom_setup <- convert_binomial_formula(formula, data, formula_components, verbose)
      clean_data <- binom_setup$data
      clean_formula <- binom_setup$formula
      
      if (verbose) {
        cat("Debug: Response range:", range(clean_data[[all.vars(clean_formula)[1]]], na.rm = TRUE), "\n")
      }
      
      model <- INLA::inla(
        formula = clean_formula,
        data = clean_data,
        family = family_name,
        control.compute = control.compute,
        verbose = FALSE,
        ...
      )
    } else {
      model <- INLA::inla(
        formula = formula,
        data = data,
        family = family_name,
        control.compute = control.compute,
        verbose = FALSE,
        ...
      )
    }
  }, error = function(e) {
    stop("INLA model fitting failed: ", e$message)
  })
  
  fit_end <- Sys.time()
  fitting_time <- as.numeric(difftime(fit_end, fit_start, units = "secs"))
  
  if (verbose) cat("Model fitting completed in", round(fitting_time, 2), "seconds.\n")
  
  list(fit = model, fitting_time = fitting_time)
}

#' Quantile Regression Using 'quantreg' Package
#' @keywords internal
quantile_regression_fit <- function(formula, data, quantile = 0.5, max_iter = 10, tol = 1e-6, verbose = TRUE) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required for quantile regression. Please install it.")
  }
  
  if (verbose) cat("Fitting quantile regression (tau = ", quantile, ")...\n", sep = "")
  
  fit_start <- Sys.time()
  qr_model <- quantreg::rq(formula, data = data, tau = quantile)
  qr_summary <- summary(qr_model, covariance = TRUE)
  
  coefs <- qr_summary$coefficients
  # Format output similar to INLA summary for compatibility
  summary_fixed <- data.frame(
    mean = coefs[, "Value"],
    sd = coefs[, "Std. Error"],
    `0.025quant` = NA,
    `0.5quant` = NA,
    `0.975quant` = NA,
    mode = NA,
    kld = NA,
    check.names = FALSE
  )
  rownames(summary_fixed) <- rownames(coefs)
  
  vcov_matrix <- qr_summary$cov
  
  fit_end <- Sys.time()
  fitting_time <- as.numeric(difftime(fit_end, fit_start, units = "secs"))
  
  if (verbose) cat("Quantile regression completed in", round(fitting_time, 2), "seconds.\n")
  
  list(
    summary.fixed = summary_fixed,
    vcov = vcov_matrix,
    cpu.used = list(total = fitting_time),
    call = qr_model$call
  )
}


#' Create A Quantile Regression Fitted Object for 'INLA'
#' @keywords internal
create_quantile_fit <- function(formula, data, quantile = 0.5, verbose = TRUE) {
  fit <- quantile_regression_fit(formula, data, quantile, verbose)
  class(fit) <- c("inla", "quantile_inla", "list")
  return(fit)
}


#' Main q-brms model fitting interface (simplified example)
#' @export
qbrms <- function(formula, data, family = gaussian(), prior = NULL, sample_prior = "no",
                  quantile = 0.5, control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                  verbose = TRUE, ...) {
  
  if (verbose) cat("Starting qbrms model fitting...\n")
  
  # Handle prior predictive sampling if requested
  if (sample_prior == "only") {
    if (verbose) cat("Sampling from prior predictive...\n")
    prior_samples <- generate_prior_samples(formula, data, family, prior, ndraws = 100, verbose = verbose, ...)
    result <- list(
      prior_samples = prior_samples,
      original_formula = formula,
      data = data,
      family = family,
      model_type = "prior_predictive"
    )
    class(result) <- c("qbrms_prior_fit", "qbrms_fit", "list")
    return(result)
  }
  
  # Parse formula and family for internal use
  parsed_formula <- parse_brms_formula(formula)
  inla_family <- convert_family_to_inla(family, quantile)
  family_name <- extract_family_name(inla_family)
  
  # If ordinal family uses augmented method, dispatch there
  if (family_name %in% c("poisson_trick_ordinal", "poisson_trick_multinomial")) {
    return(fit_ordinal_model_improved(parsed_formula$formula, data, inla_family, control.compute, verbose, ...))
  }
  
  # Quantile regression fit if requested family is asymmetric laplace
  if (family_name == "asymmetric_laplace") {
    fit_obj <- create_quantile_fit(parsed_formula$formula, data, quantile, verbose)
    result <- list(
      fit = fit_obj,
      original_formula = parsed_formula$formula,
      data = data,
      family = inla_family,
      model_type = "quantile_regression",
      quantile = quantile
    )
    class(result) <- c("qbrms_fit", "list")
    return(result)
  }
  
  # Determine if mixed or fixed effects
  formula_components <- parse_formula_components(parsed_formula$formula, data)
  if (formula_components$has_random_effects) {
    fit_result <- fit_mixed_effects_model_improved(parsed_formula$formula, data, inla_family, control.compute, verbose, ...)
  } else {
    fit_result <- fit_fixed_effects_model_improved(parsed_formula$formula, data, inla_family, control.compute, verbose, ...)
  }
  
  group_var <- NULL
  if (formula_components$has_random_effects) {
    formula_str <- deparse(parsed_formula$formula)
    group_pattern <- "\\(1 \\| ([^)]+)\\)"
    m <- regexpr(group_pattern, formula_str)
    if (m > 0) {
      group_var <- regmatches(formula_str, m)
      group_var <- gsub("\\(1 \\| |\\)", "", group_var)
      group_var <- trimws(group_var)
    }
  }
  
  result <- list(
    fit = fit_result$fit,
    original_formula = parsed_formula$formula,
    data = data,
    family = inla_family,
    model_type = ifelse(formula_components$has_random_effects, "mixed", "fixed"),
    fitting_time = fit_result$fitting_time,
    has_random_effects = formula_components$has_random_effects,
    group_var = group_var
  )
  class(result) <- c("qbrms_fit", "list")
  return(result)
}
#' Fit Ordinal Models Using Augmented Poisson Method in INLA
#' @keywords internal
fit_ordinal_model_improved <- function(formula, data, inla_family, control.compute, verbose = TRUE, ...) {
  if (verbose) cat("Fitting ordinal model via data augmentation...\n")
  fit_start <- Sys.time()
  
  data <- handle_missing_data(formula, data, verbose = verbose)
  
  response_var <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  formula_components <- parse_formula_components(formula, data)
  has_random <- formula_components$has_random_effects
  
  # Ensure response is ordered factor
  if (!is.ordered(data[[response_var]])) {
    if (is.factor(data[[response_var]])) {
      data[[response_var]] <- as.ordered(data[[response_var]])
      if (verbose) cat("Converted response factor to ordered factor.\n")
    } else {
      data[[response_var]] <- as.ordered(as.factor(data[[response_var]]))
      if (verbose) cat("Converted response to ordered factor.\n")
    }
  }
  
  y <- data[[response_var]]
  levels_y <- levels(y)
  n_cats <- nlevels(y)
  n_obs <- nrow(data)
  
  if (verbose) cat("Ordinal response:", n_cats, "categories:", paste(levels_y, collapse = " < "), "\n")
  if (n_cats < 3) stop("Need at least 3 ordinal levels for meaningful regression.")
  
  # Expand data for Poisson augmentation
  expanded_data <- data[rep(1:n_obs, each = n_cats), , drop = FALSE]
  category_indices <- rep(1:n_cats, times = n_obs)
  observed_categories <- rep(as.numeric(y), each = n_cats)
  
  expanded_data$y_poisson <- as.integer(category_indices == observed_categories)
  expanded_data$category_idx <- category_indices
  expanded_data$obs_idx <- rep(1:n_obs, each = n_cats)
  
  if (verbose) {
    cat("Expanded data from", n_obs, "to", nrow(expanded_data), "rows.\n")
    cat("Response sum check:", sum(expanded_data$y_poisson), "(should equal", n_obs, ")\n")
  }
  
  group_var <- NULL
  if (has_random) {
    formula_str <- deparse(formula, width.cutoff = 500)
    group_pattern <- "\\(\\s*1\\s*\\|\\s*([^\\)]+)\\)"
    group_matches <- regmatches(formula_str, gregexpr(group_pattern, formula_str, perl = TRUE))
    if (length(group_matches[[1]]) > 0) {
      group_var <- gsub("\\(\\s*1\\s*\\|\\s*|\\)", "", group_matches[[1]][1])
      group_var <- trimws(group_var)
      if (verbose) cat("Detected grouping variable:", group_var, "\n")
      if (is.factor(expanded_data[[group_var]])) {
        expanded_data$group_idx <- as.integer(expanded_data[[group_var]])
      } else {
        expanded_data$group_idx <- as.integer(as.factor(expanded_data[[group_var]]))
      }
    } else {
      has_random <- FALSE
      group_var <- NULL
    }
  }
  
  fixed_predictors <- if (has_random && !is.null(group_var)) {
    predictors[predictors != group_var]
  } else {
    predictors
  }
  
  predictor_part <- if (length(fixed_predictors) > 0) paste(fixed_predictors, collapse = " + ") else ""
  
  formula_parts <- c(
    "y_poisson ~ -1",
    "f(category_idx, model = 'iid', hyper = list(prec = list(initial = -10, fixed = TRUE)), constr = TRUE)"
  )
  
  if (predictor_part != "") {
    formula_parts <- c(formula_parts, predictor_part)
  }
  
  formula_parts <- c(formula_parts,
                     "f(obs_idx, model = 'iid', hyper = list(prec = list(prior = 'loggamma'), param = c(1, 0.001)))")
  
  if (has_random) {
    formula_parts <- c(formula_parts,
                       "f(group_idx, model = 'iid', hyper = list(prec = list(prior = 'loggamma'), param = c(1, 0.5)))")
  }
  
  inla_formula <- as.formula(paste(formula_parts, collapse = " + "))
  
  if (verbose) cat("INLA formula:", deparse(inla_formula), "\n")
  
  model <- NULL
  tryCatch({
    model <- INLA::inla(
      formula = inla_formula,
      data = expanded_data,
      family = "poisson",
      E = rep(1, nrow(expanded_data)),
      control.compute = control.compute,
      control.inla = list(strategy = "adaptive", int.strategy = "eb"),
      verbose = FALSE,
      ...
    )
    if (verbose) cat("Ordinal model fitting successful.\n")
  }, error = function(e) {
    if (verbose) cat("Ordinal model fitting error:", e$message, "\n")
    stop("Ordinal model fitting failed; try qrombs_ordinal_binary instead.")
  })
  
  fit_end <- Sys.time()
  fitting_time <- as.numeric(difftime(fit_end, fit_start, units = "secs"))
  
  result <- list(
    original_formula = formula,
    inla_formula = inla_formula,
    data = data,
    expanded_data = expanded_data,
    family = inla_family,
    fit = model,
    has_random = has_random,
    group_var = group_var,
    ordinal_levels = levels_y,
    n_categories = n_cats,
    model_type = "ordinal_augmented",
    fitting_time = fitting_time,
    timing = list(total_seconds = fitting_time,
                  formatted_duration = format_duration(fitting_time))
  )
  
  class(result) <- c("ordinal_augmented_qmbs_fit", "ordinal_qmbs_fit", "qmbs_fit")
  return(result)
}


#' Variance-Covariance Matrix for INLA Objects
#' @export
vcov.inla <- function(object, ...) {
  if (!is.null(object$misc$cov.fixed)) {
    return(object$misc$cov.fixed)
  } else {
    warning("Covariance matrix for fixed effects not available. Re-run with control.compute = list(config = TRUE) and use inla.vcov().")
    sds <- object$summary.fixed[, "sd"]
    vcov_matrix <- diag(sds^2)
    rownames(vcov_matrix) <- colnames(vcov_matrix) <- rownames(object$summary.fixed)
    return(vcov_matrix)
  }
}

#' Variance-Covariance Matrix for Quantile INLA Objects
#' @export
vcov.quantile_inla <- function(object, ...) {
  if (!is.null(object$vcov)) {
    return(object$vcov)
  } else {
    stop("Variance-covariance matrix not found in the quantile_inla object.")
  }
}
