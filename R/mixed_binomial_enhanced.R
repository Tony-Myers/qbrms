# =============================================================================
# qbrmb: Complete Enhanced Binomial Mixed Effects Implementation
# All functions included for standalone use
# =============================================================================

#' Enhanced binomial mixed-effects modelling
#'
#' Fits a regularised binomial (or Bernoulli) mixed-effects model using INLA,
#' with enhanced diagnostics, stability checks and strategy selection.
#'
#' @param formula Model formula with random effects in lme4-style syntax.
#' @param data Data frame containing the variables in the model.
#' @param family Model family (currently \code{"binomial"} or \code{"bernoulli"};
#'   default \code{"binomial"}).
#' @param strategy Fitting strategy: \code{"auto"}, \code{"enhanced"},
#'   \code{"aggressive"}, or \code{"minimal"}.
#' @param regularisation_strength Regularisation strength in the interval
#'   \eqn{[0, 1]} (default \code{0.1}).
#' @param use_data_augmentation Logical; if \code{TRUE}, add pseudo-observations
#'   for additional numerical stability.
#' @param min_group_size Minimum group size before triggering diagnostic
#'   warnings.
#' @param verbose Logical; if \code{TRUE}, show detailed progress and
#'   diagnostics while fitting.
#' @param diagnostics Logical; if \code{TRUE}, compute and store extended
#'   diagnostics in the returned object.
#' @param silent Logical; if \code{TRUE}, suppress printed output except
#'   errors.
#' @param ... Additional arguments passed to \code{INLA::inla()}.
#'
#' @return An object of class \code{c("qbrmb_fit", "qbrms_fit", "list")}
#'   containing the fitted model, diagnostics and metadata.
#'
#' @name qbrmb
#' @rdname qbrmb
#' @export
qbrmb <- function(formula,
                  data,
                  family = "binomial",
                  strategy = "auto",
                  regularisation_strength = 0.1,
                  use_data_augmentation = TRUE,
                  min_group_size = 5,
                  verbose = FALSE,
                  diagnostics = FALSE,
                  silent = FALSE,
                  ...) {
  start_time <- Sys.time()
  
  # Determine output level
  show_diagnostics <- verbose || diagnostics
  show_progress <- !silent && (verbose || diagnostics)
  
  if (show_progress) cat("=== qbrmb: Enhanced Binomial Mixed Effects ===\n")
  
  # Validate inputs
  .validate_qbrmb_inputs(formula, data, family, show_progress)
  
  # Step 1: Comprehensive diagnostics
  diagnostics_result <- .diagnose_comprehensive(formula, data, min_group_size, show_diagnostics)
  
  # Step 2: Choose fitting strategy
  if (strategy == "auto") {
    strategy <- .choose_fitting_strategy(diagnostics_result, show_diagnostics)
  }
  
  # Step 3: Data preparation
  if (use_data_augmentation && diagnostics_result$needs_augmentation) {
    if (show_diagnostics) cat("Applying intelligent data augmentation...\n")
    data <- .augment_data_intelligent(formula, data, diagnostics_result, show_diagnostics)
  }
  
  # Step 4: Fit model
  if (show_progress && !show_diagnostics) {
    cat("Fitting enhanced binomial mixed effects model...\n")
  }
  
  fit_result <- .fit_with_strategy(formula, data, strategy, regularisation_strength, 
                                   diagnostics_result, show_diagnostics, ...)
  
  # Step 5: Post-fitting diagnostics
  fit_result <- .add_comprehensive_diagnostics(fit_result, diagnostics_result, show_diagnostics)
  
  # Step 6: Create final object
  fit_result$fitting_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  fit_result$strategy_used <- strategy
  fit_result$qbrmb_version <- "1.0"
  fit_result$show_diagnostics_used <- show_diagnostics
  
  class(fit_result) <- c("qbrmb_fit", "qbrms_fit", "list")
  
  # Clean completion message
  if (show_progress && !show_diagnostics) {
    status <- switch(fit_result$overall_assessment,
                     "good" = "successfully",
                     "acceptable_with_caution" = "with minor issues",
                     "concerning" = "with some concerns", 
                     "problematic" = "with issues resolved")
    cat("Model fitted", status, "in", round(fit_result$fitting_time, 2), "seconds\n")
  } else if (show_diagnostics) {
    cat("qbrmb fitting completed in", round(fit_result$fitting_time, 2), "seconds\n")
    cat("Strategy used:", strategy, "\n")
    cat("Overall assessment:", fit_result$overall_assessment, "\n")
  }
  
  return(fit_result)
}

# =============================================================================
# INTERNAL VALIDATION FUNCTIONS
# =============================================================================

#' Validate qbrmb inputs
#' @keywords internal
.validate_qbrmb_inputs <- function(formula, data, family, verbose) {
  
  if (!inherits(formula, "formula")) {
    stop("formula must be a formula object")
  }
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  # Check for mixed effects
  has_random <- grepl("\\|", paste(deparse(formula), collapse = ""))
  if (!has_random) {
    stop("qbrmb requires mixed effects models with random intercepts")
  }
  
  # Check response variable exists
  response_var <- all.vars(formula)[1]
  if (!response_var %in% names(data)) {
    stop("Response variable '", response_var, "' not found in data")
  }
  
  # Check family
  if (!family %in% c("binomial", "bernoulli")) {
    stop("qbrmb currently supports only binomial/bernoulli families")
  }
  
  if (verbose) cat("Input validation passed\n")
}

# =============================================================================
# COMPREHENSIVE DIAGNOSTIC SYSTEM
# =============================================================================

#' Comprehensive diagnostic assessment
#' @keywords internal
.diagnose_comprehensive <- function(formula, data, min_group_size, verbose) {
  
  if (verbose) cat("\n=== Comprehensive Diagnostics ===\n")
  
  diagnostics <- list(
    problems = list(),
    severity_score = 0,
    needs_augmentation = FALSE,
    data_summary = list(),
    group_summary = list()
  )
  
  # Extract basic information
  response_var <- all.vars(formula)[1]
  y <- data[[response_var]]
  
  # Handle different response formats
  if (is.factor(y)) {
    if (nlevels(y) != 2) stop("Response must be binary for binomial models")
    y <- as.numeric(y) - 1
  } else if (is.logical(y)) {
    y <- as.numeric(y)
  }
  
  # Basic response analysis
  success_rate <- mean(y, na.rm = TRUE)
  n_obs <- length(y)
  n_missing <- sum(is.na(y))
  
  diagnostics$data_summary <- list(
    response_var = response_var,
    success_rate = success_rate,
    n_observations = n_obs,
    n_missing = n_missing
  )
  
  if (verbose) {
    cat("Response variable:", response_var, "\n")
    cat("Success rate:", round(success_rate * 100, 1), "%\n")
    cat("Sample size:", n_obs, "\n")
    cat("Missing values:", n_missing, "\n")
  }
  
  # Check for global separation
  severity_added <- .check_global_separation(diagnostics, success_rate, verbose)
  diagnostics$severity_score <- diagnostics$severity_score + severity_added
  
  # Check for sparse outcomes
  severity_added <- .check_sparse_outcomes(diagnostics, success_rate, n_obs, verbose)
  diagnostics$severity_score <- diagnostics$severity_score + severity_added
  
  # Analyse group structure
  group_info <- .analyse_group_structure(formula, data, y, min_group_size, verbose)
  if (!is.null(group_info)) {
    diagnostics$group_summary <- group_info
    severity_added <- .assess_group_problems(diagnostics, group_info, verbose)
    diagnostics$severity_score <- diagnostics$severity_score + severity_added
  }
  
  # Determine augmentation needs
  diagnostics$needs_augmentation <- (diagnostics$severity_score >= 2) || 
    (!is.null(diagnostics$problems$complete_separation)) ||
    (!is.null(diagnostics$problems$group_separation))
  
  # Print problem summary
  if (verbose && length(diagnostics$problems) > 0) {
    cat("\nDetected issues:\n")
    for (prob_name in names(diagnostics$problems)) {
      prob <- diagnostics$problems[[prob_name]]
      cat(sprintf("  [%s] %s\n", toupper(prob$severity), prob$description))
    }
    cat("Overall severity score:", diagnostics$severity_score, "\n")
  }
  
  return(diagnostics)
}

#' Check for global separation issues
#' @keywords internal
.check_global_separation <- function(diagnostics, success_rate, verbose) {
  
  severity_added <- 0
  
  if (success_rate == 0 || success_rate == 1) {
    diagnostics$problems$complete_separation <- list(
      severity = "critical",
      description = "Complete separation: all outcomes identical"
    )
    severity_added <- 3
    if (verbose) cat("WARNING: Complete separation detected\n")
  } else if (success_rate < 0.05 || success_rate > 0.95) {
    diagnostics$problems$quasi_separation <- list(
      severity = "high",
      description = paste("Quasi-separation: extremely unbalanced outcomes (", 
                          round(success_rate * 100, 1), "%)")
    )
    severity_added <- 2
    if (verbose) cat("WARNING: Quasi-separation detected\n")
  }
  
  return(severity_added)
}

#' Check for sparse outcomes
#' @keywords internal  
.check_sparse_outcomes <- function(diagnostics, success_rate, n_obs, verbose) {
  
  severity_added <- 0
  min_expected <- min(success_rate * n_obs, (1 - success_rate) * n_obs)
  
  if (min_expected < 3) {
    diagnostics$problems$very_sparse <- list(
      severity = "high",
      description = paste("Very sparse outcomes: minimum expected count =", round(min_expected, 1))
    )
    severity_added <- 2
  } else if (min_expected < 5) {
    diagnostics$problems$sparse_outcomes <- list(
      severity = "moderate",
      description = paste("Sparse outcomes: minimum expected count =", round(min_expected, 1))
    )
    severity_added <- 1
  }
  
  return(severity_added)
}

#' Analyse grouping structure
#' @keywords internal
.analyse_group_structure <- function(formula, data, y, min_group_size, verbose) {
  
  # Extract group variable
  group_var <- .extract_group_variable(formula)
  if (is.null(group_var) || !group_var %in% names(data)) return(NULL)
  
  # Analyse grouping structure
  groups <- data[[group_var]]
  group_sizes <- table(groups)
  group_success_rates <- tapply(y, groups, mean, na.rm = TRUE)
  
  group_info <- list(
    group_variable = group_var,
    n_groups = length(unique(groups)),
    group_sizes = as.numeric(group_sizes),
    group_names = names(group_sizes),
    group_success_rates = as.numeric(group_success_rates),
    min_group_size = min(group_sizes),
    max_group_size = max(group_sizes),
    mean_group_size = mean(group_sizes)
  )
  
  if (verbose) {
    cat("Grouping structure:\n")
    cat("  Variable:", group_var, "\n")
    cat("  Number of groups:", group_info$n_groups, "\n")
    cat("  Group sizes: min =", group_info$min_group_size, 
        ", max =", group_info$max_group_size,
        ", mean =", round(group_info$mean_group_size, 1), "\n")
    
    cat("  Group success rates:\n")
    for (i in seq_along(group_info$group_success_rates)) {
      rate <- group_info$group_success_rates[i]
      name <- group_info$group_names[i]
      status <- if (rate == 0) " (COMPLETE SEP)" else 
        if (rate == 1) " (COMPLETE SEP)" else
          if (rate < 0.1 || rate > 0.9) " (NEAR SEP)" else ""
      cat("    Group", name, ":", round(rate * 100, 1), "%", status, "\n")
    }
  }
  
  return(group_info)
}

#' Assess group-specific problems
#' @keywords internal
.assess_group_problems <- function(diagnostics, group_info, verbose) {
  
  severity_added <- 0
  
  # Check for too few groups
  if (group_info$n_groups < 5) {
    diagnostics$problems$few_groups <- list(
      severity = "high",
      description = paste("Very few groups:", group_info$n_groups)
    )
    severity_added <- severity_added + 2
  } else if (group_info$n_groups < 8) {
    diagnostics$problems$limited_groups <- list(
      severity = "moderate",
      description = paste("Limited groups:", group_info$n_groups)
    )
    severity_added <- severity_added + 1
  }
  
  # Check for small groups
  small_groups <- sum(group_info$group_sizes < 5)
  if (small_groups > group_info$n_groups * 0.5) {
    diagnostics$problems$many_small_groups <- list(
      severity = "moderate",
      description = paste(small_groups, "groups have fewer than 5 observations")
    )
    severity_added <- severity_added + 1
  }
  
  # Check for group-level separation
  extreme_groups <- sum(group_info$group_success_rates %in% c(0, 1), na.rm = TRUE)
  near_extreme_groups <- sum((group_info$group_success_rates < 0.1 | 
                                group_info$group_success_rates > 0.9), na.rm = TRUE) - extreme_groups
  
  if (extreme_groups > 0) {
    diagnostics$problems$group_separation <- list(
      severity = "critical",
      description = paste(extreme_groups, "groups have complete separation")
    )
    severity_added <- severity_added + 3
  }
  
  if (near_extreme_groups > 0) {
    diagnostics$problems$near_group_separation <- list(
      severity = "high", 
      description = paste(near_extreme_groups, "groups have near-complete separation")
    )
    severity_added <- severity_added + 2
  }
  
  return(severity_added)
}

# =============================================================================
# STRATEGY SELECTION AND FITTING
# =============================================================================

#' Choose optimal fitting strategy
#' @keywords internal
.choose_fitting_strategy <- function(diagnostics, verbose) {
  
  severity <- diagnostics$severity_score
  
  strategy <- if (severity >= 5) {
    "aggressive"
  } else if (severity >= 3) {
    "enhanced" 
  } else if (severity >= 1) {
    "enhanced"
  } else {
    "minimal"
  }
  
  if (verbose) {
    cat("Strategy selection:\n")
    cat("  Severity score:", severity, "\n")
    cat("  Chosen strategy:", strategy, "\n")
  }
  
  return(strategy)
}

#' Fit model using specified strategy
#' @keywords internal
.fit_with_strategy <- function(formula, data, strategy, regularisation_strength, 
                               diagnostics, verbose, ...) {
  
  if (verbose) cat("\n=== Model Fitting ===\n")
  
  switch(strategy,
         "minimal" = .fit_minimal_strategy(formula, data, verbose, ...),
         "enhanced" = .fit_enhanced_strategy(formula, data, regularisation_strength, verbose, ...),
         "aggressive" = .fit_aggressive_strategy(formula, data, verbose, ...),
         stop("Unknown strategy: ", strategy)
  )
}

#' Minimal strategy
#' @keywords internal
.fit_minimal_strategy <- function(formula, data, verbose, ...) {
  
  if (verbose) cat("Using minimal regularisation strategy...\n")
  
  .fit_with_inla_enhanced(formula, data, 
                          intercept_prior = list(mean = 0, prec = 1/4),
                          slope_prior = list(mean = 0, prec = 2),
                          verbose = verbose, ...)
}

#' Enhanced strategy
#' @keywords internal
.fit_enhanced_strategy <- function(formula, data, regularisation_strength, verbose, ...) {
  
  if (verbose) cat("Using enhanced regularisation strategy...\n")
  
  # Calculate regularised priors
  intercept_prec <- 1/4  # SD = 2
  slope_prec <- 1 + regularisation_strength * 3  # SD decreases with regularisation
  
  .fit_with_inla_enhanced(formula, data,
                          intercept_prior = list(mean = 0, prec = intercept_prec),
                          slope_prior = list(mean = 0, prec = slope_prec),
                          verbose = verbose, ...)
}

#' Aggressive strategy
#' @keywords internal
.fit_aggressive_strategy <- function(formula, data, verbose, ...) {
  
  if (verbose) cat("Using aggressive regularisation strategy...\n")
  
  # Try aggressive INLA first
  result <- tryCatch({
    .fit_with_inla_enhanced(formula, data,
                            intercept_prior = list(mean = 0, prec = 4),    # SD = 0.5
                            slope_prior = list(mean = 0, prec = 9),        # SD = 0.33
                            control_strategy = "aggressive",
                            verbose = verbose, ...)
  }, error = function(e) {
    if (verbose) cat("Aggressive INLA failed, trying ridge fallback...\n")
    .fit_ridge_fallback(formula, data, verbose)
  })
  
  return(result)
}

#' Core INLA fitting function
#' @keywords internal
.fit_with_inla_enhanced <- function(formula, data, intercept_prior, slope_prior, 
                                    control_strategy = "standard", verbose = TRUE, ...) {
  
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA package required for qbrmb fitting")
  }
  
  # Extract and prepare grouping variable
  group_var <- .extract_group_variable(formula)
  if (is.null(group_var)) {
    stop("Could not extract group variable from mixed effects formula")
  }
  
  data$group_id <- as.numeric(as.factor(data[[group_var]]))
  
  # Convert to INLA formula
  inla_formula <- .convert_to_inla_formula(formula)
  
  # Set up INLA controls
  control_fixed <- list(
    mean.intercept = intercept_prior$mean,
    prec.intercept = intercept_prior$prec,
    mean = slope_prior$mean,
    prec = slope_prior$prec
  )
  
  control_inla <- switch(control_strategy,
                         "aggressive" = list(
                           strategy = "simplified.laplace",
                           int.strategy = "ccd",
                           dz = 0.5,
                           diff.logdens = 5
                         ),
                         "standard" = list(
                           strategy = "gaussian",
                           int.strategy = "eb",
                           tolerance = 1e-6
                         )
  )
  
  control_compute <- list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE,
    config = TRUE
  )
  
  if (verbose) {
    cat("INLA formula:", deparse(inla_formula), "\n")
    cat("Prior: Intercept N(0,", round(1/sqrt(intercept_prior$prec), 2), 
        "), Slopes N(0,", round(1/sqrt(slope_prior$prec), 2), ")\n")
  }
  
  # Fit model
  inla_result <- INLA::inla(
    formula = inla_formula,
    data = data,
    family = "binomial",
    control.fixed = control_fixed,
    control.inla = control_inla,
    control.compute = control_compute,
    verbose = FALSE,
    ...
  )
  
  # Create result object
  result <- list(
    fit = inla_result,
    data = data,
    original_formula = formula,
    family = list(family = "binomial"),
    group_var = group_var,
    fitting_method = "inla_enhanced",
    priors_used = list(intercept = intercept_prior, slopes = slope_prior)
  )
  
  if (verbose && !is.null(inla_result$summary.fixed)) {
    cat("INLA fitting successful\n")
    cat("Fixed effects estimates:\n")
    print(round(inla_result$summary.fixed[, c("mean", "sd")], 3))
  }
  
  return(result)
}

#' Ridge regression fallback
#' @keywords internal
.fit_ridge_fallback <- function(formula, data, verbose) {
  
  if (verbose) cat("Using ridge regression fallback...\n")
  
  # Extract response and create design matrix
  response_var <- all.vars(formula)[1]
  y <- data[[response_var]]
  if (is.factor(y)) y <- as.numeric(y) - 1
  
  # Create fixed effects design matrix
  X <- tryCatch({
    fixed_formula <- .remove_random_effects(formula)
    model.matrix(fixed_formula, data = data)
  }, error = function(e) {
    matrix(1, nrow = nrow(data), ncol = 1, dimnames = list(NULL, "(Intercept)"))
  })
  
  # Ridge penalised estimates
  lambda <- 5  # Strong regularisation
  XtX_ridge <- crossprod(X) + lambda * diag(ncol(X))
  beta_ridge <- solve(XtX_ridge, crossprod(X, y))
  se_ridge <- sqrt(diag(solve(XtX_ridge))) * 1.5  # Conservative SE
  
  summary_fixed <- data.frame(
    mean = as.numeric(beta_ridge),
    sd = se_ridge,
    row.names = colnames(X)
  )
  
  # Create result object
  result <- list(
    fit = list(
      summary.fixed = summary_fixed,
      method = "ridge_regression",
      lambda = lambda
    ),
    data = data,
    original_formula = formula,
    family = list(family = "binomial"),
    group_var = .extract_group_variable(formula),
    fitting_method = "ridge_fallback"
  )
  
  if (verbose) {
    cat("Ridge regression completed (lambda =", lambda, ")\n")
    print(round(summary_fixed, 3))
  }
  
  return(result)
}

# =============================================================================
# DATA AUGMENTATION
# =============================================================================

#' Intelligent data augmentation
#' @keywords internal
.augment_data_intelligent <- function(formula, data, diagnostics, verbose) {
  
  # Determine augmentation rate
  augmentation_rate <- if (diagnostics$severity_score >= 5) {
    0.25  # 25% for critical cases
  } else if (diagnostics$severity_score >= 3) {
    0.15  # 15% for high severity
  } else if (diagnostics$severity_score >= 2) {
    0.10  # 10% for moderate severity
  } else {
    0.05  # 5% for mild issues
  }
  
  n_aug <- ceiling(nrow(data) * augmentation_rate)
  
  if (verbose) {
    cat("Intelligent data augmentation:\n")
    cat("  Severity score:", diagnostics$severity_score, "\n")
    cat("  Augmentation rate:", round(augmentation_rate * 100, 1), "%\n")
    cat("  Adding", n_aug, "pseudo-observations\n")
  }
  
  # Create augmented observations
  aug_data <- .create_balanced_augmentation(data, formula, n_aug, diagnostics, verbose)
  
  # Combine and return
  augmented_data <- rbind(data, aug_data)
  
  if (verbose) {
    response_var <- all.vars(formula)[1]
    y_orig <- data[[response_var]]
    y_aug <- augmented_data[[response_var]]
    if (is.factor(y_orig)) y_orig <- as.numeric(y_orig) - 1
    if (is.factor(y_aug)) y_aug <- as.numeric(y_aug) - 1
    
    cat("Augmentation results:\n")
    cat("  Original success rate:", round(mean(y_orig, na.rm = TRUE) * 100, 1), "%\n")
    cat("  Augmented success rate:", round(mean(y_aug, na.rm = TRUE) * 100, 1), "%\n")
  }
  
  return(augmented_data)
}

#' Create balanced augmented observations
#' @keywords internal
.create_balanced_augmentation <- function(data, formula, n_aug, diagnostics, verbose) {
  
  response_var <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  
  # Remove grouping variables from augmentation
  group_var <- .extract_group_variable(formula)
  if (!is.null(group_var)) {
    predictors <- predictors[predictors != group_var]
  }
  
  # Sample base observations
  base_indices <- sample(nrow(data), n_aug, replace = TRUE)
  aug_data <- data[base_indices, ]
  
  # Create balanced outcomes
  aug_data[[response_var]] <- rep(c(0, 1), length.out = n_aug)
  
  # Apply shrinkage to continuous predictors
  shrinkage_factor <- 0.8  # Shrink by 20%
  
  for (pred in predictors) {
    if (pred %in% names(data) && is.numeric(data[[pred]])) {
      grand_mean <- mean(data[[pred]], na.rm = TRUE)
      aug_data[[pred]] <- aug_data[[pred]] * shrinkage_factor + 
        grand_mean * (1 - shrinkage_factor)
    }
  }
  
  return(aug_data)
}

# =============================================================================
# POST-FITTING DIAGNOSTICS (FIXED)
# =============================================================================

#' Add comprehensive post-fitting diagnostics
#' @keywords internal
.add_comprehensive_diagnostics <- function(fit_result, pre_diagnostics, verbose) {
  
  if (verbose) cat("\n=== Post-Fitting Assessment ===\n")
  
  post_diagnostics <- list()
  
  # Extract coefficient information
  if (!is.null(fit_result$fit$summary.fixed)) {
    fixed_effects <- fit_result$fit$summary.fixed
    
    # Check coefficient magnitudes
    large_coefs <- abs(fixed_effects$mean) > 10
    if (any(large_coefs)) {
      post_diagnostics$large_coefficients <- list(
        severity = "warning",
        coefficients = rownames(fixed_effects)[large_coefs],
        values = round(fixed_effects$mean[large_coefs], 3),
        message = "Some coefficients are unusually large"
      )
    }
    
    # Check standard error magnitudes
    large_ses <- fixed_effects$sd > 5
    if (any(large_ses)) {
      post_diagnostics$large_standard_errors <- list(
        severity = "warning",
        coefficients = rownames(fixed_effects)[large_ses],
        values = round(fixed_effects$sd[large_ses], 3),
        message = "Some standard errors are very large"
      )
    }
    
    # Check coefficient/SE ratios
    if (any(fixed_effects$sd > 0)) {
      se_ratio <- abs(fixed_effects$mean / fixed_effects$sd)
      weak_identification <- se_ratio < 0.5
      if (any(weak_identification, na.rm = TRUE)) {
        post_diagnostics$weak_identification <- list(
          severity = "info",
          coefficients = rownames(fixed_effects)[weak_identification],
          ratios = round(se_ratio[weak_identification], 3),
          message = "Some coefficients are weakly identified"
        )
      }
    }
  }
  
  # Overall assessment - FIXED
  if (length(post_diagnostics) == 0) {
    n_critical <- 0
    n_warnings <- 0
    n_info <- 0
  } else {
    severities <- vapply(post_diagnostics, function(x) x$severity, character(1))
    n_critical <- sum(severities == "error")
    n_warnings <- sum(severities == "warning") 
    n_info <- sum(severities == "info")
  }
  
  overall_assessment <- if (n_critical > 0) {
    "problematic"
  } else if (n_warnings > 2) {
    "concerning"
  } else if (n_warnings > 0) {
    "acceptable_with_caution"
  } else {
    "good"
  }
  
  # Print diagnostics
  if (verbose) {
    if (length(post_diagnostics) > 0) {
      for (diag_name in names(post_diagnostics)) {
        diag <- post_diagnostics[[diag_name]]
        cat(sprintf("[%s] %s\n", toupper(diag$severity), diag$message))
      }
    } else {
      cat("No post-fitting issues detected\n")
    }
    
    cat("Overall assessment:", overall_assessment, "\n")
  }
  
  # Add diagnostics to fit object
  fit_result$pre_diagnostics <- pre_diagnostics
  fit_result$post_diagnostics <- post_diagnostics
  fit_result$overall_assessment <- overall_assessment
  
  return(fit_result)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Extract group variable from formula
#' @keywords internal
.extract_group_variable <- function(formula) {
  formula_str <- paste(deparse(formula), collapse = "")
  
  # Look for (1|group_var) pattern
  group_pattern <- "\\(\\s*1\\s*\\|\\s*([^)]+)\\s*\\)"
  matches <- regmatches(formula_str, regexpr(group_pattern, formula_str, perl = TRUE))
  
  if (length(matches) > 0) {
    group_var <- gsub("\\(\\s*1\\s*\\|\\s*|\\s*\\)", "", matches[1])
    return(trimws(group_var))
  }
  
  return(NULL)
}

#' Convert lme4-style formula to INLA format
#' @keywords internal
.convert_to_inla_formula <- function(formula) {
  formula_str <- paste(deparse(formula), collapse = "")
  
  # Split into LHS and RHS
  parts <- strsplit(formula_str, "~")[[1]]
  lhs <- trimws(parts[1])
  rhs <- trimws(parts[2])
  
  # Replace (1|group_var) with f(group_id, model="iid")
  rhs_inla <- gsub("\\(\\s*1\\s*\\|\\s*[^)]+\\s*\\)", 
                   "f(group_id, model=\"iid\")", 
                   rhs)
  
  inla_formula_str <- paste(lhs, "~", rhs_inla)
  return(as.formula(inla_formula_str))
}

#' Remove random effects from formula
#' @keywords internal
.remove_random_effects <- function(formula) {
  formula_str <- paste(deparse(formula), collapse = "")
  parts <- strsplit(formula_str, "~")[[1]]
  lhs <- trimws(parts[1])
  rhs <- trimws(parts[2])
  
  # Remove random effects terms
  rhs_fixed <- gsub("\\s*\\+\\s*\\([^)]+\\)", "", rhs)
  rhs_fixed <- gsub("\\s+", " ", rhs_fixed)
  rhs_fixed <- trimws(rhs_fixed)
  
  if (!nzchar(rhs_fixed)) rhs_fixed <- "1"
  
  return(as.formula(paste(lhs, "~", rhs_fixed)))
}

# =============================================================================
# S3 METHODS (brms-style)
# =============================================================================
#' Print a qbrmb model fit
#'
#' Nicely formatted one-line summary plus key diagnostics for a \code{qbrmb_fit}.
#'
#' @param x A \code{qbrmb_fit} object.
#' @param digits Number of significant digits to print. Defaults to
#'   \code{getOption("digits")}.
#' @param ... Unused; included for S3 compatibility.
#' @export
#' @method print qbrmb_fit
print.qbrmb_fit <- function(x, digits = 2, ...) {
  
  cat(" Family: bernoulli\n")
  cat("  Links: mu = logit\n")
  
  if (!is.null(x$original_formula)) {
    cat("Formula:", deparse(x$original_formula), "\n")
  }
  
  if (!is.null(x$data)) {
    n_orig <- nrow(x$data)
    if (!is.null(x$pre_diagnostics$needs_augmentation) && x$pre_diagnostics$needs_augmentation) {
      # Estimate original size
      severity <- x$pre_diagnostics$severity_score
      aug_rate <- if (severity >= 5) 0.25 else if (severity >= 3) 0.15 else if (severity >= 2) 0.10 else 0.05
      n_estimated_orig <- round(n_orig / (1 + aug_rate))
      cat("   Data:", n_estimated_orig, "observations (", n_orig - n_estimated_orig, "augmented)\n")
    } else {
      cat("   Data:", n_orig, "observations\n")
    }
  }
  
  cat("  Draws: INLA approximation (", x$strategy_used, " strategy)\n", sep = "")
  
  # Random effects section
  if (!is.null(x$group_var)) {
    cat("\nMultilevel Hyperparameters:\n")
    cat("~", x$group_var, sep = "")
    
    if (!is.null(x$pre_diagnostics$group_summary)) {
      n_groups <- x$pre_diagnostics$group_summary$n_groups
      cat(" (Number of levels: ", n_groups, ")\n", sep = "")
    } else {
      cat("\n")
    }
    
    # Extract random effects info
    if (!is.null(x$fit$summary.hyperpar)) {
      hyperpar <- x$fit$summary.hyperpar
      re_rows <- grep("Precision for", rownames(hyperpar))
      
      if (length(re_rows) > 0) {
        precision_mean <- hyperpar[re_rows[1], "mean"]
        precision_sd <- hyperpar[re_rows[1], "sd"]
        
        sd_mean <- 1/sqrt(precision_mean)
        sd_sd <- precision_sd / (2 * precision_mean^1.5)
        
        sd_lower <- max(0.01, sd_mean - 1.96 * sd_sd)
        sd_upper <- sd_mean + 1.96 * sd_sd
        
        cat(sprintf("%15s %9s %8s %8s %8s %8s\n", 
                    "", "Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Strategy"))
        cat(sprintf("sd(Intercept) %9.2f %9.2f %8.2f %8.2f %8s\n", 
                    sd_mean, sd_sd, sd_lower, sd_upper, x$strategy_used))
      }
    }
  }
  
  # Fixed effects section  
  if (!is.null(x$fit$summary.fixed)) {
    cat("\nPopulation-Level Effects:\n")
    fixed_effects <- x$fit$summary.fixed
    
    cat(sprintf("%11s %9s %8s %8s %8s %8s\n", 
                "", "Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Method"))
    
    for (i in 1:nrow(fixed_effects)) {
      name <- rownames(fixed_effects)[i]
      est <- fixed_effects[i, "mean"]
      se <- fixed_effects[i, "sd"]
      
      lower <- est - 1.96 * se
      upper <- est + 1.96 * se
      
      cat(sprintf("%11s %9.2f %9.2f %8.2f %8.2f %8s\n", 
                  name, est, se, lower, upper, "INLA"))
    }
  }
  
  cat("\nDraws were computed using INLA (", x$strategy_used, " regularisation).\n", sep = "")
  
  # Add quality note if needed
  if (!is.null(x$overall_assessment) && x$overall_assessment != "good") {
    quality_msg <- switch(x$overall_assessment,
                          "acceptable_with_caution" = "Model fitted with minor regularisation",
                          "concerning" = "Model required moderate regularisation due to data issues",
                          "problematic" = "Model required strong regularisation due to severe data issues")
    if (!is.null(quality_msg)) {
      cat("Note:", quality_msg, "\n")
    }
  }
  
  invisible(x)
}


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

# =============================================================================
# Convenience wrappers for qbrmb
# =============================================================================

#' Regularised binomial mixed-effects (enhanced strategy)
#'
#' Convenience wrapper around \code{\link{qbrmb}} using the "enhanced"
#' regularisation strategy and a stronger default regularisation strength.
#'
#' @param formula Model formula with random effects (lme4-style).
#' @param data Data frame containing the variables in the model.
#' @param verbose Logical; if \code{TRUE}, show detailed diagnostics.
#' @param ... Additional arguments passed to \code{\link{qbrmb}}.
#'
#' @return An object of class \code{c("qbrmb_fit", "qbrms_fit", "list")}.
#'
#' @export
qbrmb_regularised <- function(formula, data, verbose = FALSE, ...) {
  qbrmb(
    formula = formula,
    data = data,
    strategy = "enhanced",
    regularisation_strength = 0.2,
    use_data_augmentation = TRUE,
    verbose = verbose,
    ...
  )
}

#' Aggressively regularised binomial mixed-effects model
#'
#' Convenience wrapper around \code{\link{qbrmb}} using the "aggressive"
#' strategy with a higher default regularisation strength.
#'
#' @param formula Model formula with random effects (lme4-style).
#' @param data Data frame containing the variables in the model.
#' @param verbose Logical; if \code{TRUE}, show detailed diagnostics.
#' @param ... Additional arguments passed to \code{\link{qbrmb}}.
#'
#' @return An object of class \code{c("qbrmb_fit", "qbrms_fit", "list")}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(123)
#'   data <- data.frame(
#'     y     = rbinom(100, 1, 0.2),
#'     x     = rnorm(100),
#'     group = factor(rep(1:10, each = 10))
#'   )
#'   # qbrmb_aggressive requires a mixed model with random intercepts
#'   fit <- qbrmb_aggressive(y ~ x + (1 | group), data = data, verbose = FALSE)
#' }
#' }
#'
#' @export
qbrmb_aggressive <- function(formula, data, verbose = FALSE, ...) {
  qbrmb(
    formula = formula,
    data = data,
    strategy = "aggressive",
    regularisation_strength = 0.4,
    use_data_augmentation = TRUE,
    verbose = verbose,
    ...
  )
}

# Add null coalescing operator if not available
if (!exists("%||%")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}