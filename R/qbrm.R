# =============================================================================
# R/qbrm.R
# Main user interface function for qbrms with diagnostic capabilities
# =============================================================================

#' Quick Bayesian Regression Models
#'
#' @description
#' Enhanced interface to qbrms with all required parameters and built-in diagnostics.
#'
#' @param formula Model formula in lme4/brms style.
#' @param data Data frame containing the variables in the model.
#' @param family Model family (default: gaussian()).
#' @param prior Prior specifications (default: NULL).
#' @param sample_prior Whether to sample from priors ("no", "yes", "only"). Default: "no".
#' @param verbose Logical; print diagnostic information (default: TRUE).
#' @param ... Additional arguments passed to qbrms().
#'
#' @return A qbrms_fit object with model results.
#' @export
qbrm <- function(formula, data, family = gaussian(), prior = NULL, 
                 sample_prior = "no", verbose = TRUE, ...) {
  
  # Input validation
  if (missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' arguments are required")
  }
  
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  # Check for potential issues with binomial mixed effects models
  # Note: We use tryCatch to safely check family name in case of complex objects
  fam_name <- tryCatch(extract_family_name(family), error = function(e) "unknown")
  
  if (fam_name == "binomial") {
    diagnosis <- .diagnose_binomial_issues(formula, data, verbose)
    
    if (diagnosis$has_issues) {
      if (verbose) {
        cat("Potential identification issues detected in mixed effects binomial model:\n")
        for (issue in diagnosis$issues) {
          cat("  -", issue, "\n")
        }
        cat("Consider using qbrms_binomial_regularised() for more stable results.\n")
      }
    }
  }
  
  if (verbose) {
    cat("Starting qbrms model fitting...\n")
  }
  
  # Call qbrms directly (it handles its own safety internally)
  result <- qbrms(
    formula = formula, 
    data = data, 
    family = family, 
    prior = prior,
    sample_prior = sample_prior,
    verbose = verbose,
    ...
  )
  
  if (verbose) {
    cat("qbrms model fitting completed.\n")
  }
  
  return(result)
}

#' Internal function to diagnose binomial issues
#' @keywords internal
.diagnose_binomial_issues <- function(formula, data, verbose = TRUE) {
  
  issues <- character(0)
  has_issues <- FALSE
  
  # Extract response variable
  response_var <- all.vars(formula)[1]
  if (!response_var %in% names(data)) {
    return(list(has_issues = FALSE, issues = character(0)))
  }
  
  y <- data[[response_var]]
  
  # Check success rate for binary outcomes
  if (is.numeric(y) && all(y %in% c(0, 1))) {
    success_rate <- mean(y, na.rm = TRUE)
    # Flag if very unbalanced
    if (success_rate < 0.1 || success_rate > 0.9) {
      has_issues <- TRUE
      issues <- c(issues, "Highly unbalanced outcomes detected")
    }
    
    # Check for complete separation
    if (success_rate == 0 || success_rate == 1) {
      has_issues <- TRUE
      issues <- c(issues, "Complete separation detected")
    }
  }
  
  # Check for random effects issues
  has_random <- grepl("\\|", paste(deparse(formula), collapse = ""))
  if (has_random) {
    # Extract grouping variable
    formula_str <- paste(deparse(formula), collapse = "")
    random_match <- regmatches(formula_str, regexpr("\\([^|]*\\|[^)]*\\)", formula_str))
    
    if (length(random_match) > 0) {
      group_var <- trimws(gsub(".*\\|\\s*([^)]+)\\).*", "\\1", random_match))
      
      if (group_var %in% names(data)) {
        # Check group sizes
        group_sizes <- table(data[[group_var]])
        small_groups <- sum(group_sizes < 5)
        
        if (small_groups > 0) {
          has_issues <- TRUE
          issues <- c(issues, paste(small_groups, "groups with < 5 observations"))
        }
        
        # Check for group-level separation in binary outcomes
        if (is.numeric(y) && all(y %in% c(0, 1))) {
          group_rates <- tapply(y, data[[group_var]], mean, na.rm = TRUE)
          separated_groups <- sum(group_rates %in% c(0, 1), na.rm = TRUE)
          
          if (separated_groups > 0) {
            has_issues <- TRUE
            issues <- c(issues, paste(separated_groups, "groups with complete separation"))
          }
        }
      }
    }
  }
  
  return(list(has_issues = has_issues, issues = issues))
}

#' Fixed Regularised Binomial Mixed Effects Fitting
#'
#' @description
#' Fits binomial mixed effects models with regularisation, with all parameters handled correctly.
#'
#' @param formula Model formula with mixed effects structure.
#' @param data Data frame containing the variables.
#' @param regularise Logical; if TRUE, apply regularisation techniques.
#' @param sample_prior Whether to sample from priors ("no", "yes", "only"). Default: "no".
#' @param verbose Logical; print progress information.
#' @param ... Additional arguments passed to qbrms().
#'
#' @return A qbrms_fit object with additional regularisation metadata.
#' @export
qbrms_binomial_regularised <- function(formula, data, regularise = TRUE, 
                                       sample_prior = "no", verbose = TRUE, ...) {
  
  if (!regularise) {
    if (verbose) cat("Using standard qbrms fitting...\n")
    return(qbrms(formula = formula, data = data, family = binomial(), 
                 sample_prior = sample_prior, verbose = verbose, ...))
  }
  
  if (verbose) cat("Using regularised binomial fitting...\n")
  
  # Apply data augmentation for separation issues
  augmented_data <- .augment_data_for_stability(data, formula, verbose)
  
  # Define regularising priors
  reg_priors <- c(
    prior(normal(0, 1.5), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd")
  )
  
  if (verbose) {
    cat("Applying regularisation techniques:\n")
    cat("  - Data augmentation:", nrow(augmented_data) - nrow(data), "observations added\n")
    cat("  - Regularising priors applied\n")
  }
  
  # Fit the model with error handling
  tryCatch({
    if (verbose) cat("Starting regularised qbrms fitting...\n")
    
    fit <- qbrms(
      formula = formula,
      data = augmented_data,
      family = binomial(),
      prior = reg_priors,
      sample_prior = sample_prior,
      verbose = verbose,
      ...
    )
    
    if (verbose) cat("Regularised qbrms fitting completed.\n")
    
    # Add metadata about regularisation
    fit$regularisation_info <- list(
      applied = TRUE,
      original_n = nrow(data),
      augmented_n = nrow(augmented_data),
      augmentation_count = nrow(augmented_data) - nrow(data)
    )
    
    return(fit)
    
  }, error = function(e) {
    if (verbose) {
      cat("Regularised fitting failed, trying standard approach:\n")
      cat("Error was:", e$message, "\n")
    }
    
    # Fall back to standard approach
    return(qbrms(formula = formula, data = data, family = binomial(), 
                 sample_prior = sample_prior, verbose = verbose, ...))
  })
}


#' Data augmentation for model stability
#' @keywords internal
.augment_data_for_stability <- function(data, formula, verbose = TRUE) {
  
  response_var <- all.vars(formula)[1]
  y <- data[[response_var]]
  
  # Only augment if response is binary
  if (!all(y %in% c(0, 1))) {
    if (verbose) cat("Non-binary response, skipping augmentation\n")
    return(data)
  }
  
  # Determine augmentation size (5% of original data, minimum 2)
  n_aug <- max(2, ceiling(nrow(data) * 0.05))
  
  # Create augmented observations
  aug_indices <- sample(nrow(data), n_aug, replace = TRUE)
  aug_data <- data[aug_indices, ]
  
  # Balance the augmented responses
  aug_data[[response_var]] <- rep(c(0, 1), length.out = n_aug)
  
  # Shrink continuous predictors towards their means (mild regularisation)
  predictors <- all.vars(formula)[-1]
  for (pred in predictors) {
    if (pred %in% names(data) && is.numeric(data[[pred]])) {
      pred_mean <- mean(data[[pred]], na.rm = TRUE)
      shrinkage <- 0.9  # Mild shrinkage
      aug_data[[pred]] <- aug_data[[pred]] * shrinkage + pred_mean * (1 - shrinkage)
    }
  }
  
  # Combine original and augmented data
  rbind(data, aug_data)
}

#' @importFrom stats delete.response model.matrix terms as.formula
NULL