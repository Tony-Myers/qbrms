# =============================================================================
# R/family_conversion.R
# =============================================================================

#' Family conversion utilities for qbrms package
#' 
#' @description
#' This module provides helper functions for family specification validation
#' and prior predictive sampling. The main family conversion is handled in
#' families.R to avoid duplication.
#' 
#' @name family_conversion
#' @keywords internal
NULL

# =============================================================================
# PRIOR PREDICTIVE SAMPLING
# =============================================================================

#' Generate Prior Predictive Samples
#'
#' @description
#' Generate samples from the prior predictive distribution for qbrms models.
#' Returns an ndraws x nobs matrix for compatibility with qbrms structure.
#'
#' @param formula Model formula
#' @param data Data frame containing model variables
#' @param family Model family specification
#' @param prior Prior specifications (currently uses default priors)
#' @param ndraws Number of draws from prior predictive distribution
#' @param verbose Logical; print progress messages
#' @param ... Additional arguments (currently ignored)
#'
#' @return Matrix of prior predictive samples (ndraws x nrow(data))
#'
#' @details
#' Uses simple default priors:
#' \itemize{
#'   \item Intercept: Normal(0, 2.5)
#'   \item Other coefficients: Normal(0, 1.0)
#' }
#' 
#' For families other than gaussian, binomial, and poisson, falls back to 
#' Gaussian-like sampling.
#'
#' @keywords internal
generate_prior_predictive_samples <- function(formula, data, family = gaussian(),
                                              prior = NULL, ndraws = 100,
                                              verbose = TRUE, ...) {
  
  # Helper function to remove random effects for design matrix
  drop_random_effects <- function(fml) {
    f_txt <- paste(deparse(fml), collapse = "")
    lhs <- sub("~.*", "", f_txt)
    rhs <- sub(".*~", "", f_txt)
    
    # Remove random effects terms
    rhs_clean <- gsub("\\([^\\)]*\\|[^\\)]*\\)", "", rhs)
    rhs_clean <- gsub("\\+\\s*\\+", "+", rhs_clean)  # Fix double +
    rhs_clean <- gsub("^\\s*\\+\\s*", "", rhs_clean)  # Remove leading +
    rhs_clean <- gsub("\\s*\\+\\s*$", "", rhs_clean)  # Remove trailing +
    
    if (nchar(trimws(rhs_clean)) == 0) {
      rhs_clean <- "1"
    }
    
    as.formula(paste(lhs, "~", rhs_clean))
  }
  
  # Get design matrix  
  tryCatch({
    fixed_formula <- drop_random_effects(formula)
    X <- model.matrix(fixed_formula, data = data)
  }, error = function(e) {
    if (verbose) cat("Warning: Could not create design matrix, using intercept only\n")
    X <- matrix(1, nrow = nrow(data), ncol = 1)
    colnames(X) <- "(Intercept)"
  })
  
  n_coef <- ncol(X)
  n_obs <- nrow(data)
  
  # Default priors: intercept gets wider prior
  prior_means <- rep(0, n_coef)
  prior_sds <- ifelse(colnames(X) == "(Intercept)", 2.5, 1.0)
  
  # Generate coefficient samples
  coef_samples <- matrix(NA, nrow = ndraws, ncol = n_coef)
  for (i in seq_len(n_coef)) {
    coef_samples[, i] <- rnorm(ndraws, prior_means[i], prior_sds[i])
  }
  
  # Generate predictions based on family
  family_name <- tryCatch({
    extract_family_name(convert_family_to_inla(family))
  }, error = function(e) "gaussian")
  
  pred_samples <- matrix(NA, nrow = ndraws, ncol = n_obs)
  
  for (i in seq_len(ndraws)) {
    linear_pred <- as.numeric(X %*% coef_samples[i, ])
    
    # Family-specific response generation
    if (family_name %in% c("gaussian", "normal")) {
      sigma <- abs(rnorm(1, 0, 0.5)) + 0.1  # Ensure positive
      pred_samples[i, ] <- rnorm(n_obs, linear_pred, sigma)
      
    } else if (family_name == "binomial") {
      probs <- plogis(linear_pred)
      pred_samples[i, ] <- rbinom(n_obs, 1, probs)
      
    } else if (family_name == "poisson") {
      lambda <- exp(linear_pred)
      lambda <- pmin(lambda, 20)  # Prevent extremely large values
      pred_samples[i, ] <- rpois(n_obs, lambda)
      
    } else {
      # Fallback to Gaussian-like for other families
      sigma <- abs(rnorm(1, 0, 0.5)) + 0.1
      pred_samples[i, ] <- rnorm(n_obs, linear_pred, sigma)
    }
  }
  
  return(pred_samples)
}