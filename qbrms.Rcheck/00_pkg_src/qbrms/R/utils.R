# =============================================================================
# R/utils.R
# =============================================================================

#' Null Coalescing Operator
#'
#' @description
#' Returns the first non-NULL value.
#'
#' @name grapes-or-or-grapes
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
#' @keywords internal
extract_family_name <- function(inla_family) {
  if (is.list(inla_family)) {
    return(inla_family$family %||% "gaussian")
  } else if (is.character(inla_family) && length(inla_family) >= 1) {
    return(inla_family[1])
  } else {
    return("gaussian")
  }
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
create_dummy_data_for_priors <- function(formula, data, n_dummy, family_name, verbose = TRUE) { # MODIFIED: Added family_name
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
generate_prior_predictions_simple <- function(model, data, formula, family_name, ndraws, n_cats = NULL) { # MODIFIED: Added n_cats argument
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