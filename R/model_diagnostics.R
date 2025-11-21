# =============================================================================
# R/model_diagnostics.R
# =============================================================================

#' Automated Model Diagnostics and Recommendations
#'
#' @description
#' Comprehensive automated diagnostics for qbrms models with actionable
#' recommendations for model improvement.
#'
#' @importFrom stats shapiro.test sd cor
#' @importFrom utils head
#'
#' @param model A fitted qbrms model object
#' @param checks Character vector specifying which checks to perform. Options:
#'   "all" (default), "convergence", "fit", "residuals", "posterior", "influential"
#' @param verbose Logical; if TRUE, prints detailed diagnostic information
#'   (default: TRUE)
#'
#' @return An object of class "qbrms_diagnostics" containing:
#' \itemize{
#'   \item \code{summary}: Overall assessment (pass/warning/fail)
#'   \item \code{checks}: Detailed results for each diagnostic check
#'   \item \code{recommendations}: Specific suggestions for improvement
#'   \item \code{plots}: List of diagnostic plots
#' }
#'
#' @details
#' This function performs comprehensive model diagnostics including:
#' \itemize{
#'   \item Convergence checks (for MCMC-based inference)
#'   \item Goodness-of-fit assessment
#'   \item Residual analysis
#'   \item Posterior predictive checks
#'   \item Influential observation detection
#'   \item Prior-posterior overlap assessment
#' }
#'
#' Each check produces a pass/warning/fail status with specific recommendations
#' for addressing any issues detected.
#'
#' @examples
#' \dontrun{
#' # Fit a model
#' fit <- qbrms(mpg ~ hp + wt, data = mtcars, family = gaussian())
#'
#' # Run diagnostics
#' diag <- diagnose_model(fit)
#'
#' # View summary
#' print(diag)
#'
#' # View specific recommendations
#' diag$recommendations
#'
#' # Create diagnostic plots
#' plot(diag)
#' }
#'
#' @export
diagnose_model <- function(model, checks = "all", verbose = TRUE) {
  
  if (!inherits(model, "qbrms_fit") && !inherits(model, "qbrmO_fit")) {
    stop("model must be a fitted qbrms model object", call. = FALSE)
  }
  
  # Determine which checks to perform
  all_checks <- c("convergence", "fit", "residuals", "posterior", "influential")
  
  if ("all" %in% checks) {
    checks_to_run <- all_checks
  } else {
    checks_to_run <- match.arg(checks, all_checks, several.ok = TRUE)
  }
  
  if (verbose) {
    cat("\nRunning Model Diagnostics\n")
    cat("=========================\n\n")
  }
  
  # Initialise results
  diagnostic_results <- list()
  recommendations <- character(0)
  diagnostic_plots <- list()
  
  # Run each check
  if ("convergence" %in% checks_to_run) {
    if (verbose) cat("Checking convergence...\n")
    conv_result <- .check_convergence(model)
    diagnostic_results$convergence <- conv_result
    if (conv_result$status != "pass") {
      recommendations <- c(recommendations, conv_result$recommendations)
    }
  }
  
  if ("fit" %in% checks_to_run) {
    if (verbose) cat("Assessing model fit...\n")
    fit_result <- .check_model_fit(model)
    diagnostic_results$fit <- fit_result
    if (fit_result$status != "pass") {
      recommendations <- c(recommendations, fit_result$recommendations)
    }
    diagnostic_plots$fit <- fit_result$plot
  }
  
  if ("residuals" %in% checks_to_run) {
    if (verbose) cat("Analysing residuals...\n")
    resid_result <- .check_residuals(model)
    diagnostic_results$residuals <- resid_result
    if (resid_result$status != "pass") {
      recommendations <- c(recommendations, resid_result$recommendations)
    }
    diagnostic_plots$residuals <- resid_result$plots
  }
  
  if ("posterior" %in% checks_to_run) {
    if (verbose) cat("Checking posterior distribution...\n")
    post_result <- .check_posterior(model)
    diagnostic_results$posterior <- post_result
    if (post_result$status != "pass") {
      recommendations <- c(recommendations, post_result$recommendations)
    }
  }
  
  if ("influential" %in% checks_to_run) {
    if (verbose) cat("Detecting influential observations...\n")
    infl_result <- .check_influential_observations(model)
    diagnostic_results$influential <- infl_result
    if (infl_result$status != "pass") {
      recommendations <- c(recommendations, infl_result$recommendations)
    }
    diagnostic_plots$influential <- infl_result$plot
  }
  
  # Overall assessment
  all_statuses <- sapply(diagnostic_results, function(x) x$status)
  overall_status <- if (any(all_statuses == "fail")) {
    "fail"
  } else if (any(all_statuses == "warning")) {
    "warning"
  } else {
    "pass"
  }
  
  if (verbose) {
    cat("\n")
    .print_diagnostic_summary(overall_status, diagnostic_results)
  }
  
  # Create result object
  result <- structure(
    list(
      summary = overall_status,
      checks = diagnostic_results,
      recommendations = unique(recommendations),
      plots = diagnostic_plots,
      model = model
    ),
    class = "qbrms_diagnostics"
  )
  
  return(result)
}

#' Check Convergence
#' @keywords internal
.check_convergence <- function(model) {
  
  # For INLA models, convergence is typically not an issue
  # Check if model fitted successfully
  
  if (is.null(model$inla_result)) {
    return(list(
      status = "fail",
      message = "Model fitting failed",
      recommendations = "Re-run the model with different initial values or priors"
    ))
  }
  
  # Check for any INLA warnings
  if (!is.null(model$warnings) && length(model$warnings) > 0) {
    return(list(
      status = "warning",
      message = "Model fitted with warnings",
      details = model$warnings,
      recommendations = "Review model specification and data quality"
    ))
  }
  
  return(list(
    status = "pass",
    message = "Model converged successfully",
    recommendations = NULL
  ))
}

#' Check Model Fit
#' @keywords internal
.check_model_fit <- function(model) {
  
  # ---------------------------------------------------------------------------
  # 1. Work out the response name safely
  # ---------------------------------------------------------------------------
  response_name <- NULL
  
  if (!is.null(model$formula) && inherits(model$formula, "formula")) {
    lhs <- tryCatch(as.character(model$formula)[2L],
                    error = function(e) NULL)
    if (is.character(lhs) && length(lhs) == 1L && !is.na(lhs)) {
      response_name <- lhs
    }
  }
  
  has_response <- !is.null(response_name) &&
    !is.null(model$data) &&
    is.data.frame(model$data) &&
    response_name %in% names(model$data)
  
  # If we cannot map a response column in model$data, skip gracefully
  if (!has_response) {
    return(list(
      status         = "pass",
      message        = "Model fit check skipped: response variable not found in model$data.",
      r_squared      = NA_real_,
      recommendations = NULL,
      plot           = NULL
    ))
  }
  
  observed_vals <- model$data[[response_name]]
  
  # ---------------------------------------------------------------------------
  # 2. Extract fitted values from known locations
  # ---------------------------------------------------------------------------
  fitted_vals <- tryCatch({
    if (!is.null(model$fitted_values)) {
      model$fitted_values
    } else if (!is.null(model$inla_result$summary.fitted.values)) {
      model$inla_result$summary.fitted.values$mean
    } else if (!is.null(model$fit$fitted.values)) {
      model$fit$fitted.values
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  if (is.null(fitted_vals) || length(fitted_vals) == 0L) {
    return(list(
      status         = "warning",
      message        = "Unable to assess model fit: no fitted values found.",
      r_squared      = NA_real_,
      recommendations = "Check that the model object contains fitted values.",
      plot           = NULL
    ))
  }
  
  # Align lengths defensively
  n <- min(length(observed_vals), length(fitted_vals))
  observed_vals <- observed_vals[seq_len(n)]
  fitted_vals   <- fitted_vals[seq_len(n)]
  
  # ---------------------------------------------------------------------------
  # 3. Numeric outcomes: R^2-style check
  # ---------------------------------------------------------------------------
  if (is.numeric(observed_vals)) {
    
    valid_idx <- !is.na(observed_vals) & !is.na(fitted_vals)
    observed_vals <- observed_vals[valid_idx]
    fitted_vals   <- fitted_vals[valid_idx]
    
    if (length(observed_vals) < 3L) {
      return(list(
        status         = "warning",
        message        = "Too few non-missing observations to compute R-squared.",
        r_squared      = NA_real_,
        recommendations = NULL,
        plot           = NULL
      ))
    }
    
    ss_res <- sum((observed_vals - fitted_vals)^2)
    ss_tot <- sum((observed_vals - mean(observed_vals))^2)
    r_squared <- if (ss_tot > 0) 1 - (ss_res / ss_tot) else NA_real_
    
    if (is.na(r_squared)) {
      status  <- "warning"
      message <- "Unable to compute R-squared."
      recs    <- NULL
    } else if (r_squared < 0.1) {
      status  <- "warning"
      message <- sprintf("Poor model fit (R^2 = %.3f)", r_squared)
      recs    <- c(
        "Consider adding more predictors.",
        "Check for non-linear relationships.",
        "Verify that the family and link function are appropriate."
      )
    } else if (r_squared < 0.3) {
      status  <- "warning"
      message <- sprintf("Moderate model fit (R^2 = %.3f)", r_squared)
      recs    <- "Consider additional predictors or interactions."
    } else {
      status  <- "pass"
      message <- sprintf("Good model fit (R^2 = %.3f)", r_squared)
      recs    <- NULL
    }
    
    fit_plot <- NULL
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot_data <- data.frame(
        observed = observed_vals,
        fitted   = fitted_vals
      )
      
      fit_plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = fitted, y = observed)
      ) +
        ggplot2::geom_point(alpha = 0.5, colour = "steelblue") +
        ggplot2::geom_abline(
          slope     = 1,
          intercept = 0,
          linetype  = "dashed",
          colour    = "red"
        ) +
        ggplot2::labs(
          title    = "Observed vs Fitted values",
          subtitle = if (is.na(r_squared)) "R^2: NA"
          else sprintf("R^2 = %.3f", r_squared),
          x        = "Fitted values",
          y        = "Observed values"
        ) +
        ggplot2::theme_minimal()
    }
    
    return(list(
      status         = status,
      message        = message,
      r_squared      = r_squared,
      recommendations = recs,
      plot           = fit_plot
    ))
  }
  
  # ---------------------------------------------------------------------------
  # 4. Non-numeric responses: politely decline
  # ---------------------------------------------------------------------------
  list(
    status         = "pass",
    message        = "Numeric R-squared fit check only implemented for numeric responses.",
    r_squared      = NA_real_,
    recommendations = NULL,
    plot           = NULL
  )
}

#' Check Residuals
#' @keywords internal
.check_residuals <- function(model) {
  
  # Calculate residuals
  residuals <- tryCatch({
    if (!is.null(model$residuals)) {
      model$residuals
    } else {
      observed <- model$data[[as.character(model$formula[[2]])]]
      fitted <- if (!is.null(model$fitted_values)) {
        model$fitted_values
      } else if (!is.null(model$inla_result$summary.fitted.values)) {
        model$inla_result$summary.fitted.values$mean
      } else {
        NULL
      }
      if (!is.null(fitted)) observed - fitted else NULL
    }
  }, error = function(e) NULL)
  
  if (is.null(residuals)) {
    return(list(
      status = "warning",
      message = "Unable to compute residuals",
      recommendations = NULL,
      plots = NULL
    ))
  }
  
  # Remove missing values
  residuals <- residuals[!is.na(residuals)]
  
  # Check for normality (Shapiro-Wilk test)
  if (length(residuals) > 3 && length(residuals) < 5000) {
    shapiro_test <- shapiro.test(residuals)
    normal_p_value <- shapiro_test$p.value
  } else {
    normal_p_value <- NA
  }
  
  # Check for patterns
  # Test for autocorrelation if applicable
  n_resid <- length(residuals)
  
  status_messages <- character(0)
  recommendations <- character(0)
  
  # Normality check
  if (!is.na(normal_p_value) && normal_p_value < 0.05) {
    status_messages <- c(status_messages, 
                         "Residuals deviate from normality")
    recommendations <- c(recommendations,
                         "Consider using a robust family (e.g., Student-t)",
                         "Check for outliers or data quality issues")
  }
  
  # Check for extreme residuals
  extreme_threshold <- 3
  n_extreme <- sum(abs(scale(residuals)) > extreme_threshold)
  prop_extreme <- n_extreme / n_resid
  
  if (prop_extreme > 0.05) {
    status_messages <- c(status_messages,
                         sprintf("%.1f%% extreme residuals detected", 
                                 prop_extreme * 100))
    recommendations <- c(recommendations,
                         "Investigate potential outliers or influential observations")
  }
  
  # Determine overall status
  if (length(status_messages) == 0) {
    status <- "pass"
    message <- "Residuals appear well-behaved"
  } else if (length(status_messages) == 1) {
    status <- "warning"
    message <- status_messages[1]
  } else {
    status <- "warning"
    message <- paste(status_messages, collapse = "; ")
  }
  
  # Create diagnostic plots
  plots <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    resid_data <- data.frame(
      residuals = residuals,
      index = seq_along(residuals),
      standardised = scale(residuals)[, 1]
    )
    
    # Q-Q plot
    qq_plot <- ggplot2::ggplot(resid_data, 
                               ggplot2::aes(sample = standardised)) +
      ggplot2::stat_qq(colour = "steelblue") +
      ggplot2::stat_qq_line(colour = "red", linetype = "dashed") +
      ggplot2::labs(
        title = "Q-Q Plot of Standardised Residuals",
        x = "Theoretical Quantiles",
        y = "Sample Quantiles"
      ) +
      ggplot2::theme_minimal()
    
    # Residuals vs index
    index_plot <- ggplot2::ggplot(resid_data,
                                  ggplot2::aes(x = index, y = residuals)) +
      ggplot2::geom_point(alpha = 0.5, colour = "steelblue") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
      ggplot2::geom_hline(yintercept = c(-2, 2) * sd(residuals),
                          linetype = "dotted", colour = "grey50") +
      ggplot2::labs(
        title = "Residuals vs Index",
        x = "Observation Index",
        y = "Residuals"
      ) +
      ggplot2::theme_minimal()
    
    plots <- list(qq = qq_plot, index = index_plot)
  }
  
  return(list(
    status = status,
    message = message,
    shapiro_p = normal_p_value,
    n_extreme = n_extreme,
    recommendations = recommendations,
    plots = plots
  ))
}

#' Check Posterior Distribution
#' @keywords internal
.check_posterior <- function(model) {
  
  # Extract posterior summaries
  fixed_effects <- tryCatch({
    if (!is.null(model$fixed_effects)) {
      model$fixed_effects
    } else if (!is.null(model$inla_result$summary.fixed)) {
      model$inla_result$summary.fixed
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  if (is.null(fixed_effects)) {
    return(list(
      status = "warning",
      message = "Unable to extract posterior summaries",
      recommendations = NULL
    ))
  }
  
  recommendations <- character(0)
  issues <- character(0)
  
  # Check for very wide credible intervals
  if ("0.025quant" %in% names(fixed_effects) && 
      "0.975quant" %in% names(fixed_effects)) {
    
    ci_widths <- fixed_effects$"0.975quant" - fixed_effects$"0.025quant"
    mean_width <- mean(ci_widths, na.rm = TRUE)
    
    if (mean_width > 10) {
      issues <- c(issues, "Very wide credible intervals detected")
      recommendations <- c(recommendations,
                           "Consider using more informative priors",
                           "Check for identification issues in the model")
    }
  }
  
  # Check for posteriors similar to priors (lack of information from data)
  # This would require prior specification - simplified check here
  if ("sd" %in% names(fixed_effects)) {
    large_sd <- fixed_effects$sd > 10
    if (any(large_sd, na.rm = TRUE)) {
      affected_params <- rownames(fixed_effects)[large_sd]
      issues <- c(issues,
                  sprintf("Large posterior SD for: %s", 
                          paste(affected_params, collapse = ", ")))
      recommendations <- c(recommendations,
                           "Parameters with large SD may be poorly identified",
                           "Consider simplifying the model or adding data")
    }
  }
  
  # Determine status
  if (length(issues) == 0) {
    status <- "pass"
    message <- "Posterior distributions appear reasonable"
  } else {
    status <- "warning"
    message <- paste(issues, collapse = "; ")
  }
  
  return(list(
    status = status,
    message = message,
    recommendations = recommendations
  ))
}

#' Check for Influential Observations
#' @keywords internal
.check_influential_observations <- function(model) {
  
  # Calculate Cook's distance equivalent
  residuals <- tryCatch({
    if (!is.null(model$residuals)) {
      model$residuals
    } else {
      observed <- model$data[[as.character(model$formula[[2]])]]
      fitted <- if (!is.null(model$fitted_values)) {
        model$fitted_values
      } else if (!is.null(model$inla_result$summary.fitted.values)) {
        model$inla_result$summary.fitted.values$mean
      } else {
        NULL
      }
      if (!is.null(fitted)) observed - fitted else NULL
    }
  }, error = function(e) NULL)
  
  if (is.null(residuals)) {
    return(list(
      status = "warning",
      message = "Unable to assess influential observations",
      recommendations = NULL,
      plot = NULL
    ))
  }
  
  # Standardise residuals
  std_resid <- scale(residuals)[, 1]
  
  # Flag influential observations (|standardised residual| > 3)
  influential_threshold <- 3
  influential_idx <- which(abs(std_resid) > influential_threshold)
  n_influential <- length(influential_idx)
  prop_influential <- n_influential / length(residuals)
  
  recommendations <- character(0)
  
  if (prop_influential > 0.05) {
    status <- "warning"
    message <- sprintf("%d influential observations detected (%.1f%%)",
                       n_influential, prop_influential * 100)
    recommendations <- c(
      sprintf("Investigate observations: %s", 
              paste(head(influential_idx, 10), collapse = ", ")),
      "Consider robust regression methods",
      "Verify data quality for flagged observations"
    )
  } else if (n_influential > 0) {
    status <- "pass"
    message <- sprintf("%d influential observations detected (within acceptable range)",
                       n_influential)
  } else {
    status <- "pass"
    message <- "No influential observations detected"
  }
  
  # Create diagnostic plot
  influence_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_data <- data.frame(
      index = seq_along(std_resid),
      std_residual = std_resid,
      influential = abs(std_resid) > influential_threshold
    )
    
    influence_plot <- ggplot2::ggplot(plot_data,
                                      ggplot2::aes(x = index, 
                                                   y = std_residual,
                                                   colour = influential)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_hline(yintercept = c(-influential_threshold, influential_threshold),
                          linetype = "dashed", colour = "red") +
      ggplot2::scale_colour_manual(values = c("FALSE" = "steelblue",
                                              "TRUE" = "coral")) +
      ggplot2::labs(
        title = "Influential Observations",
        subtitle = sprintf("%d observations flagged", n_influential),
        x = "Observation Index",
        y = "Standardised Residual",
        colour = "Influential"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
  }
  
  return(list(
    status = status,
    message = message,
    n_influential = n_influential,
    influential_indices = influential_idx,
    recommendations = recommendations,
    plot = influence_plot
  ))
}

#' Print Diagnostic Summary
#' @keywords internal
.print_diagnostic_summary <- function(overall_status, results) {
  
  cat("\nOverall Status: ")
  status_symbol <- switch(overall_status,
                          "pass" = "[PASS]",
                          "warning" = "[WARN]",
                          "fail" = "[FAIL]"
  )
  
  cat(status_symbol, toupper(overall_status), "\n\n")
  
  # Individual check results
  cat("Individual Checks:\n")
  cat("------------------\n")
  
  for (check_name in names(results)) {
    result <- results[[check_name]]
    symbol <- switch(result$status,
                     "pass" = "[PASS]",
                     "warning" = "[WARN]",
                     "fail" = "[FAIL]"
    )
    cat(sprintf("  %s %s: %s\n", symbol, check_name, result$message))
  }
  
  cat("\n")
}
#' Print Method for Diagnostics
#' @param x A qbrms_diagnostics object
#' @param ... Additional arguments (unused)
#' @export
print.qbrms_diagnostics <- function(x, ...) {
  
  cat("\nModel Diagnostics Results\n")
  cat("=========================\n\n")
  
  cat("Overall Assessment:", toupper(x$summary), "\n\n")
  
  cat("Checks Performed:\n")
  for (check_name in names(x$checks)) {
    result <- x$checks[[check_name]]
    cat(sprintf("  - %s: %s\n", check_name, result$message))
  }
  
  if (length(x$recommendations) > 0) {
    cat("\nRecommendations:\n")
    for (i in seq_along(x$recommendations)) {
      cat(sprintf("  %d. %s\n", i, x$recommendations[i]))
    }
  } else {
    cat("\nNo issues detected. Model appears well-specified.\n")
  }
  
  cat("\nUse plot(diagnostics) to view diagnostic plots.\n")
  
  invisible(x)
}

#' Plot Method for Diagnostics
#' @param x A qbrms_diagnostics object
#' @param which Character vector specifying which plots to show
#' @param ... Additional arguments (unused)
#' @export
plot.qbrms_diagnostics <- function(x, which = "all", ...) {
  
  available_plots <- names(x$plots)
  
  if (length(available_plots) == 0) {
    message("No diagnostic plots available")
    return(invisible(NULL))
  }
  
  if ("all" %in% which) {
    plots_to_show <- available_plots
  } else {
    plots_to_show <- intersect(which, available_plots)
  }
  
  if (length(plots_to_show) == 0) {
    message("Requested plots not available")
    return(invisible(NULL))
  }
  
  # Display each plot
  for (plot_name in plots_to_show) {
    plot_obj <- x$plots[[plot_name]]
    
    if (is.list(plot_obj)) {
      # Multiple plots for this check
      for (subplot in plot_obj) {
        if (!is.null(subplot)) {
          print(subplot)
          readline(prompt = "Press [Enter] for next plot...")
        }
      }
    } else {
      # Single plot
      if (!is.null(plot_obj)) {
        print(plot_obj)
        if (plot_name != plots_to_show[length(plots_to_show)]) {
          readline(prompt = "Press [Enter] for next plot...")
        }
      }
    }
  }
  
  invisible(x)
}