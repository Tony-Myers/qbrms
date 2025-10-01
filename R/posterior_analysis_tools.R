# =============================================================================
# R/bayesian_analysis.R
# Bayes factors, probability of direction, ROPE, HDI, and EMMs, p_significance
# =============================================================================

#' @title Bayesian Analysis Functions (qbrms)
#' @description Posterior analysis tools for qbrms models: Bayes factors,
#'   probability of direction, ROPE, HDI, and estimated marginal means.
#' @name qbrms_bayesian_analysis
NULL


# =============================================================================
# SECTION 1: BAYES FACTORS AND HYPOTHESIS TESTING
# =============================================================================

#' Bayesian Hypothesis Testing (very simple approximations)
#'
#' @description
#' Compute a crude Bayes factor for a point, interval, or comparison hypothesis
#' using approximate posterior draws recovered from a `qbrms_fit`.
#' This is deliberately simple and intended for exploratory use.
#'
#' @param object A \code{qbrms_fit} object.
#' @param hypothesis Character string, for example "Intercept > 0",
#'   "b_x = 0", or "b_x > 0.2".
#' @param prior Optional prior information (unused here, kept for API compatibility).
#' @param null Numeric null value for point tests (default 0).
#' @param direction One of "two-sided", "greater", "less" (kept for API compatibility).
#' @param rope Optional numeric length-2 vector \code{c(lower, upper)} to
#'   define a ROPE for point tests.
#' @param nsim Number of posterior draws to simulate from the fitted summary.
#' @param verbose Logical; print progress information.
#'
#' @return An object of class \code{qbrms_bayesfactor}.
#' @export
bayesfactor <- function(object,
                        hypothesis,
                        prior = NULL,
                        null = 0,
                        direction = "two-sided",
                        rope = NULL,
                        nsim = 4000,
                        verbose = TRUE) {
  if (!inherits(object, "qbrms_fit")) {
    stop("object must be a qbrms_fit object")
  }
  
  samples <- .ba_extract_parameter_samples(object, nsim)
  info    <- .ba_parse_hypothesis(hypothesis, names(samples))
  
  if (verbose) {
    cat("Bayesian Hypothesis Test\n")
    cat("Hypothesis: ", hypothesis, "\n", sep = "")
    cat("Parameters involved: ",
        if (length(info$parameters)) paste(info$parameters, collapse = ", ") else "(none)",
        "\n", sep = "")
  }
  
  res <- switch(
    info$type,
    "point"      = .ba_compute_point_bf(samples, info, null, rope),
    "interval"   = .ba_compute_interval_bf(samples, info),
    "comparison" = .ba_compute_comparison_bf(samples, info),
    stop("Unrecognised hypothesis type")
  )
  
  out <- list(
    hypothesis     = hypothesis,
    bayes_factor   = res$bf,
    log_bf         = log(res$bf),
    evidence       = .ba_interpret_bf(res$bf),
    posterior_prob = res$posterior_prob,
    prior_prob     = res$prior_prob,
    details        = res$details,
    nsim           = nsim
  )
  class(out) <- c("qbrms_bayesfactor", "list")
  if (verbose) {
    cat("\nResults:\n")
    cat("Bayes Factor (BF10): ", round(out$bayes_factor, 3), "\n", sep = "")
    cat("Evidence: ", out$evidence, "\n", sep = "")
    cat("Posterior probability of H1: ", round(out$posterior_prob, 3), "\n", sep = "")
  }
  out
}

#' @export
#' @method print qbrms_bayesfactor
print.qbrms_bayesfactor <- function(x, ...) {
  cat("Bayesian Hypothesis Test\n")
  cat("========================\n\n")
  cat("Hypothesis: ", x$hypothesis, "\n", sep = "")
  cat("Bayes Factor (BF10): ", round(x$bayes_factor, 3), "\n", sep = "")
  cat("Log Bayes Factor: ", round(x$log_bf, 3), "\n", sep = "")
  cat("Evidence: ", x$evidence, "\n", sep = "")
  cat("Posterior probability of H1: ", round(x$posterior_prob, 3), "\n", sep = "")
  cat("Sample size: ", x$nsim, "\n", sep = "")
  invisible(x)
}

# =============================================================================
# SECTION 2: PROBABILITY OF DIRECTION, ROPE, HDI
# =============================================================================

#' Probability of Direction (pd)
#'
#' @description
#' Estimate the probability that a parameter is strictly positive
#' (or strictly negative) under the posterior, based on simulated draws
#' from a \code{qbrms_fit}.
#'
#' @param object A \code{qbrms_fit} object.
#' @param parameters Optional character vector of parameter names. If \code{NULL},
#'   all fixed-effect coefficients are used.
#' @param nsim Number of draws to simulate from the fitted summary.
#' @param null Numeric value defining the reference for direction (default 0).
#'
#' @return A data frame of class \code{qbrms_p_direction}.
#' @export
p_direction <- function(object, parameters = NULL, nsim = 4000, null = 0) {
  if (!inherits(object, "qbrms_fit")) {
    stop("object must be a qbrms_fit object")
  }
  
  draws_list <- .ba_extract_parameter_samples(object, nsim)
  if (is.null(parameters)) parameters <- names(draws_list)
  
  out <- data.frame(
    Parameter = parameters,
    PD        = NA_real_,
    Direction = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(parameters)) {
    p <- parameters[i]
    if (!is.null(draws_list[[p]])) {
      s  <- draws_list[[p]]
      nP <- sum(s >  null)
      nN <- sum(s <  null)
      pd <- max(nP, nN) / length(s)
      out$PD[i]        <- pd
      out$Direction[i] <- if (nP > nN) "positive" else "negative"
    }
  }
  
  class(out) <- c("qbrms_p_direction", "data.frame")
  out
}

#' @export
#' @method print qbrms_p_direction
print.qbrms_p_direction <- function(x, ...) {
  cat("Probability of Direction (pd)\n")
  cat("============================\n\n")
  print.data.frame(x, row.names = FALSE)
  cat("\nInterpretation:\n")
  cat("pd > 0.95: Strong directional evidence\n")
  cat("pd > 0.90: Moderate directional evidence\n")
  invisible(x)
}

#' ROPE analysis
#'
#' @description
#' Compute the posterior mass inside a Region Of Practical Equivalence (ROPE)
#' for selected parameters.
#'
#' @param object A \code{qbrms_fit}.
#' @param parameters Optional character vector of parameter names. If \code{NULL},
#'   all fixed-effect coefficients are used.
#' @param rope Numeric length-2 vector \code{c(lower, upper)}.
#' @param nsim Number of posterior draws to simulate.
#'
#' @return A data frame of class \code{qbrms_rope}.
#' @export
rope_analysis <- function(object, parameters = NULL, rope = c(-0.1, 0.1), nsim = 4000) {
  if (!inherits(object, "qbrms_fit")) {
    stop("object must be a qbrms_fit object")
  }
  if (length(rope) != 2) stop("rope must be c(lower, upper)")
  
  draws_list <- .ba_extract_parameter_samples(object, nsim)
  if (is.null(parameters)) parameters <- names(draws_list)
  
  out <- data.frame(
    Parameter      = parameters,
    ROPE_low       = rope[1],
    ROPE_high      = rope[2],
    ROPE_percentage= NA_real_,
    Outside_ROPE   = NA_real_,
    Conclusion     = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(parameters)) {
    p <- parameters[i]
    s <- draws_list[[p]]
    if (!is.null(s)) {
      in_rope <- mean(s >= rope[1] & s <= rope[2])
      out$ROPE_percentage[i] <- in_rope
      out$Outside_ROPE[i]    <- 1 - in_rope
      out$Conclusion[i]      <- if (in_rope > 0.95) {
        "Practically equivalent to null"
      } else if ((1 - in_rope) > 0.95) {
        "Practically significant"
      } else {
        "Inconclusive"
      }
    }
  }
  
  class(out) <- c("qbrms_rope", "data.frame")
  out
}

#' @export
#' @method print qbrms_rope
print.qbrms_rope <- function(x, ...) {
  cat("ROPE Analysis Results\n")
  cat("=====================\n\n")
  print.data.frame(x, row.names = FALSE, digits = 3)
  invisible(x)
}

#' Highest Density Interval (HDI)
#'
#' @description
#' Compute highest density intervals for parameters based on simulated
#' posterior draws from a \code{qbrms_fit}.
#'
#' @param object A \code{qbrms_fit}.
#' @param parameters Optional character vector; default uses all fixed effects.
#' @param prob Probability mass for the interval (default 0.95).
#' @param nsim Number of draws to simulate.
#'
#' @return A data frame of class \code{qbrms_hdi}.
#' @export
hdi <- function(object, parameters = NULL, prob = 0.95, nsim = 4000) {
  if (!inherits(object, "qbrms_fit")) {
    stop("object must be a qbrms_fit object")
  }
  
  draws_list <- .ba_extract_parameter_samples(object, nsim)
  if (is.null(parameters)) parameters <- names(draws_list)
  
  out <- data.frame(
    Parameter  = parameters,
    HDI_lower  = NA_real_,
    HDI_upper  = NA_real_,
    HDI_width  = NA_real_,
    Probability= prob,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(parameters)) {
    p <- parameters[i]
    s <- draws_list[[p]]
    if (!is.null(s)) {
      h <- .ba_compute_hdi(s, prob)
      out$HDI_lower[i] <- h[1]
      out$HDI_upper[i] <- h[2]
      out$HDI_width[i] <- h[2] - h[1]
    }
  }
  
  class(out) <- c("qbrms_hdi", "data.frame")
  out
}

#' @export
#' @method print qbrms_hdi
print.qbrms_hdi <- function(x, ...) {
  cat("Highest Density Intervals\n")
  cat("=========================\n\n")
  print.data.frame(x, row.names = FALSE, digits = 3)
  invisible(x)
}

# =============================================================================
# SECTION 3: Probability of Practical Significance
# =============================================================================
#' Probability of Practical Significance (Enhanced bayestestR-style)
#'
#' @description
#' Compute the probability that each parameter is above a threshold 
#' in the median's direction, similar to bayestestR::p_significance().
#' This represents the proportion of the posterior distribution that 
#' indicates a "significant" effect in the median's direction.
#'
#' @param object A \code{qbrms_fit} object.
#' @param parameters Optional character vector of parameter names; if \code{NULL},
#'   all fixed-effect coefficients are used.
#' @param threshold The threshold value that separates significant from negligible effect:
#'   - \code{"default"}: Uses 0.1 as threshold range around zero
#'   - A single numeric value (e.g., 0.1): Creates symmetric range around zero (-0.1, 0.1)
#'   - A numeric vector of length two (e.g., c(-0.2, 0.1)): Asymmetric threshold
#'   - A list of numeric vectors: Each vector corresponds to a parameter
#'   - A named list: Names correspond to parameter names
#' @param nsim Number of draws to simulate for the approximation.
#' @param verbose Logical; print progress information.
#'
#' @return A data frame of class \code{qbrms_p_significance} with columns
#'   \code{Parameter}, \code{ps}, \code{Median}, \code{CI_low}, \code{CI_high},
#'   \code{Threshold_low}, \code{Threshold_high}, and \code{Interpretation}.
#' @export
p_significance <- function(object,
                           parameters = NULL,
                           threshold = "default",
                           nsim = 1000,
                           verbose = TRUE) {
  if (!inherits(object, "qbrms_fit")) {
    stop("object must be a qbrms_fit object")
  }
  
  # Extract posterior samples
  if (exists(".ba_extract_parameter_samples", mode = "function")) {
    draws_list <- .ba_extract_parameter_samples(object, nsim)
  } else {
    # Fallback sampler from summary means/sds or coef()
    draws_list <- list()
    if (!is.null(object$fit$summary.fixed)) {
      m  <- object$fit$summary.fixed
      mu <- m[, "mean"]; sdv <- m[, "sd"]
      nms <- rownames(m)
      for (i in seq_along(nms)) {
        if (is.finite(mu[i]) && is.finite(sdv[i]) && sdv[i] > 0) {
          draws_list[[nms[i]]] <- stats::rnorm(nsim, mu[i], sdv[i])
        }
      }
    } else {
      cf <- tryCatch(stats::coef(object), error = function(e) NULL)
      if (!is.null(cf)) {
        for (i in seq_along(cf)) {
          nm <- names(cf)[i]
          sd <- abs(cf[i]) * 0.1 + 0.01
          draws_list[[nm]] <- stats::rnorm(nsim, cf[i], sd)
        }
      }
    }
    if (!length(draws_list)) stop("Could not extract posterior samples from model")
  }
  
  if (is.null(parameters)) parameters <- names(draws_list)
  
  # Process threshold argument
  thresholds <- .process_threshold(threshold, parameters, draws_list)
  
  if (verbose) {
    cat("Computing probability of practical significance...\n")
    cat("Parameters:", length(parameters), "\n")
  }
  
  # Initialize results data frame
  out <- data.frame(
    Parameter        = parameters,
    ps               = NA_real_,
    Median           = NA_real_,
    CI_low           = NA_real_,
    CI_high          = NA_real_,
    Threshold_low    = NA_real_,
    Threshold_high   = NA_real_,
    Interpretation   = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Compute for each parameter
  for (i in seq_along(parameters)) {
    param <- parameters[i]
    samples <- draws_list[[param]]
    
    if (!is.null(samples) && length(samples) > 0) {
      # Get threshold for this parameter
      thresh <- thresholds[[param]]
      
      # Calculate median and credible interval
      median_val <- stats::median(samples)
      ci_vals <- stats::quantile(samples, c(0.025, 0.975))
      
      # Determine direction based on median
      if (median_val > 0) {
        # Positive direction: probability above upper threshold
        ps_val <- mean(samples > thresh[2])
      } else {
        # Negative direction: probability below lower threshold  
        ps_val <- mean(samples < thresh[1])
      }
      
      # Fill results
      out$ps[i]             <- ps_val
      out$Median[i]         <- median_val
      out$CI_low[i]         <- ci_vals[1]
      out$CI_high[i]        <- ci_vals[2]
      out$Threshold_low[i]  <- thresh[1]
      out$Threshold_high[i] <- thresh[2]
      out$Interpretation[i] <- .interpret_ps(ps_val)
      
      if (verbose) {
        cat("  ", param, ": ps =", round(ps_val, 3), 
            "(median =", round(median_val, 3), ")\n")
      }
    }
  }
  
  class(out) <- c("qbrms_p_significance", "data.frame")
  return(out)
}

# Helper function to process threshold argument
.process_threshold <- function(threshold, parameters, draws_list) {
  n_params <- length(parameters)
  
  # Handle different threshold specifications
  if (identical(threshold, "default")) {
    # Default: symmetric +/-0.1 for all parameters
    default_thresh <- c(-0.1, 0.1)
    thresholds <- replicate(n_params, default_thresh, simplify = FALSE)
    names(thresholds) <- parameters
    
  } else if (is.numeric(threshold)) {
    if (length(threshold) == 1) {
      # Single value: create symmetric range
      thresh_val <- abs(threshold)
      thresh_range <- c(-thresh_val, thresh_val)
      thresholds <- replicate(n_params, thresh_range, simplify = FALSE)
      names(thresholds) <- parameters
      
    } else if (length(threshold) == 2) {
      # Two values: use as lower and upper bounds for all parameters
      thresholds <- replicate(n_params, threshold, simplify = FALSE)
      names(thresholds) <- parameters
      
    } else {
      stop("threshold must be 'default', a single value, or a vector of length 2")
    }
    
  } else if (is.list(threshold)) {
    # List specification
    thresholds <- list()
    
    if (is.null(names(threshold))) {
      # Unnamed list: assume same order as parameters
      if (length(threshold) != n_params) {
        stop("Length of threshold list must match number of parameters")
      }
      for (i in seq_along(parameters)) {
        thresholds[[parameters[i]]] <- .validate_threshold_vector(threshold[[i]])
      }
      
    } else {
      # Named list: match by names
      for (param in parameters) {
        if (param %in% names(threshold)) {
          thresholds[[param]] <- .validate_threshold_vector(threshold[[param]])
        } else {
          # Use default for unspecified parameters
          thresholds[[param]] <- c(-0.1, 0.1)
        }
      }
    }
    
  } else {
    stop("threshold must be 'default', numeric, or a list")
  }
  
  return(thresholds)
}

# Helper to validate individual threshold vectors
.validate_threshold_vector <- function(thresh) {
  if (length(thresh) == 1) {
    thresh_val <- abs(thresh)
    return(c(-thresh_val, thresh_val))
  } else if (length(thresh) == 2) {
    return(sort(thresh))  # Ensure lower <= upper
  } else {
    stop("Each threshold vector must have length 1 or 2")
  }
}

# Helper to interpret ps values
.interpret_ps <- function(ps) {
  if (is.na(ps)) return("Unknown")
  
  if (ps >= 0.99) {
    "Very high practical significance"
  } else if (ps >= 0.95) {
    "High practical significance"  
  } else if (ps >= 0.90) {
    "Moderate practical significance"
  } else if (ps >= 0.80) {
    "Weak practical significance"
  } else {
    "No practical significance"
  }
}

#' Print Method for Enhanced p_significance
#'
#' @param x A \code{qbrms_p_significance} object from \code{p_significance()}.
#' @param digits Number of decimal places to display (default 3).
#' @param ... Additional arguments passed to \code{print.data.frame()}.
#' 
#' @return Invisibly returns the input object.
#' 
#' @method print qbrms_p_significance
#' @export
print.qbrms_p_significance <- function(x, digits = 3, ...) {
  cat("Probability of Practical Significance\n")
  cat("=====================================\n\n")
  
  # Show threshold info
  if (length(unique(x$Threshold_low)) == 1 && length(unique(x$Threshold_high)) == 1) {
    cat("Threshold: [", x$Threshold_low[1], ", ", x$Threshold_high[1], "]\n", sep = "")
  } else {
    cat("Thresholds: Parameter-specific (see Threshold columns)\n")
  }
  cat("\n")
  
  # Format for printing
  x_print <- x[, c("Parameter", "ps", "Median", "CI_low", "CI_high", 
                   "Threshold_low", "Threshold_high", "Interpretation")]
  names(x_print) <- c("Parameter", "p(sig)", "Median", "95% CI Low", "95% CI High",
                      "Thresh Low", "Thresh High", "Interpretation")
  
  # Round numeric columns
  numeric_cols <- sapply(x_print, is.numeric)
  x_print[numeric_cols] <- lapply(x_print[numeric_cols], round, digits = digits)
  
  print.data.frame(x_print, row.names = FALSE, ...)
  
  cat("\nInterpretation:\n")
  cat("p(sig) represents the probability that the effect exceeds the threshold\n")
  cat("in the median's direction (i.e., is practically significant).\n\n")
  cat("p(sig) > 0.99: Very high practical significance\n")
  cat("p(sig) > 0.95: High practical significance\n")
  cat("p(sig) > 0.90: Moderate practical significance\n")
  cat("p(sig) > 0.80: Weak practical significance\n")
  cat("p(sig) < 0.80: No practical significance\n")
  
  invisible(x)
}

#' Summary Method for Enhanced p_significance
#'
#' @param object A \code{qbrms_p_significance} object from \code{p_significance()}.
#' @param ... Additional arguments (currently unused).
#' 
#' @return Invisibly returns the input object.
#' 
#' @method summary qbrms_p_significance  
#' @export
summary.qbrms_p_significance <- function(object, ...) {
  cat("Summary of Practical Significance Analysis\n")
  cat("=========================================\n\n")
  
  # Count parameters by interpretation
  interp_table <- table(object$Interpretation)
  cat("Distribution of significance levels:\n")
  for (interp in names(interp_table)) {
    cat("  ", interp, ": ", interp_table[interp], " parameter(s)\n", sep = "")
  }
  
  cat("\n")
  
  # Show parameters with high significance
  high_sig <- object[object$ps >= 0.95 & !is.na(object$ps), ]
  if (nrow(high_sig) > 0) {
    cat("Parameters with high practical significance (ps >= 0.95):\n")
    for (i in 1:nrow(high_sig)) {
      direction <- if (high_sig$Median[i] > 0) "positive" else "negative"
      cat("  - ", high_sig$Parameter[i], 
          " (ps = ", round(high_sig$ps[i], 3), 
          ", ", direction, " effect)\n", sep = "")
    }
  } else {
    cat("No parameters show high practical significance.\n")
  }
  
  invisible(object)
}

#' Plot Method for Enhanced p_significance
#'
#' @description Create a visual plot of probability of practical significance results.
#' 
#' @param x A \code{qbrms_p_significance} object from \code{p_significance()}.
#' @param ... Additional arguments passed to ggplot2 functions.
#' 
#' @return A ggplot2 object.
#' 
#' @method plot qbrms_p_significance
#' @export
plot.qbrms_p_significance <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }
  
  # Prepare data for plotting
  plot_data <- x[!is.na(x$ps), ]
  plot_data$Parameter <- factor(plot_data$Parameter, 
                                levels = plot_data$Parameter[order(plot_data$ps)])
  
  # Add direction information
  plot_data$Direction <- ifelse(plot_data$Median > 0, "Positive", "Negative")
  
  # Create significance level categories
  plot_data$sig_level <- cut(plot_data$ps,
                             breaks = c(0, 0.8, 0.9, 0.95, 0.99, 1),
                             labels = c("None", "Weak", "Moderate", "High", "Very High"),
                             include.lowest = TRUE)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ps, y = Parameter)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = ps, 
                                       y = Parameter, yend = Parameter,
                                       color = sig_level),
                          linewidth = 2) +
    ggplot2::geom_point(ggplot2::aes(color = sig_level, shape = Direction), 
                        size = 4) +
    ggplot2::scale_color_manual(
      values = c("None" = "#E74C3C",
                 "Weak" = "#F39C12", 
                 "Moderate" = "#F1C40F",
                 "High" = "#27AE60",
                 "Very High" = "#16A085"),
      name = "Practical\nSignificance"
    ) +
    ggplot2::scale_shape_manual(
      values = c("Positive" = 16, "Negative" = 17),
      name = "Effect\nDirection"
    ) +
    ggplot2::geom_vline(xintercept = c(0.8, 0.9, 0.95, 0.99), 
                        linetype = "dashed", alpha = 0.3) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = scales::percent_format()
    )
  
  # Create subtitle dynamically based on threshold values
  # Handle both uniform and parameter-specific thresholds
  thresh_low_vals <- unique(plot_data$Threshold_low)
  thresh_high_vals <- unique(plot_data$Threshold_high)
  
  if (length(thresh_low_vals) == 1 && length(thresh_high_vals) == 1) {
    # All parameters have same threshold
    subtitle_text <- paste0("Threshold: [", round(thresh_low_vals[1], 3), 
                            ", ", round(thresh_high_vals[1], 3), "]")
  } else {
    # Parameter-specific thresholds
    subtitle_text <- "Parameter-specific thresholds (see data for details)"
  }
  
  p <- p +
    ggplot2::labs(
      title = "Probability of Practical Significance",
      subtitle = subtitle_text,
      x = "Probability of Practical Significance",
      y = "Parameter"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )
  
  return(p)
}
# =============================================================================
# SECTION 4: INTERNAL HELPERS (namespaced with .ba_)
# =============================================================================

# Simulate draws for fixed effects from summaries in qbrms_fit
.ba_extract_parameter_samples <- function(object, nsim) {
  out <- list()
  if (!is.null(object$fit$summary.fixed)) {
    m <- object$fit$summary.fixed
    nms <- rownames(m)
    mu  <- m[, "mean"]
    sdv <- m[, "sd"]
    for (i in seq_along(nms)) {
      if (is.finite(mu[i]) && is.finite(sdv[i]) && sdv[i] > 0) {
        out[[nms[i]]] <- stats::rnorm(nsim, mu[i], sdv[i])
      }
    }
  } else {
    cf <- tryCatch(stats::coef(object), error = function(e) NULL)
    if (!is.null(cf)) {
      for (i in seq_along(cf)) {
        nm <- names(cf)[i]
        sd <- abs(cf[i]) * 0.1 + 0.01
        out[[nm]] <- stats::rnorm(nsim, cf[i], sd)
      }
    }
  }
  if (!length(out)) stop("Could not extract posterior samples from model")
  out
}

# Very simple hypothesis parser
.ba_parse_hypothesis <- function(h, param_names) {
  typ <- if (grepl("=", h, fixed = TRUE)) "point" else if (grepl(">", h, fixed = TRUE) || grepl("<", h, fixed = TRUE)) "interval" else "comparison"
  params <- param_names[vapply(param_names, function(p) grepl(p, h, fixed = TRUE), logical(1))]
  list(type = typ, parameters = params, expression = h)
}

.ba_compute_point_bf <- function(samples, info, null, rope) {
  p <- info$parameters[1]
  s <- samples[[p]]
  if (is.null(s)) stop("Parameter '", p, "' not found in posterior draws")
  
  if (!is.null(rope) && length(rope) == 2) {
    post <- mean(s >= rope[1] && s <= rope[2])
    prior <- min(0.99, max(0.01, (rope[2] - rope[1]) / 10))
  } else {
    post  <- mean(abs(s - null) < 0.01)
    prior <- 0.01
  }
  bf <- (post / (1 - post)) / (prior / (1 - prior))
  list(bf = max(bf, .Machine$double.eps), posterior_prob = post, prior_prob = prior,
       details = list(param = p, null = null, rope = rope))
}

.ba_compute_interval_bf <- function(samples, info) {
  p <- info$parameters[1]
  s <- samples[[p]]
  if (is.null(s)) stop("Parameter '", p, "' not found in posterior draws")
  
  if (grepl(">", info$expression, fixed = TRUE)) {
    thr <- as.numeric(sub(".*>\\s*", "", info$expression))
    post <- mean(s > thr)
  } else if (grepl("<", info$expression, fixed = TRUE)) {
    thr <- as.numeric(sub(".*<\\s*", "", info$expression))
    post <- mean(s < thr)
  } else {
    thr  <- NA_real_
    post <- 0.5
  }
  prior <- 0.5
  bf <- (post / (1 - post)) / (prior / (1 - prior))
  list(bf = max(bf, .Machine$double.eps), posterior_prob = post, prior_prob = prior,
       details = list(param = p, threshold = thr))
}

.ba_compute_comparison_bf <- function(samples, info) {
  list(bf = 1.0, posterior_prob = 0.5, prior_prob = 0.5,
       details = list(type = "comparison", note = "Not implemented"))
}

.ba_interpret_bf <- function(bf) {
  if (bf >= 100)      "Extreme evidence for H1"
  else if (bf >= 30) "Very strong evidence for H1"
  else if (bf >= 10) "Strong evidence for H1"
  else if (bf >= 3)  "Moderate evidence for H1"
  else if (bf >= 1)  "Weak evidence for H1"
  else if (bf >= 1/3)"Weak evidence for H0"
  else if (bf >= 1/10)"Moderate evidence for H0"
  else if (bf >= 1/30)"Strong evidence for H0"
  else if (bf >= 1/100)"Very strong evidence for H0"
  else               "Extreme evidence for H0"
}

.ba_compute_hdi <- function(samples, prob) {
  s <- sort(samples)
  n <- length(s)
  k <- max(1L, floor(n * prob))
  starts <- 1:(n - k + 1)
  widths <- s[starts + k - 1] - s[starts]
  j <- which.min(widths)
  c(s[j], s[j + k - 1])
}

# Build an EMM grid from specs and data
.ba_create_emm_grid <- function(data, specs, at) {
  specs <- trimws(specs)
  specs <- specs[nzchar(specs)]
  if (!length(specs)) stop("No 'specs' terms provided for EMMs")
  
  # Domains for target terms
  domains <- lapply(specs, function(v) {
    if (!v %in% names(data)) stop("Variable '", v, "' not found in data")
    if (is.factor(data[[v]])) {
      levels(data[[v]])
    } else {
      unique(stats::na.omit(data[[v]]))
    }
  })
  names(domains) <- specs
  
  grid <- do.call(expand.grid, c(domains, stringsAsFactors = FALSE))
  
  # Fill other variables
  for (v in names(data)) {
    if (!(v %in% names(grid))) {
      if (!is.null(at) && v %in% names(at)) {
        grid[[v]] <- at[[v]]
      } else if (is.numeric(data[[v]])) {
        grid[[v]] <- mean(data[[v]], na.rm = TRUE)
      } else if (is.factor(data[[v]])) {
        tab <- table(data[[v]])
        lvl <- names(tab)[which.max(tab)]
        grid[[v]] <- factor(lvl, levels = levels(data[[v]]))
      } else {
        grid[[v]] <- data[[v]][1]
      }
    }
  }
  grid
}

# Draw nsim coefficient vectors from N(coef, V)
.ba_sample_posterior_coefs <- function(coefs, V, nsim) {
  p <- length(coefs)
  if (!is.matrix(V) || nrow(V) != p || ncol(V) != p) {
    V <- diag(rep(0.01, p))
  }
  # Ensure positive-definiteness
  V <- V + diag(1e-8, p)
  L <- try(chol(V), silent = TRUE)
  if (inherits(L, "try-error")) {
    sds <- sqrt(pmax(1e-8, diag(V)))
    mat <- replicate(p, stats::rnorm(nsim, 0, 1))
    sweep(mat, 2, sds, `*`) + matrix(rep(coefs, each = nsim), nrow = nsim)
  } else {
    Z <- matrix(stats::rnorm(p * nsim), nrow = p, ncol = nsim)
    t(matrix(coefs, nrow = p, ncol = nsim) + L %*% Z) # nsim x p
  }
}


##' Print Method for p_significance
#'
#' @description Print results from probability of practical significance analysis.
#' 
#' @param x A \code{qbrms_p_significance} object from \code{p_significance()}.
#' @param digits Number of decimal places to display (default 3).
#' @param ... Additional arguments passed to \code{print.data.frame()}.
#' 
#' @return Invisibly returns the input object.
#' 
#' @method print qbrms_p_significance
#' @export
print.qbrms_p_significance <- function(x, digits = 3, ...) {
  cat("Probability of Practical Significance\n")
  cat("=====================================\n\n")
  cat("Threshold: [", x$Threshold_low[1], ", ", x$Threshold_high[1], "]\n\n", sep = "")
  
  # Format for printing
  x_print <- x[, c("Parameter", "ps", "Median", "CI_low", "CI_high", "Interpretation")]
  names(x_print) <- c("Parameter", "p(sig)", "Median", "95% CI Low", "95% CI High", "Interpretation")
  
  # Round numeric columns
  numeric_cols <- sapply(x_print, is.numeric)
  x_print[numeric_cols] <- lapply(x_print[numeric_cols], round, digits = digits)
  
  print.data.frame(x_print, row.names = FALSE, ...)
  
  cat("\nInterpretation:\n")
  cat("p(sig) > 0.99: Very high practical significance\n")
  cat("p(sig) > 0.95: High practical significance\n")
  cat("p(sig) > 0.90: Moderate practical significance\n")
  cat("p(sig) > 0.80: Weak practical significance\n")
  cat("p(sig) < 0.80: No practical significance\n")
  
  invisible(x)
}

#' Summary Method for p_significance
#'
#' @description Provide a summary of probability of practical significance results.
#' 
#' @param object A \code{qbrms_p_significance} object from \code{p_significance()}.
#' @param ... Additional arguments (currently unused).
#' 
#' @return Invisibly returns the input object.
#' 
#' @method summary qbrms_p_significance
#' @export
summary.qbrms_p_significance <- function(object, ...) {
  cat("Summary of Practical Significance\n")
  cat("=================================\n\n")
  
  # Count parameters by interpretation
  interp_table <- table(object$Interpretation)
  
  for (interp in names(interp_table)) {
    cat(interp, ": ", interp_table[interp], " parameter(s)\n", sep = "")
  }
  
  cat("\n")
  
  # Show parameters with high significance
  high_sig <- object[object$ps >= 0.95 & !is.na(object$ps), ]
  if (nrow(high_sig) > 0) {
    cat("Parameters with high practical significance (ps >= 0.95):\n")
    for (i in 1:nrow(high_sig)) {
      cat("  - ", high_sig$Parameter[i], 
          " (ps = ", round(high_sig$ps[i], 3), ")\n", sep = "")
    }
  }
  
  invisible(object)
}

#' Plot Method for p_significance
#'
#' @description Create a visual plot of probability of practical significance results.
#' 
#' @param x A \code{qbrms_p_significance} object from \code{p_significance()}.
#' @param ... Additional arguments passed to ggplot2 functions.
#' 
#' @return A ggplot2 object.
#' 
#' @method plot qbrms_p_significance
#' @export
plot.qbrms_p_significance <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }
  
  # Prepare data for plotting
  plot_data <- x[!is.na(x$ps), ]
  plot_data$Parameter <- factor(plot_data$Parameter, 
                                levels = plot_data$Parameter[order(plot_data$ps)])
  
  # Create color scale based on ps value
  plot_data$sig_level <- cut(plot_data$ps,
                             breaks = c(0, 0.8, 0.9, 0.95, 0.99, 1),
                             labels = c("None", "Weak", "Moderate", "High", "Very High"),
                             include.lowest = TRUE)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ps, y = Parameter)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = ps, 
                                       y = Parameter, yend = Parameter,
                                       color = sig_level),
                          linewidth = 2) +
    ggplot2::geom_point(ggplot2::aes(color = sig_level), size = 4) +
    ggplot2::scale_color_manual(
      values = c("None" = "#E74C3C",
                 "Weak" = "#F39C12", 
                 "Moderate" = "#F1C40F",
                 "High" = "#27AE60",
                 "Very High" = "#16A085"),
      name = "Practical\nSignificance"
    ) +
    ggplot2::geom_vline(xintercept = c(0.8, 0.9, 0.95, 0.99), 
                        linetype = "dashed", alpha = 0.3) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = scales::percent_format()
    ) +
    ggplot2::labs(
      title = "Probability of Practical Significance",
      subtitle = paste0("Threshold: [", round(x$Threshold_low[1], 3), 
                        ", ", round(x$Threshold_high[1], 3), "]"),
      x = "Probability of Practical Significance",
      y = "Parameter"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )
  
  return(p)
}
#' Compare Significance Across Multiple Models
#'
#' @description
#' Compare the probability of practical significance for parameters
#' across multiple qbrms models.
#'
#' @param ... qbrms_fit objects to compare
#' @param parameters Character vector of parameters to compare
#' @param threshold Threshold specification (same as p_significance)
#' @param model_names Character vector of model names
#'
#' @return Data frame with comparison results
#' @export
compare_significance <- function(..., 
                                 parameters = NULL,
                                 threshold = "default",  # Changed from rope
                                 model_names = NULL) {
  
  models <- list(...)
  
  if (is.null(model_names)) {
    model_names <- paste0("Model_", seq_along(models))
  }
  
  # Check all are qbrms_fit objects
  if (!all(sapply(models, inherits, "qbrms_fit"))) {
    stop("All objects must be qbrms_fit objects")
  }
  
  # Get p_significance for each model
  # Use match.fun to explicitly resolve the function for R CMD check
  p_sig_fun <- match.fun("p_significance")
  ps_list <- lapply(models, function(m) {
    p_sig_fun(m, parameters = parameters, threshold = threshold, verbose = FALSE)
  })
  
  # Combine results
  combined_results <- NULL
  
  for (i in seq_along(ps_list)) {
    ps_df <- ps_list[[i]]
    ps_df$Model <- model_names[i]
    
    if (is.null(combined_results)) {
      combined_results <- ps_df
    } else {
      combined_results <- rbind(combined_results, ps_df)
    }
  }
  
  # Reshape for comparison
  comparison <- reshape(combined_results[, c("Parameter", "ps", "Model")],
                        idvar = "Parameter",
                        timevar = "Model",
                        direction = "wide")
  
  names(comparison) <- gsub("ps\\.", "", names(comparison))
  
  class(comparison) <- c("qbrms_significance_comparison", "data.frame")
  return(comparison)
}
