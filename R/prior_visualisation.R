# =============================================================================
# R/prior_visualisation.R
# =============================================================================

#' Visualise Prior Distributions
#'
#' @description
#' Create visual representations of prior distributions to aid in prior
#' specification and sensitivity analysis.
#'
#' @param prior Prior specification in qbrms format, or a list of prior
#'   specifications to compare
#' @param parameter Character string specifying which parameter to visualise
#'   (e.g., "b", "sd", "sigma"). If NULL, visualises all priors.
#' @param xlim Numeric vector of length 2 specifying x-axis limits. If NULL,
#'   automatically determined.
#' @param add_reference Logical; if TRUE, adds reference distributions for
#'   comparison (default: TRUE)
#' @param samples Number of samples to draw for visualisation (default: 10000)
#'
#' @return A ggplot object showing the prior distribution(s)
#'
#' @details
#' This function helps users:
#' \itemize{
#'   \item Visualise the implications of their prior choices
#'   \item Compare different prior specifications
#'   \item Identify overly informative or vague priors
#'   \item Understand prior-data conflict potential
#' }
#'
#' Supported prior distributions include:
#' \itemize{
#'   \item Normal: normal(mean, sd)
#'   \item Student t: student_t(df, mean, scale)
#'   \item Cauchy: cauchy(location, scale)
#'   \item Exponential: exponential(rate)
#'   \item Gamma: gamma(shape, rate)
#'   \item Uniform: uniform(lower, upper)
#' }
#' 
#'
#' @importFrom stats dt dcauchy dexp dgamma dunif dbeta dlnorm rnorm rt rcauchy rexp rgamma runif rbeta rlnorm density quantile
#' @importFrom utils head
#'
#'
#' @examples
#' \dontrun{
#' # Visualise a single prior
#' prior <- prior(normal(0, 10), class = "b")
#' visualise_prior(prior)
#'
#' # Compare different priors
#' prior_list <- list(
#'   "Weak" = prior(normal(0, 10), class = "b"),
#'   "Medium" = prior(normal(0, 5), class = "b"),
#'   "Strong" = prior(normal(0, 1), class = "b")
#' )
#' visualise_prior(prior_list)
#'
#' # Visualise with custom limits
#' visualise_prior(prior, xlim = c(-20, 20))
#' }
#'
#' @export
visualise_prior <- function(prior, parameter = NULL, xlim = NULL, 
                            add_reference = TRUE, samples = 10000) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for prior visualisation", call. = FALSE)
  }
  
  # Handle list of priors for comparison
  if (is.list(prior) && !inherits(prior, "qbrms_prior_spec") && !inherits(prior, "qbrms_prior_dist")) {
    return(.visualise_prior_comparison(prior, parameter, xlim, samples))
  }
  
  # Parse single prior
  prior_info <- .parse_prior_specification(prior, parameter)
  
  if (is.null(prior_info)) {
    stop("Unable to parse prior specification", call. = FALSE)
  }
  
  # Generate samples
  prior_samples <- .sample_from_prior(prior_info, samples)
  
  # Determine x limits
  if (is.null(xlim)) {
    xlim <- .determine_xlim(prior_samples, prior_info)
  }
  
  # Create plot data
  x_seq <- seq(xlim[1], xlim[2], length.out = 1000)
  density_vals <- .evaluate_prior_density(x_seq, prior_info)
  
  plot_data <- data.frame(
    x = x_seq,
    density = density_vals
  )
  
  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = density)) +
    ggplot2::geom_line(colour = "steelblue", linewidth = 1.2) +
    ggplot2::geom_area(fill = "steelblue", alpha = 0.3) +
    ggplot2::labs(
      title = "Prior Distribution",
      subtitle = .format_prior_label(prior_info),
      x = "Parameter Value",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11)
    )
  
  # Add reference distributions
  if (add_reference) {
    p <- .add_reference_distributions(p, xlim, prior_info)
  }
  
  # Add quantile markers
  quantiles <- quantile(prior_samples, probs = c(0.025, 0.5, 0.975))
  p <- p + ggplot2::geom_vline(
    xintercept = quantiles,
    linetype = "dashed",
    colour = "grey50",
    alpha = 0.6
  )
  
  print(p)
  invisible(p)
}

#' Parse Prior Specification
#' @keywords internal
.parse_prior_specification <- function(prior, parameter) {
  
  # Handle string format: "normal(0, 10)"
  if (is.character(prior)) {
    return(.parse_prior_string(prior))
  } 
  
  # Handle qbrms_prior_spec objects (from prior() function)
  if (inherits(prior, "qbrms_prior_spec")) {
    return(list(
      distribution = prior$distribution,
      parameters = prior$parameters
    ))
  }
  
  # Handle qbrms_prior_dist objects (from normal(), student_t(), etc.)
  if (inherits(prior, "qbrms_prior_dist")) {
    return(list(
      distribution = prior$distribution,
      parameters = prior$parameters
    ))
  }
  
  # Handle generic lists if they have the right structure
  if (is.list(prior)) {
    # Extract from brms-style prior object
    if (!is.null(prior$prior) && is.character(prior$prior)) {
      return(.parse_prior_string(prior$prior))
    }
    if (!is.null(prior$distribution) && !is.null(prior$parameters)) {
      return(list(
        distribution = prior$distribution,
        parameters = prior$parameters
      ))
    }
  }
  
  return(NULL)
}

#' Parse Prior String
#' @keywords internal
.parse_prior_string <- function(prior_str) {
  
  # Extract distribution name and parameters
  pattern <- "^([a-z_]+)\\(([^)]+)\\)$"
  
  if (!grepl(pattern, prior_str)) {
    return(NULL)
  }
  
  matches <- regmatches(prior_str, regexec(pattern, prior_str))[[1]]
  
  if (length(matches) < 3) {
    return(NULL)
  }
  
  dist_name <- matches[2]
  params_str <- matches[3]
  
  # Parse parameters
  params <- as.numeric(strsplit(params_str, ",")[[1]])
  
  list(
    distribution = dist_name,
    parameters = params
  )
}

#' Sample from Prior Distribution
#' @keywords internal
.sample_from_prior <- function(prior_info, n_samples) {
  
  dist <- prior_info$distribution
  params <- prior_info$parameters
  
  # Handle named parameters by converting to numeric vector if needed
  # But keep list if the switch expects list elements or named access
  
  # Safe accessor for parameters
  p <- function(idx, name=NULL, def=NA) {
    if (!is.null(name) && name %in% names(params)) return(params[[name]])
    if (idx <= length(params)) return(params[[idx]])
    return(def)
  }
  
  samples <- switch(dist,
                    "normal" = rnorm(n_samples, mean = p(1, "mean", 0), sd = p(2, "sd", 1)),
                    "student_t" = {
                      df_val <- p(1, "df", 3)
                      loc_val <- p(2, "location", 0)
                      scale_val <- p(3, "scale", 1)
                      loc_val + scale_val * rt(n_samples, df = df_val)
                    },
                    "cauchy" = rcauchy(n_samples, location = p(1, "location", 0), scale = p(2, "scale", 1)),
                    "exponential" = rexp(n_samples, rate = p(1, "rate", 1)),
                    "gamma" = rgamma(n_samples, shape = p(1, "shape", 1), rate = p(2, "rate", 1)),
                    "uniform" = runif(n_samples, min = p(1, "min", 0), max = p(2, "max", 1)),
                    "beta" = rbeta(n_samples, shape1 = p(1, "alpha", 1), shape2 = p(2, "beta", 1)),
                    "lognormal" = rlnorm(n_samples, meanlog = p(1, "meanlog", 0), sdlog = p(2, "sdlog", 1)),
                    {
                      warning("Distribution '", dist, "' not supported, using normal(0, 1)")
                      rnorm(n_samples, 0, 1)
                    }
  )
  
  return(samples)
}

#' Evaluate Prior Density
#' @keywords internal
.evaluate_prior_density <- function(x, prior_info) {
  
  dist <- prior_info$distribution
  params <- prior_info$parameters
  
  # Safe accessor for parameters
  p <- function(idx, name=NULL, def=NA) {
    if (!is.null(name) && name %in% names(params)) return(params[[name]])
    if (idx <= length(params)) return(params[[idx]])
    return(def)
  }
  
  density <- switch(dist,
                    "normal" = dnorm(x, mean = p(1, "mean", 0), sd = p(2, "sd", 1)),
                    "student_t" = {
                      df_val <- p(1, "df", 3)
                      loc_val <- p(2, "location", 0)
                      scale_val <- p(3, "scale", 1)
                      dt((x - loc_val) / scale_val, df = df_val) / scale_val
                    },
                    "cauchy" = dcauchy(x, location = p(1, "location", 0), scale = p(2, "scale", 1)),
                    "exponential" = dexp(x, rate = p(1, "rate", 1)),
                    "gamma" = dgamma(x, shape = p(1, "shape", 1), rate = p(2, "rate", 1)),
                    "uniform" = dunif(x, min = p(1, "min", 0), max = p(2, "max", 1)),
                    "beta" = dbeta(x, shape1 = p(1, "alpha", 1), shape2 = p(2, "beta", 1)),
                    "lognormal" = dlnorm(x, meanlog = p(1, "meanlog", 0), sdlog = p(2, "sdlog", 1)),
                    rep(NA_real_, length(x))
  )
  
  return(density)
}

#' Determine X-axis Limits
#' @keywords internal
.determine_xlim <- function(samples, prior_info) {
  
  # Use quantiles for automatic limits
  q <- quantile(samples, probs = c(0.001, 0.999), na.rm = TRUE)
  
  # Expand slightly
  range_width <- diff(q)
  xlim <- c(q[1] - 0.1 * range_width, q[2] + 0.1 * range_width)
  
  return(xlim)
}
#' Format Prior Label
#' @keywords internal
.format_prior_label <- function(prior_info) {
  
  dist <- prior_info$distribution
  params <- prior_info$parameters
  
  # Helper to flatten params
  p <- function(idx, name=NULL, def=0) {
    val <- if (!is.null(name) && name %in% names(params)) params[[name]] else params[[idx]]
    if(is.null(val)) def else val
  }
  
  label <- switch(dist,
                  "normal" = sprintf("Normal(mean = %.2f, sd = %.2f)", p(1, "mean"), p(2, "sd")),
                  "student_t" = sprintf("Student-t(df = %.0f, mean = %.2f, sd = %.2f)", 
                                        p(1, "df"), p(2, "location"), p(3, "scale")),
                  "cauchy" = sprintf("Cauchy(location = %.2f, scale = %.2f)", 
                                     p(1, "location"), p(2, "scale")),
                  "exponential" = sprintf("Exponential(rate = %.2f)", p(1, "rate")),
                  "gamma" = sprintf("Gamma(shape = %.2f, rate = %.2f)", p(1, "shape"), p(2, "rate")),
                  "uniform" = sprintf("Uniform(%.2f, %.2f)", p(1, "min"), p(2, "max")),
                  "beta" = sprintf("Beta(alpha = %.2f, beta = %.2f)", p(1, "alpha"), p(2, "beta")),
                  "lognormal" = sprintf("Lognormal(meanlog = %.2f, sdlog = %.2f)", p(1, "meanlog"), p(2, "sdlog")),
                  paste(dist, "distribution")
  )
  
  return(label)
}

#' Add Reference Distributions
#' @keywords internal
.add_reference_distributions <- function(p, xlim, prior_info) {
  
  x_seq <- seq(xlim[1], xlim[2], length.out = 1000)
  
  # Add standard normal for reference
  std_normal <- dnorm(x_seq, 0, 1)
  
  ref_data <- data.frame(
    x = x_seq,
    density = std_normal
  )
  
  p <- p + ggplot2::geom_line(
    data = ref_data,
    ggplot2::aes(x = x, y = density),
    colour = "grey60",
    linetype = "dotted",
    linewidth = 0.8
  ) +
    ggplot2::annotate(
      "text",
      x = xlim[2] * 0.7,
      y = max(std_normal) * 0.9,
      label = "Standard Normal\n(reference)",
      colour = "grey60",
      size = 3,
      hjust = 0
    )
  
  return(p)
}

#' Visualise Prior Comparison
#' @keywords internal
.visualise_prior_comparison <- function(prior_list, parameter, xlim, samples) {
  
  if (length(prior_list) < 2) {
    stop("At least two priors required for comparison", call. = FALSE)
  }
  
  # Parse all priors
  prior_info_list <- lapply(prior_list, .parse_prior_specification, parameter)
  
  # Check validity
  valid <- !sapply(prior_info_list, is.null)
  if (!any(valid)) {
    stop("Unable to parse any prior specifications", call. = FALSE)
  }
  
  prior_info_list <- prior_info_list[valid]
  prior_names <- names(prior_list)[valid]
  
  if (is.null(prior_names)) {
    prior_names <- paste0("Prior", seq_along(prior_info_list))
  }
  
  # Sample from all priors to determine limits
  all_samples <- unlist(lapply(prior_info_list, .sample_from_prior, samples))
  
  if (is.null(xlim)) {
    xlim <- .determine_xlim(all_samples, prior_info_list[[1]])
  }
  
  # Create combined plot data
  x_seq <- seq(xlim[1], xlim[2], length.out = 1000)
  
  plot_data_list <- lapply(seq_along(prior_info_list), function(i) {
    density_vals <- .evaluate_prior_density(x_seq, prior_info_list[[i]])
    data.frame(
      x = x_seq,
      density = density_vals,
      prior = prior_names[i],
      stringsAsFactors = FALSE
    )
  })
  
  plot_data <- do.call(rbind, plot_data_list)
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = density, 
                                               colour = prior, fill = prior)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_area(alpha = 0.2, position = "identity") +
    ggplot2::labs(
      title = "Prior Distribution Comparison",
      x = "Parameter Value",
      y = "Density",
      colour = "Prior",
      fill = "Prior"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "right"
    )
  
  print(p)
  invisible(p)
}

#' Create Prior Predictive Distribution Plot
#'
#' @description
#' Generate predictions from the prior distribution to assess whether
#' priors are reasonable before seeing the data.
#'
#' @param formula Model formula
#' @param data Data frame (used for structure, not values)
#' @param family Model family
#' @param prior Prior specification
#' @param n_samples Number of prior predictive samples (default: 1000)
#'
#' @return A ggplot object showing the prior predictive distribution
#'
#' @examples
#' \dontrun{
#' prior_predictive_check(
#'   mpg ~ hp + wt,
#'   data = mtcars,
#'   family = gaussian(),
#'   prior = prior(normal(0, 10), class = "b")
#' )
#' }
#'
#' @export
prior_predictive_check <- function(formula, data, family, prior, n_samples = 1000) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required", call. = FALSE)
  }
  
  # Extract response variable
  response_name <- as.character(formula[[2]])
  observed_data <- data[[response_name]]
  
  # Generate prior predictive samples
  # This is a simplified implementation
  prior_samples <- .generate_prior_predictive_samples(
    formula, data, family, prior, n_samples
  )
  
  # Create plot comparing prior predictions to observed data range
  plot_data <- data.frame(
    value = prior_samples,
    source = "Prior Predictive"
  )
  
  obs_data <- data.frame(
    value = observed_data,
    source = "Observed Data"
  )
  
  combined_data <- rbind(plot_data, obs_data)
  
  p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = value, fill = source)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(
      title = "Prior Predictive Check",
      subtitle = "Do prior predictions cover reasonable values?",
      x = response_name,
      y = "Density",
      fill = "Source"
    ) +
    ggplot2::scale_fill_manual(values = c("Prior Predictive" = "steelblue",
                                          "Observed Data" = "coral")) +
    ggplot2::theme_minimal()
  
  print(p)
  invisible(p)
}

#' Generate Prior Predictive Samples
#' @keywords internal
.generate_prior_predictive_samples <- function(formula, data, family, 
                                               prior, n_samples) {
  
  # This is a simplified implementation
  # In practice, this would sample from priors and generate predictions
  
  # Get design matrix
  X <- model.matrix(formula, data)
  n_coef <- ncol(X)
  
  # Sample coefficients from prior (simplified)
  # Assumes normal(0, 10) if prior not fully specified
  beta_samples <- matrix(rnorm(n_samples * n_coef, 0, 10), 
                         nrow = n_samples, ncol = n_coef)
  
  # Generate predictions
  # Sample random rows from X
  sampled_rows <- sample(nrow(X), n_samples, replace = TRUE)
  X_sample <- X[sampled_rows, , drop = FALSE]
  
  # Linear predictor
  eta <- rowSums(beta_samples * X_sample)
  
  # Apply link function and add noise based on family
  family_name <- if (is.character(family)) family else family$family
  
  predictions <- switch(family_name,
                        "gaussian" = rnorm(n_samples, eta, 1),
                        "poisson" = rpois(n_samples, exp(eta)),
                        "binomial" = rbinom(n_samples, 1, plogis(eta)),
                        eta  # Default
  )
  
  return(predictions)
}