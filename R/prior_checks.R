# =============================================================================
# R/prior_checks.R - Prior checks and density plotting (qbrms)
# =============================================================================

#' @title Prior Predictive Checks and Density Plotting
#' @description Functions for prior checks and density plotting with ggplot2
#'   extensibility.
#' @details This file provides functionality for prior predictive checks and
#'   density plotting that integrates with ggplot2 for full customisation.
#' @keywords internal
#' @name prior_checks
NULL

# Fallback for %||% if not provided elsewhere ---------------------------------
#' @keywords internal

# =============================================================================
# SECTION 1: MAIN PRIOR-ONLY FUNCTIONS
# =============================================================================

#' @title Prior Predictive Checks Without Data
#'
#' @description
#' Generate prior predictive samples for a model defined by \code{formula},
#' without requiring observed data, and return a ggplot object to visualise
#' the implied distribution.
#'
#' @param formula Model formula.
#' @param family Model family (default \code{gaussian()}).
#' @param prior Prior specifications (default \code{NULL}).
#' @param n_obs Number of observations to simulate (default \code{100}).
#' @param ndraws Number of prior draws (default \code{100}).
#' @param type Plot type, one of \code{"dens_overlay"} or \code{"hist"}.
#' @param seed Optional random seed.
#' @param predictor_values Named list of fixed predictor values (default \code{NULL}).
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#'
#' @return A \code{ggplot2} object.
#' @export
pp_check_prior <- function(formula,
                           family = gaussian(),
                           prior = NULL,
                           n_obs = 100,
                           ndraws = 100,
                           type = "dens_overlay",
                           seed = NULL,
                           predictor_values = NULL,
                           verbose = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  if (is.character(formula)) formula <- as.formula(formula)
  
  if (verbose) cat("Generating prior predictive samples (no observed data required)...\n")
  
  synthetic_data <- .create_synthetic_data(formula, n_obs, predictor_values, verbose)
  
  X <- tryCatch(
    stats::model.matrix(formula, data = synthetic_data),
    error = function(e) {
      if (verbose) cat("Warning: Could not create design matrix; using simple structure\n")
      cbind(`(Intercept)` = 1, x = stats::rnorm(n_obs))
    }
  )
  
  coef_names <- colnames(X)
  if (verbose) cat("Design matrix created:", n_obs, "obs,", ncol(X), "coefficients\n")
  
  prior_specs <- .extract_prior_specs_standalone(prior, coef_names, verbose)
  yrep <- .generate_prior_only_samples(X, prior_specs, family, ndraws, verbose)
  
  if (verbose) {
    rng <- range(yrep)
    cat("Prior samples generated:\n")
    cat("  Range:", paste(round(rng, 2), collapse = " to "), "\n")
    cat("  Mean:", round(mean(yrep), 2), "\n")
  }
  
  .create_prior_only_plot(yrep, type, ndraws, verbose)
}

#' Create Prior-Only Object for pp_check
#'
#' Construct a small \code{qbrms_prior_only} object that contains simulated data
#' and prior draws, suitable for passing to \code{pp_check()}.
#'
#' @param formula Model formula.
#' @param family Model family (default \code{gaussian()}).
#' @param prior Prior specifications (default \code{NULL} uses defaults).
#' @param n_obs Number of observations to simulate (default \code{100}).
#' @param predictor_values Named list of fixed predictor values (default \code{NULL}).
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#'
#' @return An object of class \code{qbrms_prior_only}.
#' @export
create_prior_object <- function(formula,
                                family = gaussian(),
                                prior = NULL,
                                n_obs = 100,
                                predictor_values = NULL,
                                verbose = TRUE) {
  
  if (verbose) cat("Creating prior-only object for pp_check...\n")
  if (is.character(formula)) formula <- as.formula(formula)
  
  synthetic_data <- .create_synthetic_data(formula, n_obs, predictor_values, verbose)
  X <- stats::model.matrix(formula, data = synthetic_data)
  coef_names <- colnames(X)
  prior_specs <- .extract_prior_specs_standalone(prior, coef_names, verbose)
  prior_samples <- .generate_prior_only_samples(X, prior_specs, family, 100, verbose)
  
  result <- list(
    prior_samples = prior_samples,
    original_formula = formula,
    data = synthetic_data,
    family = family,
    model_type = "prior_only",
    prior_specs = prior,
    is_prior_only = TRUE
  )
  
  class(result) <- c("qbrms_prior_only", "qbrms_prior_fit", "qbrms_fit", "list")
  if (verbose) cat("Prior-only object created successfully.\n")
  result
}

# =============================================================================
# SECTION 2: ENHANCED DENSITY PLOTTING
# =============================================================================

#' Density Plot for qbrms Models
#'
#' Create density plots of posterior distributions with optional prior and
#' observed-data overlays. Returns a \code{ggplot2} object that can be modified
#' with standard ggplot2 syntax.
#'
#' @param object A \code{qbrms_fit} object.
#' @param parameter Parameter name to plot. If \code{NULL}, plots the response distribution.
#' @param show_prior Logical; if \code{TRUE}, overlay the prior density.
#' @param show_data Logical; if \code{TRUE}, overlay the observed data density.
#' @param ndraws Number of posterior draws to use (default \code{100}).
#' @param prior_ndraws Number of prior draws to use (default \code{100}).
#' @param alpha_levels Named list controlling transparency for layers.
#' @param colours Named list of colours for layers.
#' @param seed Optional random seed.
#' @param verbose Logical; print progress messages.
#'
#' @return A \code{ggplot2} object.
#' @export
density_plot <- function(object,
                         parameter = NULL,
                         show_prior = FALSE,
                         show_data = FALSE,
                         ndraws = 100,
                         prior_ndraws = 100,
                         alpha_levels = list(posterior = 0.8, prior = 0.6, data = 1.0),
                         colours = list(posterior = "#1F78B4", prior = "#E31A1C", data = "#000000"),
                         seed = NULL,
                         verbose = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  if (!inherits(object, "qbrms_fit")) {
    stop("object must be a qbrms_fit object")
  }
  
  if (is.null(parameter)) {
    plot_type <- "response"
    if (verbose) cat("Creating posterior predictive density plot...\n")
  } else {
    plot_type <- "parameter"
    if (verbose) cat("Creating density plot for parameter:", parameter, "\n")
  }
  
  if (plot_type == "response") {
    plot_data <- .extract_response_densities(object, show_prior, show_data,
                                             ndraws, prior_ndraws, verbose)
  } else {
    plot_data <- .extract_parameter_densities(object, parameter, show_prior,
                                              ndraws, prior_ndraws, verbose)
  }
  
  p <- .create_density_base_plot(plot_data, colours, alpha_levels, parameter, verbose)
  
  class(p) <- c("qbrms_density_plot", class(p))
  attr(p, "qbrms_plot_info") <- list(
    parameter = parameter,
    plot_type = plot_type,
    show_prior = show_prior,
    show_data = show_data
  )
  
  p
}

# Former duplicate of plot_parameters(): private wrapper to forward old defaults.
# Not exported and no Rd topic.

#' @keywords internal
#' @noRd
plot_parameters_prior_checks <- function(object,
                                         pars = NULL,
                                         show_prior = FALSE,
                                         ndraws = 200,
                                         prior_ndraws = 200,
                                         ncol = 2,
                                         alpha_levels = list(posterior = 0.8, prior = 0.5),
                                         colours = list(posterior = "#1F78B4", prior = "#E31A1C"),
                                         verbose = TRUE,
                                         ...) {
  alpha_vec  <- if (is.list(alpha_levels)) as.numeric(unlist(alpha_levels, use.names = FALSE)) else alpha_levels
  colour_vec <- if (is.list(colours))      as.character(unlist(colours,      use.names = FALSE)) else colours
  
  # Resolve our own exported function at runtime from this package's namespace
  fun <- utils::getFromNamespace("plot_parameters", "qbrms")
  
  fun(object       = object,
      pars         = pars,
      show_prior   = show_prior,
      ndraws       = ndraws,
      prior_ndraws = prior_ndraws,
      ncol         = ncol,
      alpha_levels = alpha_vec,
      colours      = colour_vec,
      verbose      = verbose,
      ...)
}


#' Quick Density Comparison
#' @param object A \code{qbrms_fit} object.
#' @param parameter Optional parameter name to focus the comparison.
#' @param ... Additional arguments forwarded to \code{density_plot()}.
#' @return A \code{ggplot2} object.
#' @export
quick_density_comparison <- function(object, parameter = NULL, ...) {
  density_plot(object,
               parameter   = parameter,
               show_prior  = TRUE,
               show_data   = TRUE,
               ...)
}

# =============================================================================
# SECTION 3: S3 METHODS FOR PRIOR-ONLY OBJECTS
# =============================================================================

#' @method pp_check qbrms_prior_only
#' @export
pp_check.qbrms_prior_only <- function(object,
                                      type = "dens_overlay",
                                      ndraws = 100,
                                      seed = NULL,
                                      show_observed = FALSE,
                                      ...) {
  .pp_check_core(object,
                 type          = type,
                 ndraws        = ndraws,
                 seed          = seed,
                 is_prior      = TRUE,
                 show_observed = FALSE,
                 ...)
}

#' @method print qbrms_prior_only
#' @export
print.qbrms_prior_only <- function(x, ...) {
  cat("qbrms Prior-Only Object\n\n")
  cat("Formula:", deparse(x$original_formula), "\n")
  cat("Family: ", extract_family_name(x$family), "\n")
  cat("Data:   ", nrow(x$data), "simulated observations\n")
  
  if (!is.null(x$prior_specs)) {
    cat("Custom priors: Yes\n")
  } else {
    cat("Custom priors: No (using defaults)\n")
  }
  
  cat("\nUse pp_check(object) to visualise the prior predictive distribution.\n")
  invisible(x)
}

# =============================================================================
# SECTION 4: INTERNAL HELPER FUNCTIONS
# =============================================================================

#' Create Synthetic Data for Prior Checks
#' @keywords internal
.create_synthetic_data <- function(formula, n_obs, predictor_values = NULL, verbose = TRUE) {
  all_vars <- all.vars(formula)
  response_var <- all_vars[1]
  predictor_vars <- all_vars[-1]
  
  synthetic_data <- data.frame(row.names = seq_len(n_obs))
  synthetic_data[[response_var]] <- stats::rnorm(n_obs, 0, 1)
  
  for (var in predictor_vars) {
    if (!is.null(predictor_values) && var %in% names(predictor_values)) {
      value <- predictor_values[[var]]
      if (is.numeric(value) && length(value) == 1) {
        synthetic_data[[var]] <- rep(value, n_obs)
        if (verbose) cat("  ", var, ": fixed at", value, "\n")
      } else if (is.character(value) && length(value) == 1) {
        synthetic_data[[var]] <- factor(rep(value, n_obs))
        if (verbose) cat("  ", var, ": fixed at", value, "(factor)\n")
      } else {
        synthetic_data[[var]] <- stats::rnorm(n_obs, 0, 1)
      }
    } else {
      if (grepl("^(age|year|time|count|income|score)", var, ignore.case = TRUE)) {
        synthetic_data[[var]] <- stats::rnorm(n_obs, 0, 1)
      } else if (grepl("^(group|treatment|category|type|gender)", var, ignore.case = TRUE)) {
        synthetic_data[[var]] <- factor(sample(c("A", "B"), n_obs, replace = TRUE))
      } else {
        synthetic_data[[var]] <- stats::rnorm(n_obs, 0, 1)
      }
    }
  }
  
  synthetic_data
}

#' Extract Prior Specifications (Standalone Version)
#' @keywords internal
.extract_prior_specs_standalone <- function(prior, coef_names, verbose = TRUE) {
  prior_specs <- list()
  
  for (coef_name in coef_names) {
    if (coef_name == "(Intercept)") {
      prior_specs[[coef_name]] <- list(distribution = "normal", parameters = list(mean = 0, sd = 2.5))
    } else {
      prior_specs[[coef_name]] <- list(distribution = "normal", parameters = list(mean = 0, sd = 1))
    }
  }
  
  if (!is.null(prior)) {
    if (inherits(prior, "qbrms_prior_spec")) {
      prior_specs <- .apply_single_prior_spec_standalone(prior, prior_specs, coef_names)
    } else if (inherits(prior, "qbrms_prior_list")) {
      for (p in prior) {
        prior_specs <- .apply_single_prior_spec_standalone(p, prior_specs, coef_names)
      }
    }
  }
  
  if (verbose) {
    cat("Prior specifications:\n")
    for (name in names(prior_specs)) {
      spec <- prior_specs[[name]]
      params_str <- paste(names(spec$parameters), "=", spec$parameters, collapse = ", ")
      cat(sprintf("  %s: %s(%s)\n", name, spec$distribution, params_str))
    }
  }
  
  prior_specs
}

#' Apply Single Prior Specification
#' @keywords internal
.apply_single_prior_spec_standalone <- function(prior_spec, prior_specs, coef_names) {
  if (prior_spec$class == "Intercept") {
    if ("(Intercept)" %in% coef_names) {
      prior_specs[["(Intercept)"]] <- list(
        distribution = prior_spec$distribution,
        parameters   = prior_spec$parameters
      )
    }
  } else if (prior_spec$class == "b") {
    if (!is.null(prior_spec$coef)) {
      target_coef <- prior_spec$coef
      if (target_coef %in% coef_names) {
        prior_specs[[target_coef]] <- list(
          distribution = prior_spec$distribution,
          parameters   = prior_spec$parameters
        )
      }
    } else {
      slope_coefs <- coef_names[coef_names != "(Intercept)"]
      for (coef in slope_coefs) {
        prior_specs[[coef]] <- list(
          distribution = prior_spec$distribution,
          parameters   = prior_spec$parameters
        )
      }
    }
  }
  
  prior_specs
}

#' Generate Prior-Only Samples
#' @keywords internal
.generate_prior_only_samples <- function(X, prior_specs, family, ndraws, verbose = TRUE) {
  n_obs <- nrow(X)
  n_coefs <- ncol(X)
  coef_names <- colnames(X)
  
  yrep <- matrix(NA_real_, nrow = ndraws, ncol = n_obs)
  family_name <- extract_family_name(family)
  
  for (i in seq_len(ndraws)) {
    beta_sample <- numeric(n_coefs)
    names(beta_sample) <- coef_names
    
    for (j in seq_len(n_coefs)) {
      coef_name <- coef_names[j]
      prior_spec <- prior_specs[[coef_name]]
      beta_sample[j] <- .sample_from_prior_safe(prior_spec)
    }
    
    linear_pred <- as.numeric(X %*% beta_sample)
    yrep[i, ] <- .generate_response_from_family(linear_pred, family_name, n_obs)
  }
  
  yrep
}

#' Safe Prior Sampling
#' @keywords internal  
.sample_from_prior_safe <- function(prior_spec) {
  switch(prior_spec$distribution,
         # Existing distributions (unchanged)
         "normal" = stats::rnorm(1, 
                                 prior_spec$parameters$mean %||% 0, 
                                 prior_spec$parameters$sd %||% 1),
         "student_t" = {
           df <- prior_spec$parameters$df %||% 3
           location <- prior_spec$parameters$location %||% 0
           scale <- prior_spec$parameters$scale %||% 1
           location + scale * stats::rt(1, df)
         },
         "cauchy" = {
           location <- prior_spec$parameters$location %||% 0
           scale <- prior_spec$parameters$scale %||% 1
           stats::rcauchy(1, location, scale)
         },
         "uniform" = {
           min_val <- prior_spec$parameters$min %||% -10
           max_val <- prior_spec$parameters$max %||% 10
           stats::runif(1, min_val, max_val)
         },
         
         "gamma" = {
           shape <- prior_spec$parameters$shape %||% 2
           rate <- prior_spec$parameters$rate %||% 1
           stats::rgamma(1, shape = shape, rate = rate)
         },
         "beta" = {
           alpha <- prior_spec$parameters$alpha %||% 1
           beta <- prior_spec$parameters$beta %||% 1
           stats::rbeta(1, shape1 = alpha, shape2 = beta)
         },
         "exponential" = {
           rate <- prior_spec$parameters$rate %||% 1
           stats::rexp(1, rate = rate)
         },
         "lognormal" = {
           meanlog <- prior_spec$parameters$meanlog %||% 0
           sdlog <- prior_spec$parameters$sdlog %||% 1
           stats::rlnorm(1, meanlog = meanlog, sdlog = sdlog)
         },
         
         # Default fallback (unchanged)
         stats::rnorm(1, 0, 1)
  )
}

#' Generate Response from Family Distribution
#' @keywords internal
.generate_response_from_family <- function(linear_pred, family_name, n_obs) {
  switch(family_name,
         "gaussian" = {
           error_sd <- abs(stats::rnorm(1, 0, 1))
           if (error_sd < 0.1) error_sd <- 0.1
           stats::rnorm(n_obs, linear_pred, error_sd)
         },
         "binomial" = {
           probs <- stats::plogis(linear_pred)
           stats::rbinom(n_obs, 1, probs)
         },
         "poisson" = {
           lambda <- exp(pmax(pmin(linear_pred, 10), -10))
           stats::rpois(n_obs, lambda)
         },
         "asymmetric_laplace" = {
           scale <- abs(stats::rnorm(1, 0, 1))
           stats::rnorm(n_obs, linear_pred, pmax(scale, 0.1))
         },
         {
           error_sd <- abs(stats::rnorm(1, 0, 1))
           stats::rnorm(n_obs, linear_pred, pmax(error_sd, 0.1))
         }
  )
}

#' Create Prior-Only Plot
#' @keywords internal
.create_prior_only_plot <- function(yrep, type = "dens_overlay", ndraws = NULL, verbose = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }
  
  title_text <- "Prior Predictive Check"
  col_rep <- "#6497b1"
  
  p <- switch(type,
              "dens_overlay" = {
                k <- min(20, nrow(yrep))
                rep_list <- lapply(seq_len(k), function(i) {
                  di <- stats::density(as.numeric(yrep[i, ]))
                  data.frame(x = di$x, d = di$y, draw = i, stringsAsFactors = FALSE)
                })
                df_rep <- do.call(rbind, rep_list)
                
                ggplot2::ggplot(df_rep, ggplot2::aes(x = x, y = d, group = draw)) +
                  ggplot2::geom_line(linewidth = 0.5, alpha = 0.3, colour = col_rep) +
                  ggplot2::labs(title = title_text, x = "Value", y = "Density") +
                  ggplot2::theme_minimal()
              },
              "hist" = {
                df <- data.frame(value = as.numeric(yrep), stringsAsFactors = FALSE)
                ggplot2::ggplot(df, ggplot2::aes(x = value)) +
                  ggplot2::geom_histogram(bins = 30, fill = col_rep, colour = "white", alpha = 0.7) +
                  ggplot2::labs(title = title_text, x = "Value", y = "Count") +
                  ggplot2::theme_minimal()
              }
  )
  
  class(p) <- c("qbrms_prior_check", class(p))
  p
}

# =============================================================================
# SECTION 5: DENSITY PLOT HELPERS
# =============================================================================

#' Extract Response Distribution Densities
#' @keywords internal
.extract_response_densities <- function(object, show_prior, show_data,
                                        ndraws, prior_ndraws, verbose) {
  
  plot_densities <- list()
  
  if (verbose) cat("Extracting posterior predictive samples...\n")
  yrep_posterior <- tryCatch({
    generate_posterior_predictive_samples(object, ndraws)
  }, error = function(e) {
    if (verbose) cat("Warning: Could not generate posterior samples\n")
    NULL
  })
  
  if (!is.null(yrep_posterior)) {
    post_values <- as.numeric(yrep_posterior)
    plot_densities$posterior <- data.frame(
      value = post_values,
      type  = "Posterior",
      stringsAsFactors = FALSE
    )
  }
  
  if (show_prior) {
    if (verbose) cat("Extracting prior predictive samples...\n")
    yrep_prior <- if (!is.null(object$prior_samples)) {
      n_use <- min(prior_ndraws, nrow(object$prior_samples))
      object$prior_samples[seq_len(n_use), , drop = FALSE]
    } else {
      tryCatch({
        generate_prior_predictive_samples(
          formula = object$original_formula,
          data    = object$data,
          family  = object$family,
          prior   = object$prior_specs,
          ndraws  = prior_ndraws,
          verbose = FALSE
        )
      }, error = function(e) NULL)
    }
    
    if (!is.null(yrep_prior)) {
      prior_values <- as.numeric(yrep_prior)
      plot_densities$prior <- data.frame(
        value = prior_values,
        type  = "Prior",
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (show_data) {
    if (verbose) cat("Extracting observed data...\n")
    response_var <- all.vars(object$original_formula)[1]
    if (response_var %in% names(object$data)) {
      y_obs <- object$data[[response_var]]
      if (is.factor(y_obs)) y_obs <- as.numeric(y_obs)
      plot_densities$data <- data.frame(
        value = as.numeric(y_obs),
        type  = "Observed Data",
        stringsAsFactors = FALSE
      )
    }
  }
  
  plot_densities
}

#' Extract Parameter Distribution Densities
#' @keywords internal
.extract_parameter_densities <- function(object, parameter, show_prior,
                                         ndraws, prior_ndraws, verbose) {
  
  plot_densities <- list()
  
  if (!is.null(object$fit$summary.fixed) && parameter %in% rownames(object$fit$summary.fixed)) {
    param_row <- object$fit$summary.fixed[parameter, ]
    mean_val <- param_row[["mean"]]
    sd_val   <- param_row[["sd"]]
    
    if (!is.na(mean_val) && !is.na(sd_val) && sd_val > 0) {
      post_samples <- stats::rnorm(ndraws, mean_val, sd_val)
      plot_densities$posterior <- data.frame(
        value = post_samples,
        type  = "Posterior",
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (show_prior && !is.null(object$prior_specs)) {
    prior_specs <- .extract_prior_specs_standalone(object$prior_specs, parameter, verbose = FALSE)
    
    if (parameter %in% names(prior_specs)) {
      prior_spec <- prior_specs[[parameter]]
      prior_samples <- replicate(prior_ndraws, .sample_from_prior_safe(prior_spec))
      
      plot_densities$prior <- data.frame(
        value = prior_samples,
        type  = "Prior",
        stringsAsFactors = FALSE
      )
    }
  }
  
  plot_densities
}

#' Create Base Density Plot
#' @keywords internal
.create_density_base_plot <- function(plot_densities, colours, alpha_levels, parameter, verbose) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for density plotting")
  }
  
  all_data <- do.call(rbind, plot_densities)
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No data available for plotting")
  }
  
  p <- ggplot2::ggplot(all_data, ggplot2::aes(x = value, colour = type, fill = type)) +
    ggplot2::geom_density(alpha = 0.2, linewidth = 1)
  
  colour_values <- c()
  fill_values   <- c()
  
  if ("Posterior" %in% unique(all_data$type)) {
    colour_values["Posterior"] <- colours$posterior %||% "#1F78B4"
    fill_values["Posterior"]   <- colours$posterior %||% "#1F78B4"
  }
  if ("Prior" %in% unique(all_data$type)) {
    colour_values["Prior"] <- colours$prior %||% "#E31A1C"
    fill_values["Prior"]   <- colours$prior %||% "#E31A1C"
  }
  if ("Observed Data" %in% unique(all_data$type)) {
    colour_values["Observed Data"] <- colours$data %||% "#000000"
    fill_values["Observed Data"]   <- colours$data %||% "#000000"
  }
  
  p <- p +
    ggplot2::scale_colour_manual(values = colour_values, name = "Distribution") +
    ggplot2::scale_fill_manual(values   = fill_values,   name = "Distribution")
  
  x_label <- if (is.null(parameter)) "Response Value" else parameter
  title_text <- if (is.null(parameter)) "Response Distribution" else paste("Parameter:", parameter)
  
  p <- p + ggplot2::labs(
    title = title_text,
    x     = x_label,
    y     = "Density"
  ) + ggplot2::theme_minimal()
  
  if (length(unique(all_data$type)) > 1) {
    p <- p + ggplot2::theme(legend.position = "bottom")
  }
  
  if (verbose) cat("Base density plot created; ready for ggplot2 customisation\n")
  
  p
}

# =============================================================================
# SECTION 6: WORKFLOW HELPER FUNCTIONS
# =============================================================================

#' Complete Prior-to-Posterior Workflow
#'
#' Fit a model with priors sampled, then produce a comparison density plot that
#' overlays posterior, prior, and observed distributions where available.
#'
#' @param formula Model formula.
#' @param data Data frame.
#' @param family Model family (default \code{gaussian()}).
#' @param prior Prior specification (default \code{NULL}).
#' @param verbose Logical; print progress messages.
#' @param ... Additional arguments forwarded to \code{qbrms()}.
#'
#' @return A list of class \code{qbrms_workflow_result} with elements \code{fit} and \code{plot}.
#' @export
prior_to_posterior_workflow <- function(formula,
                                        data,
                                        family = gaussian(),
                                        prior = NULL,
                                        verbose = TRUE,
                                        ...) {
  
  if (verbose) cat("=== Prior-to-Posterior Workflow ===\n")
  
  if (verbose) cat("1. Fitting model with prior samples...\n")
  fit <- qbrms(formula = formula,
               data    = data,
               family  = family,
               prior   = prior,
               sample_prior = "yes",
               verbose = verbose,
               ...)
  
  if (verbose) cat("2. Creating comparison plot...\n")
  plot <- density_plot(fit, show_prior = TRUE, show_data = TRUE, verbose = verbose)
  
  if (verbose) cat("3. Workflow complete!\n")
  
  result <- list(
    fit     = fit,
    plot    = plot,
    formula = formula,
    family  = family,
    prior   = prior
  )
  
  class(result) <- c("qbrms_workflow_result", "list")
  result
}

#' @method print qbrms_workflow_result
#' @export
print.qbrms_workflow_result <- function(x, ...) {
  cat("qbrms Prior-to-Posterior Workflow Result\n\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Family: ", extract_family_name(x$family), "\n")
  cat("Data:   ", nrow(x$fit$data), "observations\n\n")
  cat("Components:\n")
  cat("  $fit  - Fitted qbrms model\n")
  cat("  $plot - Comparison plot\n\n")
  cat("Use summary(result$fit) to see model results.\n")
  cat("Use result$plot to display the comparison plot.\n")
  invisible(x)
}

#' Plot Parameters with Prior/Posterior Comparison
#'
#' @description
#' Create density plots for multiple model parameters, optionally comparing
#' posterior estimates with their priors. Returns a ggplot2 object with
#' faceted parameter plots.
#'
#' @param object A \code{qbrms_fit} object.
#' @param pars Optional character vector of parameter names to plot. If \code{NULL},
#'   plots all fixed-effect parameters.
#' @param show_prior Logical; if \code{TRUE}, overlay prior distributions.
#' @param ndraws Number of posterior draws to use for plotting.
#' @param prior_ndraws Number of prior draws to use if \code{show_prior = TRUE}.
#' @param ncol Number of columns for faceting (default 2).
#' @param alpha_levels Numeric vector of length 2 giving alpha levels for 
#'   c(posterior, prior). Default c(0.8, 0.5).
#' @param colours Character vector of length 2 giving colours for 
#'   c(posterior, prior). Default c("#1F78B4", "#E31A1C").
#' @param verbose Logical; print progress information.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object with faceted parameter density plots.
#' @export
#' @examples
#' \dontrun{
#' fit <- qbrms(y ~ x1 + x2, data = my_data, sample_prior = "yes")
#' 
#' # Plot all parameters
#' plot_parameters(fit)
#' 
#' # Plot specific parameters with priors
#' plot_parameters(fit, pars = c("x1", "x2"), show_prior = TRUE)
#' 
#' # Customize appearance
#' plot_parameters(fit, show_prior = TRUE) + 
#'   theme_bw() + 
#'   labs(title = "My Parameter Estimates")
#' }
plot_parameters <- function(object,
                            pars = NULL,
                            show_prior = FALSE,
                            ndraws = 200,
                            prior_ndraws = 200,
                            ncol = 2,
                            alpha_levels = c(0.8, 0.5),
                            colours = c("#1F78B4", "#E31A1C"),
                            verbose = TRUE,
                            ...) {
  
  if (!inherits(object, "qbrms_fit")) {
    stop("object must be a qbrms_fit object")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for parameter plotting")
  }
  
  # Extract parameter samples
  posterior_samples <- .ba_extract_parameter_samples(object, ndraws)
  
  # Determine which parameters to plot
  if (is.null(pars)) {
    pars <- names(posterior_samples)
  } else {
    # Check that requested parameters exist
    missing_pars <- setdiff(pars, names(posterior_samples))
    if (length(missing_pars) > 0) {
      warning("Parameters not found: ", paste(missing_pars, collapse = ", "))
      pars <- intersect(pars, names(posterior_samples))
    }
  }
  
  if (length(pars) == 0) {
    stop("No valid parameters found for plotting")
  }
  
  if (verbose) {
    cat("Plotting", length(pars), "parameters:", paste(pars, collapse = ", "), "\n")
  }
  
  # Prepare posterior data
  post_data_list <- lapply(pars, function(p) {
    data.frame(
      Parameter = p,
      Value = posterior_samples[[p]],
      Type = "Posterior",
      stringsAsFactors = FALSE
    )
  })
  plot_data <- do.call(rbind, post_data_list)
  
  # Add prior data if requested
  if (show_prior) {
    if (verbose) cat("Adding prior distributions...\n")
    
    # Generate prior samples
    prior_data_list <- lapply(pars, function(p) {
      # Try to get prior samples from object first
      if (!is.null(object$prior_samples)) {
        # Use existing prior samples if available
        X <- try(.qbrms_model_matrix_fixed(object, object$data), silent = TRUE)
        if (!inherits(X, "try-error") && p %in% colnames(X)) {
          # This is a rough approximation - in practice you'd want more sophisticated prior extraction
          prior_spec <- list(distribution = "normal", parameters = list(mean = 0, sd = 1))
          if (p == "(Intercept)") {
            prior_spec$parameters$sd <- 2.5
          }
          prior_samples <- replicate(prior_ndraws, .sample_from_prior_safe(prior_spec))
        } else {
          # Fallback to weakly informative prior
          prior_samples <- stats::rnorm(prior_ndraws, 0, 1)
        }
      } else {
        # Generate default prior samples
        prior_spec <- if (p == "(Intercept)") {
          list(distribution = "normal", parameters = list(mean = 0, sd = 2.5))
        } else {
          list(distribution = "normal", parameters = list(mean = 0, sd = 1))
        }
        prior_samples <- replicate(prior_ndraws, .sample_from_prior_safe(prior_spec))
      }
      
      data.frame(
        Parameter = p,
        Value = prior_samples,
        Type = "Prior",
        stringsAsFactors = FALSE
      )
    })
    
    prior_data <- do.call(rbind, prior_data_list)
    plot_data <- rbind(plot_data, prior_data)
  }
  
  # Create the plot
  plot_data$Parameter <- factor(plot_data$Parameter, levels = pars)
  plot_data$Type <- factor(plot_data$Type, levels = if (show_prior) c("Posterior", "Prior") else "Posterior")
  
  # Set up colours and alpha
  if (length(colours) < 2) colours <- rep(colours[1], 2)
  if (length(alpha_levels) < 2) alpha_levels <- rep(alpha_levels[1], 2)
  
  color_values <- c("Posterior" = colours[1])
  alpha_values <- c("Posterior" = alpha_levels[1])
  
  if (show_prior) {
    color_values["Prior"] <- colours[2]
    alpha_values["Prior"] <- alpha_levels[2]
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, fill = Type, colour = Type)) +
    ggplot2::geom_density(alpha = 0.3, linewidth = 1) +
    ggplot2::facet_wrap(~ Parameter, scales = "free", ncol = ncol) +
    ggplot2::scale_fill_manual(values = color_values, name = "Distribution") +
    ggplot2::scale_colour_manual(values = color_values, name = "Distribution") +
    ggplot2::labs(
      title = if (show_prior) "Parameter Estimates: Prior vs Posterior" else "Parameter Estimates",
      x = "Parameter Value",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Add reference line at zero
  p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", 
                               colour = "grey50", alpha = 0.7)
  
  if (verbose) cat("Parameter plot created successfully\n")
  
  class(p) <- c("qbrms_parameter_plot", class(p))
  return(p)
}

# Helper function for safe prior sampling (if not already defined)
if (!exists(".sample_from_prior_safe")) {
  .sample_from_prior_safe <- function(prior_spec) {
    switch(prior_spec$distribution,
           "normal" = stats::rnorm(1, 
                                   prior_spec$parameters$mean %||% 0, 
                                   prior_spec$parameters$sd %||% 1),
           "student_t" = {
             df <- prior_spec$parameters$df %||% 3
             location <- prior_spec$parameters$location %||% 0
             scale <- prior_spec$parameters$scale %||% 1
             location + scale * stats::rt(1, df)
           },
           "cauchy" = {
             location <- prior_spec$parameters$location %||% 0
             scale <- prior_spec$parameters$scale %||% 1
             stats::rcauchy(1, location, scale)
           },
           "uniform" = {
             min_val <- prior_spec$parameters$min %||% -10
             max_val <- prior_spec$parameters$max %||% 10
             stats::runif(1, min_val, max_val)
           },
           stats::rnorm(1, 0, 1)  # Default fallback
    )
  }
}
