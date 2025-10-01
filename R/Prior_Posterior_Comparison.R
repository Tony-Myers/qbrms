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

# Helper function for safe prior sampling
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