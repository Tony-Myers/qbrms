# =============================================================================
# Prior Specification Functions for qbrms
# =============================================================================

#' Prior Distribution Specifications
#' @name priors
#' @description Functions to specify prior distributions for qbrms models
NULL

#' Specify Normal Prior Distribution
#' @param mean Mean of the normal distribution (default 0)
#' @param sd Standard deviation of the normal distribution (default 1)
#' @return A prior distribution object
#' @export
normal <- function(mean = 0, sd = 1) {
  structure(list(
    distribution = "normal",
    parameters = list(mean = mean, sd = sd)
  ), class = "qbrms_prior_dist")
}

#' Specify Student's t Prior Distribution
#' @param df Degrees of freedom (default 3)
#' @param location Location parameter (default 0)
#' @param scale Scale parameter (default 1)
#' @return A prior distribution object
#' @export
student_t_prior <- function(df = 3, location = 0, scale = 1) {
  structure(list(
    distribution = "student_t",
    parameters = list(df = df, location = location, scale = scale)
  ), class = "qbrms_prior_dist")
}

#' Specify Cauchy Prior Distribution
#' @param location Location parameter (default 0)
#' @param scale Scale parameter (default 1)
#' @return A prior distribution object
#' @export
cauchy <- function(location = 0, scale = 1) {
  structure(list(
    distribution = "cauchy",
    parameters = list(location = location, scale = scale)
  ), class = "qbrms_prior_dist")
}

#' Specify Uniform Prior Distribution
#' @param min Lower bound (default -Inf)
#' @param max Upper bound (default Inf)
#' @return A prior distribution object
#' @export
uniform <- function(min = -Inf, max = Inf) {
  structure(list(
    distribution = "uniform",
    parameters = list(min = min, max = max)
  ), class = "qbrms_prior_dist")
}

#' Specify Exponential Prior Distribution
#' @param rate Rate parameter for the exponential distribution (default 1)
#' @return A prior distribution object
#' @export
prior_exponential <- function(rate = 1) {
  if (rate <= 0) {
    stop("Rate parameter must be positive")
  }
  
  structure(list(
    distribution = "exponential",
    parameters = list(rate = rate)
  ), class = "qbrms_prior_dist")
}

#' Specify Gamma Prior Distribution  
#' @param shape Shape parameter (must be positive)
#' @param rate Rate parameter (must be positive) 
#' @return A prior distribution object
#' @export
gamma_prior <- function(shape = 2, rate = 1) {
  if (shape <= 0 || rate <= 0) {
    stop("Gamma parameters must be positive")
  }
  structure(list(
    distribution = "gamma",
    parameters = list(shape = shape, rate = rate)
  ), class = "qbrms_prior_dist")
}

#' Specify Beta Prior Distribution
#' @param alpha First shape parameter (must be positive)
#' @param beta Second shape parameter (must be positive)
#' @return A prior distribution object  
#' @export
beta_prior <- function(alpha = 1, beta = 1) {
  if (alpha <= 0 || beta <= 0) {
    stop("Beta parameters must be positive")
  }
  structure(list(
    distribution = "beta", 
    parameters = list(alpha = alpha, beta = beta)
  ), class = "qbrms_prior_dist")
}

#' Specify Lognormal Prior Distribution  
#' @param meanlog Mean of log scale
#' @param sdlog Standard deviation of log scale (must be positive)
#' @return A prior distribution object
#' @export
lognormal_prior <- function(meanlog = 0, sdlog = 1) {
  if (sdlog <= 0) {
    stop("Lognormal sdlog must be positive")
  }
  structure(list(
    distribution = "lognormal",
    parameters = list(meanlog = meanlog, sdlog = sdlog)
  ), class = "qbrms_prior_dist")
}

#' Specify Prior for Model Parameters
#' @param prior A prior distribution object created by functions like normal(), student_t(), etc.
#' @param class Parameter class ("Intercept", "b" for fixed effects, "sd" for random effects, etc.)
#' @param coef Specific coefficient name (optional)
#' @param group Specific group name for random effects (optional)
#' @return A prior specification object
#' @export
#' @examples
#' \dontrun{
#' # Normal prior for all fixed effects
#' prior(normal(0, 1), class = "b")
#' 
#' # Student's t prior for intercept
#' prior(student_t(3, 0, 2.5), class = "Intercept")
#' 
#' # Specific prior for a coefficient
#' prior(normal(0, 0.5), class = "b", coef = "age")
#' }
prior <- function(prior, class = "b", coef = NULL, group = NULL) {
  if (!inherits(prior, "qbrms_prior_dist")) {
    stop("prior must be a distribution object created by normal(), student_t(), etc.")
  }
  
  structure(list(
    distribution = prior$distribution,
    parameters = prior$parameters,
    class = class,
    coef = coef,
    group = group
  ), class = "qbrms_prior_spec")
}

#' Combine Multiple Prior Specifications
#' @param ... Prior specification objects created by prior()
#' @return A combined prior object
#' @export
#' @examples
#' \dontrun{
#' # Multiple priors
#' priors <- c(
#'   prior(normal(0, 2.5), class = "Intercept"),
#'   prior(normal(0, 1), class = "b"),
#'   prior(cauchy(0, 1), class = "sd")
#' )
#' }
c.qbrms_prior_spec <- function(...) {
  priors <- list(...)
  structure(priors, class = "qbrms_prior_list")
}

#' Print Prior Distribution Objects - Enhanced for all distributions
#' @param x A qbrms_prior_dist object
#' @param ... Unused
#' @export
print.qbrms_prior_dist <- function(x, ...) {
  params <- x$parameters
  
  # Format parameter display based on distribution for better readability
  param_str <- switch(x$distribution,
                      "gamma" = paste0("shape=", params$shape, ", rate=", params$rate),
                      "beta" = paste0("alpha=", params$alpha, ", beta=", params$beta), 
                      "exponential" = paste0("rate=", params$rate),
                      "lognormal" = paste0("meanlog=", params$meanlog, ", sdlog=", params$sdlog),
                      "student_t" = paste0("df=", params$df, ", location=", params$location, ", scale=", params$scale),
                      "normal" = paste0("mean=", params$mean, ", sd=", params$sd),
                      "cauchy" = paste0("location=", params$location, ", scale=", params$scale),
                      "uniform" = paste0("min=", params$min, ", max=", params$max),
                      # Default fallback for any distribution not explicitly handled
                      paste(names(params), "=", params, collapse = ", ")
  )
  
  cat(sprintf("%s(%s)\n", x$distribution, param_str))
}

#' Print Prior Specification Objects  
#' @param x A qbrms_prior_spec object
#' @param ... Unused
#' @export
print.qbrms_prior_spec <- function(x, ...) {
  params <- x$parameters
  param_str <- paste(names(params), "=", params, collapse = ", ")
  
  coef_str <- if (!is.null(x$coef)) paste0(", coef = ", x$coef) else ""
  group_str <- if (!is.null(x$group)) paste0(", group = ", x$group) else ""
  
  cat(sprintf("prior(%s(%s), class = %s%s%s)\n", 
              x$distribution, param_str, x$class, coef_str, group_str))
}

#' Print Prior List Objects
#' @param x A qbrms_prior_list object  
#' @param ... Unused
#' @export
print.qbrms_prior_list <- function(x, ...) {
  cat("Prior specifications:\n")
  for (i in seq_along(x)) {
    cat(sprintf("%d. ", i))
    print(x[[i]])
  }
}

#' Default Priors for qbrms Models
#' @return A default prior list
#' @export
default_priors <- function() {
  c(
    prior(normal(0, 2.5), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(cauchy(0, 1), class = "sd")
  )
}

#' Get Default Prior for Parameter Class
#' @param class Parameter class
#' @return A default prior for that class
#' @keywords internal
get_default_prior <- function(class) {
  switch(class,
         "Intercept" = prior(normal(0, 2.5), class = "Intercept"),
         "b" = prior(normal(0, 1), class = "b"),
         "sd" = prior(cauchy(0, 1), class = "sd"),
         prior(normal(0, 1), class = class)
  )
}

#' Apply Priors to INLA Model (placeholder)
#' @param priors Prior specifications
#' @param formula Model formula
#' @param family Model family
#' @return Modified control settings for INLA (placeholder implementation)
#' @keywords internal
apply_priors_to_inla <- function(priors, formula, family) {
  # Placeholder - INLA doesn't directly support arbitrary priors like brms
  # This would need custom implementation for each prior type
  
  if (is.null(priors)) return(list())
  
  # For now, just return default INLA settings
  # Full implementation would translate qbrms priors to INLA hyperprior specifications
  warning("Custom priors not fully implemented yet - using INLA defaults")
  
  list(
    control.fixed = list(mean.intercept = 0, prec.intercept = 1, mean = 0, prec = 1),
    control.family = list()
  )
}