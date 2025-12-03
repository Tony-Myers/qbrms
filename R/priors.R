# =============================================================================
# R/priors.R
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

#' Specify Beta Prior Distribution
#' @param alpha First shape parameter
#' @param beta Second shape parameter
#' @return A prior distribution object
#' @export
beta_prior <- function(alpha = 1, beta = 1) {
  if (alpha <= 0 || beta <= 0) stop("Beta parameters must be positive")
  structure(list(
    distribution = "beta",
    parameters = list(alpha = alpha, beta = beta)
  ), class = "qbrms_prior_dist")
}

# =============================================================================
# HYBRID CONSTRUCTORS (Family / Prior Polymorphism)
# =============================================================================

#' Student's t Distribution (Prior or Family)
#'
#' @description
#' Functions that act as both family constructors (for `qbrm`) and prior
#' specifications (for `prior`), depending on arguments.
#'
#' @param link_or_df For family: link function (character). For prior: degrees of freedom (numeric).
#' @param location Location parameter (prior only).
#' @param scale Scale parameter (prior only).
#' @param link Optional link function (if acting as family).
#' @param link.sigma Link for sigma (family only).
#' @param link.nu Link for nu (family only).
#' @param ... Additional arguments.
#' @return A family object or prior object depending on inputs.
#' @export
student_t <- function(link_or_df = "identity", location = 0, scale = 1, link = NULL, link.sigma = "log", link.nu = "log", ...) {
  if (!is.null(link)) {
    return(structure(list(family = "student_t", link = link, link.sigma = link.sigma, link.nu = link.nu), class = "family"))
  }
  if (is.character(link_or_df)) {
    return(structure(list(family = "student_t", link = link_or_df, link.sigma = link.sigma, link.nu = link.nu), class = "family"))
  }
  structure(list(
    distribution = "student_t",
    parameters = list(df = link_or_df, location = location, scale = scale)
  ), class = "qbrms_prior_dist")
}

#' @rdname student_t
#' @export
student_t_prior <- function(link_or_df = 3, location = 0, scale = 1, link.sigma = "log", link.nu = "log") {
  student_t(link_or_df, location, scale, link.sigma, link.nu)
}

#' Lognormal Distribution (Prior or Family)
#' @param meanlog_or_link Mean on log scale (numeric) or link function (character).
#' @param sdlog SD on log scale (numeric).
#' @param link Optional link function (if acting as family).
#' @param ... Additional arguments.
#' @return A family object or prior object depending on inputs.
#' @export
lognormal <- function(meanlog_or_link = "identity", sdlog = 1, link = NULL, ...) {
  if (!is.null(link)) {
    return(structure(list(family = "lognormal", link = link), class = "family"))
  }
  if (is.character(meanlog_or_link)) {
    return(structure(list(family = "lognormal", link = meanlog_or_link), class = "family"))
  }
  if (sdlog <= 0) stop("Lognormal sdlog must be positive")
  structure(list(
    distribution = "lognormal",
    parameters = list(meanlog = meanlog_or_link, sdlog = sdlog)
  ), class = "qbrms_prior_dist")
}

#' @rdname lognormal
#' @export
lognormal_prior <- function(meanlog_or_link = 0, sdlog = 1) {
  lognormal(meanlog_or_link, sdlog)
}

#' Exponential Distribution (Prior or Family)
#' @param rate_or_link Rate parameter (numeric) or link function (character).
#' @param link Optional link function (if acting as family).
#' @param ... Additional arguments.
#' @return A family object or prior object depending on inputs.
#' @export
exponential <- function(rate_or_link = "log", link = NULL, ...) {
  if (!is.null(link)) {
    return(structure(list(family = "exponential", link = link), class = "family"))
  }
  if (is.character(rate_or_link)) {
    return(structure(list(family = "exponential", link = rate_or_link), class = "family"))
  }
  if (rate_or_link <= 0) stop("Rate parameter must be positive")
  structure(list(
    distribution = "exponential",
    parameters = list(rate = rate_or_link)
  ), class = "qbrms_prior_dist")
}

#' @rdname exponential
#' @export
prior_exponential <- function(rate_or_link = 1) {
  exponential(rate_or_link)
}

#' Gamma Distribution (Prior or Family)
#' @param shape_or_link Shape parameter (numeric) or link function (character).
#' @param rate Rate parameter.
#' @param link Optional link function (if acting as family).
#' @param ... Additional arguments.
#' @return A family object or prior object depending on inputs.
#' @export
gamma <- function(shape_or_link = "log", rate = 1, link = NULL, ...) {
  if (!is.null(link)) {
    return(structure(list(family = "gamma", link = link), class = "family"))
  }
  if (is.character(shape_or_link)) {
    return(structure(list(family = "gamma", link = shape_or_link), class = "family"))
  }
  gamma_prior(shape_or_link, rate)
}

#' @rdname gamma
#' @export
gamma_prior <- function(shape_or_link = 2, rate = 1) {
  if (shape_or_link <= 0 || rate <= 0) stop("Gamma parameters must be positive")
  structure(list(
    distribution = "gamma",
    parameters = list(shape = shape_or_link, rate = rate)
  ), class = "qbrms_prior_dist")
}

# =============================================================================
# PRIOR SPECIFICATION WRAPPER
# =============================================================================

#' Specify Prior for Model Parameters
#' @param prior A prior distribution object.
#' @param class Parameter class ("Intercept", "b", "sd", etc.)
#' @param coef Specific coefficient name (optional)
#' @param group Specific group name for random effects (optional)
#' @return A prior specification object
#' @export
prior <- function(prior, class = "b", coef = NULL, group = NULL) {
  if (is.character(prior)) {
    return(structure(list(prior = prior, class = class, coef = coef, group = group), 
                     class = "qbrms_prior_spec"))
  }
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
c.qbrms_prior_spec <- function(...) {
  priors <- list(...)
  structure(priors, class = "qbrms_prior_list")
}

#' Print Prior Distribution Objects
#' @param x A qbrms_prior_dist object
#' @param ... Unused
#' @return Invisibly returns \code{x}.
#' @export
print.qbrms_prior_dist <- function(x, ...) {
  params <- x$parameters
  param_str <- paste(names(params), "=", params, collapse = ", ")
  cat(sprintf("%s(%s)\n", x$distribution, param_str))
  invisible(x)
}

#' Print Prior Specification Objects  
#' @param x A qbrms_prior_spec object
#' @param ... Unused
#' @return Invisibly returns \code{x}.
#' @export
print.qbrms_prior_spec <- function(x, ...) {
  if (!is.null(x$prior) && is.character(x$prior)) {
    cat(sprintf("prior(prior = \"%s\", class = \"%s\")\n", x$prior, x$class))
    return(invisible(x))
  }
  params <- x$parameters
  param_str <- paste(names(params), "=", params, collapse = ", ")
  coef_str <- if (!is.null(x$coef)) paste0(", coef = ", x$coef) else ""
  cat(sprintf("prior(%s(%s), class = %s%s)\n", 
              x$distribution, param_str, x$class, coef_str))
  invisible(x)
}

#' Print Prior List Objects
#' @param x A qbrms_prior_list object
#' @param ... Unused
#' @return Invisibly returns \code{x}.
#' @export
print.qbrms_prior_list <- function(x, ...) {
  cat("Prior specifications:\n")
  for (i in seq_along(x)) {
    cat(sprintf("%d. ", i))
    print(x[[i]])
  }
  invisible(x)
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
#' @export
get_default_prior <- function(class) {
  switch(class,
         "Intercept" = prior(normal(0, 2.5), class = "Intercept"),
         "b" = prior(normal(0, 1), class = "b"),
         "sd" = prior(cauchy(0, 1), class = "sd"),
         prior(normal(0, 1), class = class)
  )
}

# Simple null coalescing operator if not defined elsewhere
`%||%` <- function(x, y) if (is.null(x)) y else x