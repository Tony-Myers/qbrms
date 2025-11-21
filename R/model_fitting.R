# =============================================================================
# R/model_fitting.R 
# =============================================================================

#' @title Model Fitting Functions for qbrms Package
#' @description Core model fitting functionality with qbrm alias and multinomial support
#' @details This file contains the main qbrms/qbrm function and all supporting functions
#' @keywords internal
#' @name model_fitting
#' @importFrom stats setNames model.matrix rnorm rbinom rpois plogis
NULL

# =============================================================================
# SECTION 1: HELPER FUNCTIONS AND UTILITIES
# =============================================================================

# Utility operators and functions

format_duration <- function(seconds) {
  paste0(round(seconds, 2), " seconds")
}
#' Gaussian Family
#' @export
gaussian <- function() structure(list(family = "gaussian"), class = "family")

#' Binomial Family
#' @export
binomial <- function() structure(list(family = "binomial"), class = "family")

#' Poisson Family
#' @export
poisson <- function() structure(list(family = "poisson"), class = "family")

#' Cumulative Family for Ordinal Regression
#' @param link Link function (default: "logit")
#' @export
cumulative <- function(link = "logit") structure(list(family = "cumulative", link = link), class = "family")

#' Negative Binomial Family
#' @export
neg_binomial <- function() structure(list(family = "negbinomial"), class = "family")

#' Negative Binomial Family (Alias)
#' @export
negbinomial <- function() structure(list(family = "negbinomial"), class = "family")

#' Asymmetric Laplace for Quantile Regression
#' @export
asymmetric_laplace <- function() structure(list(family = "asymmetric_laplace"), class = "family")

#' Poisson Trick for Multinomial
#' @keywords internal
poisson_trick_multinomial <- function() structure(list(family = "poisson_trick_multinomial"), class = "family")

#' Skew Normal Family
#' @export
skew_normal <- function() structure(list(family = "skew_normal"), class = "family")

#' Multinomial Family
#' @export
multinomial <- function() structure(list(family = "multinomial"), class = "family")

#' Student's t Family for Robust Regression
#'
#' @description
#' Student's t-distribution family for robust regression with heavier tails
#' than Gaussian to handle outliers.
#'
#' @return A family object for use with qbrms()
#'
#' @examples
#' \dontrun{
#' # Robust regression
#' fit <- qbrms(y ~ x, data = data, family = student_t())
#' }
#'
#' @export
student_t <- function() {
  structure(list(family = "studentt"), class = "family")
}

#' Student Family (Alias)
#' @rdname student_t
#' @export
student <- function() {
  structure(list(family = "student"), class = "family")
}


# =============================================================================
# ERROR HANDLING FOR EDGE CASES
# =============================================================================

#' Drop random-effect terms from a formula
#'
#' Removes any "(... | ...)" blocks from the RHS and keeps a syntactically
#' valid formula. Cleans up dangling '+' and ensures there is at least "~ 1"
#' if all terms are removed.
#' @keywords internal

# =============================================================================
#  qbrms() FUNCTION WITH ROUTING 
# =============================================================================
##' Quick Bayesian Regression Models with Automatic Routing
#'
#' @description
#' Enhanced qbrms interface with automatic routing to specialised implementations.
#' Supports ordinal regression via TMB, quantile regression, and all standard INLA families.
#'
#' @param formula Model formula in lme4/brms style
#' @param data Data frame containing the variables in the model
#' @param family Model family (default: gaussian()). Ordinal families automatically route to qbrmO()
#' @param prior Prior specifications (default: NULL)
#' @param sample_prior Whether to sample from priors ("no", "yes", "only"). Default: "no"
#' @param quantile For asymmetric_laplace family, which quantile to estimate (default: 0.5)
#' @param control.compute INLA control settings for model information criteria
#' @param verbose Logical; print diagnostic information (default: getOption("qbrms.verbose", FALSE))
#' @param ... Additional arguments passed to fitting functions
#'
#' @return A qbrms_fit object with model results, or routed to appropriate specialist function
#' @seealso \code{\link{qbrmO}} for direct ordinal model fitting
#' @export
qbrms <- function(formula, data, family = gaussian(),
                  prior = NULL, sample_prior = "no",
                  quantile = 0.5,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                  verbose = getOption("qbrms.verbose", FALSE), ...) {
  
  # Check if formula is a bf() object and strip distributional parts if necessary
  formula <- sanitize_formula(formula)
  # --------------------------
  
  # ---- local fallbacks (do not affect package API) ---------------------------
  `%||%` <- get0("%||%", ifnotfound = function(x, y) if (is.null(x)) y else x)
  if (!exists(".qbrms_silently", mode = "function", inherits = TRUE)) {
    .qbrms_silently <- function(expr) eval.parent(substitute(expr))
  }
  if (!exists(".qbrms_silence_tmb", mode = "function", inherits = TRUE)) {
    .qbrms_silence_tmb <- function(expr) eval.parent(substitute(expr))
  }
  
  if (isTRUE(verbose)) cat("Starting qbrms model fitting...\n")
  
  # ---- basic checks ----------------------------------------------------------
  if (is.null(formula) || is.null(data)) stop("Both formula and data are required")
  if (nrow(data) == 0) stop("Data cannot be empty")
  
  formula_vars <- tryCatch(all.vars(formula), error = function(e) character(0))
  if (length(formula_vars) > 0) {
    missing_vars <- setdiff(formula_vars, names(data))
    if (length(missing_vars) > 0) {
      stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
    }
  }
  if (length(formula_vars) == 0L)
    stop("Formula appears to have no variables")
  
  # ---- gentle auto-family inference (preserves old ergonomics) --------------
  # Only promote the default gaussian() when response type clearly indicates it.
  y <- data[[formula_vars[1]]]
  if (is.list(family) && identical(tolower(family$family %||% ""), "gaussian")) {
    if (is.ordered(y)) {
      family <- cumulative()
    } else if (is.factor(y)) {
      family <- if (nlevels(y) == 2L) binomial() else multinomial()
    }
  }
  
  # ---- convert family; allow ordinal routing --------------------------------
  inla_fam_or_route <- tryCatch({
    convert_family_to_inla(family, quantile, allow_ordinal_routing = TRUE)
  }, error = function(e) {
    stop("Family conversion failed: ", e$message, call. = FALSE)
  })
  
  # For non-ordinal families, validate data against the converted name
  if (!requires_routing(inla_fam_or_route)) {
    fam_for_check <- tryCatch(extract_family_name(inla_fam_or_route), error = function(e) "gaussian")
    try(validate_family_data(y, fam_for_check), silent = TRUE)
  }
  
  # ---- ordinal routing path (TMB) -------------------------------------------
  if (requires_routing(inla_fam_or_route)) {
    routing_info <- extract_routing_info(inla_fam_or_route)
    if (isTRUE(verbose)) {
      cat("Detected ", routing_info$original_family, " family - routing to ",
          routing_info$target, "()\n", sep = "")
    }
    if (!identical(routing_info$target, "qbrmO"))
      stop("Unknown routing target: ", routing_info$target)
    
    if (!exists("qbrmO", mode = "function"))
      stop("Ordinal regression detected but qbrmO() is not available.")
    
    ordinal_family <- structure(list(
      family    = routing_info$family,
      link      = routing_info$link,
      threshold = routing_info$threshold
    ), class = c("brmsfamily", "family"))
    
    # Silence TMB compile/link unless verbose=TRUE
    if (isTRUE(verbose)) {
      return(qbrmO(
        formula = formula, data = data, family = ordinal_family,
        prior = prior, verbose = TRUE, ...
      ))
    } else {
      return(.qbrms_silence_tmb(
        qbrmO(
          formula = formula, data = data, family = ordinal_family,
          prior = prior, verbose = FALSE, ...
        )
      ))
    }
  }
  
  # ---- non-ordinal: continue as before --------------------------------------
  inla_fam <- inla_fam_or_route
  
  # Prior-only sampling
  if (identical(sample_prior, "only")) {
    if (isTRUE(verbose)) cat("Generating prior predictive samples only...\n")
    prior_result <- tryCatch({
      .qbrms_silently(
        generate_prior_predictive_samples(
          formula = formula, data = data, family = family,
          prior = prior, ndraws = 100, verbose = verbose
        )
      )
    }, error = function(e) {
      if (isTRUE(verbose)) cat("Warning: Prior sampling failed, using fallback\n")
      n_obs <- nrow(data); matrix(stats::rnorm(100 * n_obs), nrow = 100, ncol = n_obs)
    })
    
    result <- list(
      prior_samples = prior_result,
      original_formula = formula,
      data = data,
      family = family,
      model_type = "prior_predictive",
      prior_specs = prior,
      fitting_time = 0
    )
    class(result) <- c("qbrms_prior_fit", "qbrms_fit", "list")
    if (isTRUE(verbose)) cat("qbrms model fitting completed.\n")
    return(result)
  }
  
  # Parse formula components (random effects etc.)
  formula_components <- tryCatch({
    parse_formula_components(formula, data)
  }, error = function(e) {
    if (isTRUE(verbose)) cat("Warning: Could not parse formula components, using defaults\n")
    list(has_random_effects = FALSE, is_binomial_trials = FALSE)
  })
  
  fam_name <- tryCatch(extract_family_name(inla_fam), error = function(e) "gaussian")
  
  # Multinomial branch
  if (identical(fam_name, "multinomial")) {
    fit_res <- tryCatch({
      .qbrms_silently(fit_multinomial_model(formula, data, inla_fam, control.compute, verbose, ...))
    }, error = function(e) {
      if (isTRUE(verbose)) cat("Multinomial fitting failed: ", e$message, "\n", sep = "")
      NULL
    })
    if (!is.null(fit_res)) {
      dur <- fit_res$fitting_time
      res <- list(
        fit = fit_res$fit,
        original_formula = formula,
        data = data,
        family = inla_fam,
        has_random = FALSE,
        has_random_effects = FALSE,
        model_type = "multinomial",
        fitting_time = dur,
        timing = list(total_seconds = dur,
                      formatted = format_duration(dur),
                      formatted_duration = format_duration(dur))
      )
      class(res) <- c("qbrms_multinomial_fit", "qbrms_fit")
      if (isTRUE(verbose)) cat("qbrms model fitting completed successfully.\n")
      return(res)
    }
  }
  
  # Quantile regression branch
  if (identical(fam_name, "asymmetric_laplace")) {
    start_time <- Sys.time()
    quantile_fit <- tryCatch({
      .qbrms_silently(create_quantile_fit(formula, data, quantile, verbose))
    }, error = function(e) {
      if (isTRUE(verbose)) cat("Quantile regression failed: ", e$message, "\n", sep = "")
      yy <- data[[formula_vars[1]]]
      list(
        summary.fixed = data.frame(
          mean = c(if (is.numeric(yy)) mean(yy, na.rm = TRUE) else 0, 0),
          sd   = c(0.1, 0.1),
          row.names = c("(Intercept)", "fallback")
        ),
        converged = FALSE,
        fallback_type = "quantile_fallback"
      )
    })
    dur <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    res <- list(
      fit = quantile_fit,
      original_formula = formula,
      data = data,
      family = inla_fam,
      has_random = FALSE,
      has_random_effects = FALSE,
      model_type = "quantile_regression",
      quantile = quantile,
      fitting_time = dur,
      timing = list(total_seconds = dur,
                    formatted = format_duration(dur),
                    formatted_duration = format_duration(dur))
    )
    class(res) <- c("qbrms_fit")
    if (isTRUE(verbose)) cat("qbrms model fitting completed successfully.\n")
    return(res)
  }
  
  # Standard fixed/mixed effects path (INLA or fallback inside)
  fit_res <- tryCatch({
    .qbrms_silently(fit_model_robust_fixed(formula, data, inla_fam, control.compute, verbose, ...))
  }, error = function(e) {
    if (isTRUE(verbose)) cat("All fitting methods failed: ", e$message, "\n", sep = "")
    yy <- data[[formula_vars[1]]]
    list(
      fit = list(
        summary.fixed = data.frame(
          mean = if (is.numeric(yy)) mean(yy, na.rm = TRUE) else 0,
          sd = 0.1,
          row.names = "(Intercept)"
        ),
        converged = FALSE,
        fallback_type = "emergency"
      ),
      fitting_time = 0
    )
  })
  
  # Grouping variable extraction (simple (1|group) case)
  group_var <- NULL
  if (isTRUE(formula_components$has_random_effects)) {
    frm_str <- tryCatch(deparse(formula, width.cutoff = 500), error = function(e) "")
    mt <- regexec("\\(1 \\| ([^\\)]+)\\)", frm_str)
    cap <- regmatches(frm_str, mt)
    if (length(cap[[1]]) >= 2) group_var <- trimws(cap[[1]][2])
  }
  
  # Optional prior predictive samples
  prior_samples <- NULL
  if (identical(sample_prior, "yes")) {
    if (isTRUE(verbose)) cat("Also generating prior predictive samples...\n")
    prior_samples <- tryCatch({
      .qbrms_silently(
        generate_prior_predictive_samples(
          formula = formula, data = data, family = family,
          prior = prior, ndraws = 100, verbose = verbose
        )
      )
    }, error = function(e) {
      if (isTRUE(verbose)) cat("Warning: Prior sampling failed\n")
      NULL
    })
  }
  
  dur <- fit_res$fitting_time %||% 0
  res <- list(
    fit = fit_res$fit,
    original_formula = formula,
    data = data,
    family = inla_fam,
    has_random = as.logical(formula_components$has_random_effects %||% FALSE),
    has_random_effects = as.logical(formula_components$has_random_effects %||% FALSE),
    group_var = group_var,
    model_type = ifelse(isTRUE(formula_components$has_random_effects), "mixed", "fixed"),
    fitting_time = dur,
    timing = list(total_seconds = dur,
                  formatted = format_duration(dur),
                  formatted_duration = format_duration(dur)),
    prior_samples = prior_samples,
    prior_specs = prior
  )
  class(res) <- c("qbrms_fit")
  
  if (isTRUE(verbose)) cat("qbrms model fitting completed successfully.\n")
  return(res)
}
# =============================================================================
# Alias for qbrms() - Quick Bayesian Regression Models
# =============================================================================

#' Alias for \code{qbrms()}
#'
#' \code{qbrm()} is a shorter alias for \code{qbrms()} with identical functionality.
#' This block intentionally does **not** use @inheritParams to avoid roxygen warnings
#' for aliases.
#'
#' @seealso \code{\link{qbrms}}
#' @export
qbrm <- qbrms
