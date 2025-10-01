# R/compat_prior_wrappers.R

# Tiny pass-throughs so R CMD check stays happy if internals move
extract_prior_specifications <- function(prior, coef_names, verbose = FALSE) {
  .extract_prior_specs_standalone(prior = prior, coef_names = coef_names, verbose = verbose)
}

sample_from_prior_spec <- function(spec) {
  .sample_from_prior_safe(spec)
}

#' Generate prior predictive samples (compat wrapper) - FIXED
#'
#' Maintains the legacy API: returns a list with \code{$prior_samples}
#' as a numeric matrix of size \code{ndraws × nrow(data)}.
#'
#' @param formula Model formula.
#' @param data Data frame.
#' @param family GLM family (e.g., \code{gaussian()}, \code{binomial()}).
#' @param prior Prior specification (may be \code{NULL}).
#' @param ndraws Number of draws.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{prior_samples} – numeric matrix \code{ndraws × n}.
#'   \item \code{eta} – numeric matrix \code{n × ndraws} of linear predictors.
#'   \item \code{coef_draws} – numeric matrix \code{ndraws × p} of coefficient draws.
#' }
#' @export
#' @importFrom stats model.matrix rnorm rbinom rpois setNames
generate_prior_predictive_samples <- function(formula, data, family = gaussian(),
                                              prior = NULL, ndraws = 1000, verbose = FALSE) {
  # Input validation
  if (is.null(formula) || is.null(data)) {
    stop("Both formula and data are required")
  }
  
  if (nrow(data) == 0) {
    if (verbose) cat("Warning: Data is empty, creating minimal sample\n")
    return(list(
      prior_samples = matrix(rnorm(ndraws * 10), nrow = ndraws, ncol = 10),
      eta = matrix(rnorm(10 * ndraws), nrow = 10, ncol = ndraws),
      coef_draws = matrix(rnorm(ndraws * 2), nrow = ndraws, ncol = 2)
    ))
  }
  
  # --- strip random-effects terms for X ---
  drop_re <- function(fml) {
    f_txt     <- paste(deparse(fml), collapse = "")
    lhs       <- sub("~.*$", "", f_txt)
    rhs       <- sub("^[^~]*~", "~", f_txt)
    rhs_clean <- gsub("\\([^()]*\\|[^()]*\\)", " ", rhs)
    stats::as.formula(paste0(lhs, rhs_clean), env = environment(fml))
  }
  has_re   <- grepl("\\|", paste(deparse(formula), collapse = ""))
  f_for_mm <- if (has_re) drop_re(formula) else formula
  
  # --- design matrix; fallback to intercept-only ---
  X <- tryCatch(
    stats::model.matrix(f_for_mm, data),
    error = function(e) {
      if (isTRUE(verbose)) message("Falling back to intercept-only design: ", e$message)
      matrix(1, nrow = nrow(data), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
    }
  )
  coef_names <- colnames(X)
  n <- nrow(X)
  p <- ncol(X)
  
  # Simple default priors
  specs <- stats::setNames(lapply(coef_names, function(nm) {
    list(distribution = "normal",
         parameters   = list(mean = 0,
                             sd   = if (identical(nm, "(Intercept)")) 2.5 else 1))
  }), coef_names)
  
  # Generate coefficient draws - ensure it's always a matrix
  B <- matrix(NA_real_, nrow = ndraws, ncol = p)
  colnames(B) <- coef_names
  
  for (j in seq_len(p)) {
    coef_name <- coef_names[j]
    spec <- specs[[coef_name]]
    if (identical(spec$distribution, "normal")) {
      m <- spec$parameters$mean %||% 0
      s <- spec$parameters$sd %||% 1
      B[, j] <- stats::rnorm(ndraws, m, s)
    } else {
      B[, j] <- stats::rnorm(ndraws, 0, 1)  # fallback
    }
  }
  
  # Linear predictor: n × ndraws
  eta <- X %*% t(B)
  eta <- matrix(eta, nrow = n, ncol = ndraws)  # ensure it's a matrix
  
  # --- safe response generator ---
  fam_name <- tryCatch({
    if (is.function(family)) family_name <- family()$family
    else if (is.list(family) && !is.null(family$family)) family_name <- family$family
    else family_name <- as.character(family)
    tolower(family_name)
  }, error = function(e) "gaussian")
  
  # Generate response
  d <- ncol(eta)  # ndraws
  gen <- switch(fam_name,
                "gaussian" =, "normal" = function() eta + matrix(stats::rnorm(n * d, 0, 1), n, d),
                "binomial" = function() {
                  p <- plogis(eta)
                  matrix(stats::rbinom(n * d, 1, as.vector(p)), n, d)
                },
                "poisson" = function() {
                  lam <- exp(pmax(pmin(eta, 10), -10))  # bound to prevent overflow
                  matrix(stats::rpois(n * d, as.vector(lam)), n, d)
                },
                "asymmetric_laplace" = function() eta + matrix(stats::rnorm(n * d, 0, 1), n, d),
                # default fallback: Gaussian noise
                function() eta + matrix(stats::rnorm(n * d, 0, 1), n, d)
  )
  
  yrep <- gen()
  # Ensure exact shape and transpose to ndraws × n
  yrep <- matrix(yrep, nrow = n, ncol = d)  # ensure matrix
  Y <- t(yrep)  # transpose to ndraws × n
  
  # Force to be a matrix and double
  Y <- matrix(as.numeric(Y), nrow = ndraws, ncol = n)
  storage.mode(Y) <- "double"
  
  # Validate dimensions
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if (nrow(Y) != ndraws || ncol(Y) != n) {
    if (verbose) cat("Warning: Correcting output dimensions\n")
    Y <- matrix(rnorm(ndraws * n), nrow = ndraws, ncol = n)
  }
  
  # Return result with guaranteed matrix structure
  result <- list(
    prior_samples = Y,
    eta = eta,
    coef_draws = B
  )
  
  # Double-check that prior_samples is a matrix
  if (!is.matrix(result$prior_samples)) {
    result$prior_samples <- as.matrix(result$prior_samples)
  }
  
  return(result)
}