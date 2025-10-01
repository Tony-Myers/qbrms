#' Ordinal regression via binary decomposition (fallback)
#'
#' Splits an ordered response with K levels into K-1 binary problems
#' (thresholds y > c_j) and fits a simple binomial GLM for each split.
#'
#' @param formula Model formula with an ordered response on the LHS.
#' @param data Data frame.
#' @param verbose Logical; print progress messages.
#' @param ... Ignored (compat).
#'
#' @return An object of class \code{c("ordinal_binary_qbrms_fit","qbrms_fit")} with:
#' \itemize{
#'   \item \code{binary_models}: list of length K-1 of fitted \code{glm} objects
#'   \item \code{thresholds}: character vector of thresholds used
#'   \item \code{response}: response variable name
#'   \item \code{levels}: factor levels of the ordered response
#'   \item \code{ordinal_levels}: factor levels (for test compatibility)
#' }
#' @export
qbrms_ordinal_binary <- function(formula, data, verbose = FALSE, ...) {
  drop_re <- function(fml) {
    f_txt    <- paste(deparse(fml), collapse = "")
    lhs      <- sub("~.*$", "", f_txt)
    rhs      <- sub("^[^~]*~", "~", f_txt)
    rhs_clean <- gsub("\\([^()]*\\|[^()]*\\)", " ", rhs)
    stats::as.formula(paste0(lhs, rhs_clean), env = environment(fml))
  }
  
  resp <- all.vars(formula)[1L]
  y <- data[[resp]]
  if (!is.ordered(y)) y <- ordered(y)
  K <- nlevels(y)
  if (K < 3L) stop("qbrms_ordinal_binary() needs an ordered response with at least 3 levels.")
  
  levs <- levels(y)
  thr  <- levs[-length(levs)]     # K-1 thresholds
  fe_form <- if (grepl("\\|", paste(deparse(formula), collapse = ""))) drop_re(formula) else formula
  
  binmods <- vector("list", length = K - 1L)
  names(binmods) <- paste0("cut_", seq_len(K - 1L))
  
  for (j in seq_len(K - 1L)) {
    zj <- as.integer(y > thr[j])
    df <- data
    df[[resp]] <- zj
    rhs_terms <- paste(attr(stats::terms(fe_form), "term.labels"), collapse = " + ")
    glm_form  <- if (nzchar(rhs_terms)) {
      stats::as.formula(paste(resp, "~", rhs_terms), env = environment(fe_form))
    } else {
      stats::as.formula(paste(resp, "~ 1"), env = environment(fe_form))
    }
    binmods[[j]] <- stats::glm(glm_form, data = df, family = stats::binomial())
  }
  
  out <- list(
    call           = match.call(),
    binary_models  = binmods,        # Use binmods, not binary_models
    thresholds     = thr,            # Use thr, not thresholds
    levels         = levs,           # Factor levels
    ordinal_levels = levs,           # Add this for test compatibility
    response       = resp,           # Response variable name
    formula        = formula,
    data           = data,
    model_type     = "ordinal_binary"
  )
  
  class(out) <- c("ordinal_binary_qbrms_fit", "qbrms_fit")
  out
}

# Also update the fit_ordinal_model function to ensure ordinal_levels is included:

#' Fit ordinal model with fallbacks (UPDATED)
#' @keywords internal
fit_ordinal_model <- function(formula, data, inla_family, control.compute, verbose = TRUE, ...) {
  if (verbose) cat("Fitting ordinal model...\n")
  
  # Extract response and check levels
  response_var <- all.vars(formula)[1]
  y <- data[[response_var]]
  
  if (!is.factor(y)) {
    y <- as.factor(y)
    data[[response_var]] <- y
  }
  
  if (!is.ordered(y)) {
    y <- ordered(y)
    data[[response_var]] <- y
  }
  
  if (nlevels(y) < 3) {
    stop("Need at least 3 ordinal levels")
  }
  
  ordinal_levels <- levels(y)  # Store for return
  
  # Try ordinal package first
  if (requireNamespace("ordinal", quietly = TRUE)) {
    model <- tryCatch({
      ordinal::clm(formula, data = data)
    }, error = function(e) {
      if (verbose) cat("ordinal::clm failed:", e$message, "\n")
      NULL
    })
    
    if (!is.null(model)) {
      coefs <- stats::coef(model)
      se <- sqrt(diag(stats::vcov(model)))
      
      summary_fixed <- data.frame(
        mean = coefs,
        sd = se,
        row.names = names(coefs)
      )
      
      result <- list(
        fit = list(
          summary.fixed = summary_fixed,
          converged = model$convergence == 0,
          fallback_type = "ordinal"
        ),
        original_formula = formula,
        data = data,
        family = inla_family,
        model_type = "ordinal",
        ordinal_levels = ordinal_levels,  # Add this for tests
        fitting_time = 0
      )
      
      class(result) <- c("ordinal_qbrms_fit", "qbrms_fit")
      return(result)
    }
  }
  
  # Fallback to qbrms_ordinal_binary
  if (verbose) cat("Using binary decomposition fallback...\n")
  
  result <- qbrms_ordinal_binary(formula, data, verbose = verbose)
  result$family <- inla_family
  result$fitting_time <- 0
  
  return(result)
}