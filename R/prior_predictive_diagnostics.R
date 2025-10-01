# =============================================================================
# R/prior_predictive_diagnostics.R
# =============================================================================

#' Prior predictive diagnostics and sensibility report
#'
#' Summarise prior predictive draws to check basic support, scale and shape,
#' and (optionally) how simple statistics of the observed data compare with
#' the prior-predictive distribution. Returns an object with a concise verdict.
#'
#' @param object A qbrms prior object: qbrms_prior_fit, qbrms_prior_only, or a qbrms_fit
#'   that contains \code{prior_samples}.
#' @param level Credible level for central intervals (default 0.95). Reserved.
#' @param support Optional override of the implied support: one of "real",
#'   "positive", "proportion", or "bounded". If \code{NULL}, an attempt is made
#'   to infer from the family.
#' @param lower,upper Optional numeric bounds used when \code{support = "bounded"}.
#'   If \code{support = "proportion"}, the default is \code{c(0, 1)}.
#' @param trials Optional integer vector for binomial data (bounds helper).
#' @param plausible_lower,plausible_upper Optional numeric bounds defining a user-declared
#'   “plausible range” for the outcome on the response scale. When both are supplied,
#'   the function reports the fraction of prior-predictive mass that falls in
#'   \code{[plausible_lower, plausible_upper]} and incorporates this into the verdict.
#' @param include_observed Logical; if \code{TRUE} and the object contains data,
#'   the report compares simple statistics of \code{y} to their prior-predictive
#'   reference distributions.
#' @param seed Optional seed for reproducibility.
#' @return An object of class \code{qbrms_prior_diagnostics}.
#' @export
prior_pp_diagnostics <- function(object,
                                 level = 0.95,
                                 support = NULL,
                                 lower = NULL,
                                 upper = NULL,
                                 trials = NULL,
                                 plausible_lower = NULL,
                                 plausible_upper = NULL,
                                 include_observed = TRUE,
                                 seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  fam_name <- tryCatch(extract_family_name(object$family), error = function(e) "unknown")
  
  # Obtain yrep
  yrep <- NULL
  if (!is.null(object$prior_samples) && is.matrix(object$prior_samples)) {
    yrep <- object$prior_samples
  } else if (inherits(object, "qbrms_fit") && !is.null(object$prior_samples)) {
    yrep <- object$prior_samples
  }
  if (is.null(yrep) || !is.matrix(yrep) || nrow(yrep) < 1L || ncol(yrep) < 1L) {
    stop("No prior predictive samples found in 'object'.")
  }
  
  # Extract y (if present) with alignment mask
  y <- NULL
  idx <- NULL
  response_var <- tryCatch(all.vars(object$original_formula)[1], error = function(e) NULL)
  if (!is.null(response_var) && !is.null(object$data) && response_var %in% names(object$data)) {
    y_raw <- object$data[[response_var]]
    y_num <- if (is.ordered(y_raw)) {
      as.numeric(y_raw)
    } else if (is.factor(y_raw)) {
      suppressWarnings(as.numeric(as.character(y_raw)))
    } else {
      as.numeric(y_raw)
    }
    idx <- is.finite(y_num)
    if (any(idx)) y <- y_num[idx]
  }
  if (!is.null(idx) && is.matrix(yrep) && ncol(yrep) == length(idx)) {
    yrep <- yrep[, idx, drop = FALSE]
  }
  
  # Infer default support
  if (is.null(support)) {
    support <- switch(fam_name,
                      "binomial"   = "bounded",
                      "beta"       = "proportion",
                      "poisson"    = "positive",
                      "nbinomial"  = "positive",
                      "gamma"      = "positive",
                      "lognormal"  = "positive",
                      "real"
    )
  }
  if (identical(support, "proportion")) {
    lower <- if (is.null(lower)) 0 else lower
    upper <- if (is.null(upper)) 1 else upper
  }
  if (identical(support, "bounded") && is.null(lower) && is.null(upper) && fam_name == "binomial") {
    if (!is.null(trials) && length(trials) == ncol(yrep)) {
      lower <- 0; upper <- max(trials, na.rm = TRUE)
    } else {
      vals <- as.numeric(yrep[1L, ])
      if (all(vals %in% c(0, 1))) { lower <- 0; upper <- 1 }
    }
  }
  
  n_draws <- nrow(yrep)
  n_obs   <- ncol(yrep)
  
  # Pooled quantiles of yrep
  probs <- c(0.005, 0.025, 0.05, 0.5, 0.95, 0.975, 0.995)
  marginal_quantiles <- as.numeric(stats::quantile(as.numeric(yrep), probs = probs, na.rm = TRUE))
  names(marginal_quantiles) <- paste0(formatC(probs * 100, format = "f", digits = 1), "%")
  
  # Draw-wise summaries
  draw_means <- rowMeans(yrep)
  draw_sds   <- apply(yrep, 1L, stats::sd)
  draw_mins  <- apply(yrep, 1L, min)
  draw_maxs  <- apply(yrep, 1L, max)
  draw_p0    <- apply(yrep == 0, 1L, mean)
  
  qfun <- function(x) stats::quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  by_draw_summaries <- rbind(
    mean = qfun(draw_means),
    sd   = qfun(draw_sds),
    min  = qfun(draw_mins),
    max  = qfun(draw_maxs),
    prop_zero = qfun(draw_p0)
  )
  by_draw_summaries <- as.data.frame(by_draw_summaries)
  by_draw_summaries$statistic <- rownames(by_draw_summaries)
  rownames(by_draw_summaries) <- NULL
  names(by_draw_summaries)[1:3] <- c("q2.5", "q50", "q97.5")
  
  # Support violations
  prop_neg <- mean(yrep < 0, na.rm = TRUE)
  prop_gt1 <- mean(yrep > 1, na.rm = TRUE)
  prop_oob <- NA_real_
  if (!is.null(lower) || !is.null(upper)) {
    lo <- if (is.null(lower)) -Inf else lower
    hi <- if (is.null(upper))  Inf else upper
    prop_oob <- mean(yrep < lo | yrep > hi, na.rm = TRUE)
  }
  support_violations <- list(
    proportion_negative      = prop_neg,
    proportion_above_one     = prop_gt1,
    proportion_out_of_bounds = prop_oob
  )
  
  # Plausible range coverage (if provided)
  plausible_coverage <- NA_real_
  if (!is.null(plausible_lower) && !is.null(plausible_upper)) {
    lo_pl <- min(plausible_lower, plausible_upper)
    hi_pl <- max(plausible_lower, plausible_upper)
    plausible_coverage <- mean(yrep >= lo_pl & yrep <= hi_pl, na.rm = TRUE)
  }
  
  # Observed comparison
  observed_stats <- NULL
  prior_ref_dists <- NULL
  pp_values <- NULL
  diag_plot <- NULL
  
  if (isTRUE(include_observed) && !is.null(y)) {
    stat_obs <- c(
      mean = mean(y), sd = stats::sd(y), min = min(y), max = max(y),
      prop_zero = mean(y == 0)
    )
    ref_dists <- list(
      mean = draw_means,
      sd   = draw_sds,
      min  = draw_mins,
      max  = draw_maxs,
      prop_zero = draw_p0
    )
    tail_prob <- function(draw_dist, obs) {
      p <- mean(draw_dist <= obs)
      p <- min(max(p, 0), 1)
      2 * min(p, 1 - p)
    }
    pvals <- vapply(names(ref_dists), function(nm) tail_prob(ref_dists[[nm]], stat_obs[[nm]]), numeric(1L))
    
    observed_stats  <- stat_obs
    prior_ref_dists <- lapply(ref_dists, function(x) stats::quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
    pp_values       <- pvals
    
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      dfm <- data.frame(value = draw_means)
      diag_plot <- ggplot2::ggplot(dfm, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(bins = 30, fill = "#b3cde0", colour = "#6f93ad", linewidth = 0.25, alpha = 0.9) +
        ggplot2::geom_vline(xintercept = stat_obs[["mean"]], colour = "#1a1a1a", linewidth = 0.8, linetype = "dashed") +
        ggplot2::labs(title = "Prior predictive: distribution of draw-wise means",
                      x = "mean(yrep)", y = "Count") +
        ggplot2::theme_minimal()
    }
  }
  
  # Verdict rules (now include plausible coverage if supplied)
  verdict <- .prior_diag_verdict(
    support, lower, upper,
    support_violations,
    include_observed, pp_values,
    plausible_coverage
  )
  
  res <- list(
    sample_sizes        = list(ndraws = nrow(yrep), nobs = ncol(yrep)),
    support_used        = list(support = support, lower = lower, upper = upper, family = fam_name),
    marginal_quantiles  = marginal_quantiles,
    by_draw_summaries   = by_draw_summaries,
    support_violations  = support_violations,
    plausible_bounds    = c(lower = plausible_lower, upper = plausible_upper),
    plausible_coverage  = plausible_coverage,
    observed_stats      = observed_stats,
    prior_ref_dists     = prior_ref_dists,
    pp_values           = pp_values,
    plot                = diag_plot,
    verdict             = verdict
  )
  class(res) <- "qbrms_prior_diagnostics"
  res
}

# Internal: verdict logic (incorporates plausible coverage if provided)
.prior_diag_verdict <- function(support, lower, upper, sv,
                                include_observed, pvals,
                                plausible_coverage) {
  status  <- "OK"
  reasons <- character(0)
  
  # Thresholds
  t_rev <- 0.05   # 5% out-of-support is severe
  t_cau <- 0.005  # 0.5% out-of-support is a caution
  p_rev <- 0.01   # very small two-sided p
  p_cau <- 0.05
  c_rev <- 0.10   # plausible range coverage < 10% is severe
  c_cau <- 0.50   # plausible range coverage < 50% is a caution
  
  # Support-based checks
  if (identical(support, "positive") && is.finite(sv$proportion_negative)) {
    if (sv$proportion_negative > t_rev) { status <- "Revise"; reasons <- c(reasons, "non-positive mass under a positive-only outcome") }
    else if (sv$proportion_negative > t_cau && status != "Revise") { status <- "Caution"; reasons <- c(reasons, "trace non-positive mass under a positive-only outcome") }
  }
  if (identical(support, "proportion") && is.finite(sv$proportion_above_one)) {
    if (sv$proportion_above_one > t_rev) { status <- "Revise"; reasons <- c(reasons, "mass above 1 for a proportion") }
    else if (sv$proportion_above_one > t_cau && status != "Revise") { status <- "Caution"; reasons <- c(reasons, "trace mass above 1 for a proportion") }
  }
  if (!is.null(lower) || !is.null(upper)) {
    if (is.finite(sv$proportion_out_of_bounds)) {
      if (sv$proportion_out_of_bounds > t_rev) { status <- "Revise"; reasons <- c(reasons, "substantial mass outside stated bounds") }
      else if (sv$proportion_out_of_bounds > t_cau && status != "Revise") { status <- "Caution"; reasons <- c(reasons, "trace mass outside stated bounds") }
    }
  }
  
  # Plausible range coverage (optional)
  if (is.finite(plausible_coverage)) {
    if (plausible_coverage < c_rev) { status <- "Revise"; reasons <- c(reasons, "very low coverage of user-declared plausible range") }
    else if (plausible_coverage < c_cau && status != "Revise") { status <- "Caution"; reasons <- c(reasons, "low coverage of user-declared plausible range") }
  }
  
  # Observed comparison (optional)
  if (isTRUE(include_observed) && !is.null(pvals)) {
    minp <- min(pvals, na.rm = TRUE)
    if (is.finite(minp)) {
      if (minp < p_rev) { status <- "Revise"; reasons <- c(reasons, "observed summary lies in extreme prior-predictive tail") }
      else if (minp < p_cau && status != "Revise") { status <- "Caution"; reasons <- c(reasons, "observed summary is somewhat surprising under the prior") }
    }
  }
  
  list(
    status = status,
    reasons = unique(reasons),
    thresholds = list(t_rev = t_rev, t_cau = t_cau, p_rev = p_rev, p_cau = p_cau, c_rev = c_rev, c_cau = c_cau)
  )
}

#' A convenience wrapper mirroring pp_check's show_observed flag
#'
#' @param object A qbrms prior object.
#' @param show_observed Logical; compare observed summaries when available.
#' @param plausible_lower,plausible_upper Optional plausible range bounds to score coverage.
#' @param ... Passed to \code{prior_pp_diagnostics()}.
#' @return The diagnostics object, invisibly, after printing a summary.
#' @export
prior_pp_summary <- function(object,
                             show_observed = FALSE,
                             plausible_lower = NULL,
                             plausible_upper = NULL,
                             ...) {
  x <- prior_pp_diagnostics(
    object,
    include_observed = isTRUE(show_observed),
    plausible_lower  = plausible_lower,
    plausible_upper  = plausible_upper,
    ...
  )
  print(x)
  invisible(x)
}

#' Print method for qbrms_prior_diagnostics objects
#' @param x A qbrms_prior_diagnostics object
#' @param digits Number of decimal places to display (default 3)  # <-- ADD THIS LINE
#' @param ... Additional arguments (unused)
#' @return The input object, returned invisibly
#' @method print qbrms_prior_diagnostics
#' @export
print.qbrms_prior_diagnostics <- function(x, digits = 3, ...) {
  ss <- x$sample_sizes
  su <- x$support_used
  sv <- x$support_violations
  
  cat("Prior predictive diagnostics\n")
  cat(sprintf("  Draws: %d   Observations per draw: %d\n", ss$ndraws, ss$nobs))
  cat(sprintf("  Family: %s\n", su$family %||% "unknown"))
  cat(sprintf("  Support: %s", su$support))
  if (!is.null(su$lower) || !is.null(su$upper)) {
    cat(sprintf("   Bounds: [%s, %s]",
                if (is.null(su$lower)) "-Inf" else format(su$lower),
                if (is.null(su$upper)) "Inf"  else format(su$upper)))
  }
  cat("\n")
  
  v <- x$verdict
  cat(sprintf("\nVerdict: %s", v$status))
  if (length(v$reasons)) cat(sprintf("  (%s)", paste(v$reasons, collapse = "; ")))
  cat("\n")
  
  cat("\nPooled prior-predictive quantiles of yrep:\n")
  q <- x$marginal_quantiles
  qstr <- paste(sprintf("%s: %s", names(q), format(round(q, digits), nsmall = digits)), collapse = " | ")
  cat("  ", qstr, "\n", sep = "")
  
  if (identical(su$support, "real")) {
    cat("\nSupport checks: skipped (outcome has full real-line support)\n")
  } else {
    cat("\nSupport checks (proportion of values):\n")
    cat(sprintf("  negative = %s   above_one = %s   out_of_bounds = %s\n",
                format(round(sv$proportion_negative, digits), nsmall = digits),
                format(round(sv$proportion_above_one, digits), nsmall = digits),
                if (is.na(sv$proportion_out_of_bounds)) "NA" else
                  format(round(sv$proportion_out_of_bounds, digits), nsmall = digits)))
  }
  
  # Plausible range coverage, if provided
  if (all(!is.null(x$plausible_bounds))) {
    pb <- x$plausible_bounds
    cov_str <- if (is.finite(x$plausible_coverage)) format(round(x$plausible_coverage, 3), nsmall = 3) else "NA"
    cat(sprintf("\nCoverage of plausible range [%s, %s]: %s\n",
                as.character(pb[["lower"]]), as.character(pb[["upper"]]), cov_str))
  }
  
  cat("\nDraw-wise summaries (quantiles across draws):\n")
  df <- x$by_draw_summaries
  for (i in seq_len(nrow(df))) {
    cat(sprintf("  %-10s  q2.5=%s  q50=%s  q97.5=%s\n",
                df$statistic[i],
                format(round(df$`q2.5`[i],  digits), nsmall = digits),
                format(round(df$`q50`[i],   digits), nsmall = digits),
                format(round(df$`q97.5`[i], digits), nsmall = digits)))
  }
  
  if (!is.null(x$observed_stats)) {
    cat("\nObserved vs prior-predictive reference (two-sided tail probabilities):\n")
    for (nm in names(x$observed_stats)) {
      obs <- x$observed_stats[[nm]]
      pv  <- x$pp_values[[nm]]
      cat(sprintf("  %-10s  observed = %s   p(two-sided) = %s\n",
                  nm,
                  format(round(obs, digits), nsmall = digits),
                  format(round(pv, 3), nsmall = 3)))
    }
    if (!is.null(x$plot)) {
      cat("  Use $plot for a quick visual of mean(yrep) vs observed mean.\n")
    }
  }
  invisible(x)
}

# A tiny infix helper to match the rest of your codebase
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
