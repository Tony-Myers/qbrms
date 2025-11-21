# =============================================================================
# R/loo_compare.R
# =============================================================================

#' Compare models by LOO (default) or WAIC
#'
#' @description
#' Compare multiple fitted models and rank them by out-of-sample fit.
#' If you pass qbrms/qbrmO fit objects, this uses the package's \code{loo()}
#' / \code{waic()} wrappers under the hood. If you pass actual \code{loo}
#' objects (from the \pkg{loo} package), it will delegate to
#' \code{loo::loo_compare()} automatically.
#'
#' @param ... One or more fitted models (qbrms/qbrmO), or \code{loo} objects; you
#'   can also pass a single named list of models.
#' @param criterion Character, "loo" (default) or "waic".
#' @param sort Logical; if TRUE (default) the best model is first.
#'
#' @return A data.frame with model names, estimate on the ELPD scale
#'   (higher is better), standard error (if available), differences vs best,
#'   and ranks.
#' @export
loo_compare <- function(..., criterion = c("loo", "waic"), sort = TRUE) {
  criterion <- match.arg(criterion)
  
  # Collect models
  mods <- list(...)
  if (length(mods) == 1L && is.list(mods[[1L]]) && !inherits(mods[[1L]], "qbrms_fit")) {
    mods <- mods[[1L]]
  }
  if (length(mods) < 2L) stop("Provide at least two models.", call. = FALSE)
  
  # Names
  nm <- names(mods)
  if (is.null(nm)) nm <- paste0("model", seq_along(mods))
  names(mods) <- nm
  
  # If all inputs look like 'loo' objects, delegate
  all_loo_like <- all(vapply(mods, inherits, logical(1), what = c("loo", "psis_loo")))
  if (all_loo_like) {
    if (requireNamespace("loo", quietly = TRUE)) {
      return(loo::loo_compare(mods))
    } else {
      stop("Inputs are 'loo' objects but the 'loo' package is not available.", call. = FALSE)
    }
  }
  
  # Otherwise assume qbrms/qbrmO fits and extract metrics
  elpd <- numeric(length(mods))
  se_elpd <- rep(NA_real_, length(mods))
  
  for (i in seq_along(mods)) {
    obj <- mods[[i]]
    is_qbrms_like <- inherits(obj, "qbrms_fit") || inherits(obj, "qbrms_multinomial_fit") ||
      inherits(obj, "qbrmO_fit")
    if (!is_qbrms_like) {
      stop("Object '", nm[i], "' is not a qbrms/qbrmO fit (and not a 'loo' object).", call. = FALSE)
    }
    
    if (criterion == "loo") {
      li <- loo(obj)
      if (!is.null(li$elpd_loo)) {
        elpd[i] <- li$elpd_loo
        se_elpd[i] <- li$se_loo_elpd %||% li$se_elpd_loo %||% li$elpd_loo_se %||% NA_real_
      } else if (!is.null(li$looic)) {
        elpd[i] <- -0.5 * li$looic
        se_elpd[i] <- (li$se_looic %||% NA_real_) * 0.5
      } else {
        stop("`loo()` did not return elpd_loo or looic for model '", nm[i], "'.", call. = FALSE)
      }
    } else {
      wi <- waic(obj)
      if (!is.null(wi$waic)) {
        elpd[i] <- -0.5 * wi$waic
        se_elpd[i] <- (wi$se_waic %||% wi$waic_se %||% NA_real_) * 0.5
      } else {
        stop("`waic()` did not return $waic for model '", nm[i], "'.", call. = FALSE)
      }
    }
  }
  
  # Order (higher ELPD is better)
  ord <- if (isTRUE(sort)) order(elpd, decreasing = TRUE) else seq_along(elpd)
  
  elpd_sorted  <- elpd[ord]
  se_sorted    <- se_elpd[ord]
  names_sorted <- nm[ord]
  
  # Differences vs best
  elpd_diff <- elpd_sorted - elpd_sorted[1L]
  se_diff <- rep(NA_real_, length(elpd_sorted))
  if (!is.na(se_sorted[1L])) {
    for (j in seq_along(se_sorted)) {
      if (!is.na(se_sorted[j])) {
        se_diff[j] <- sqrt(se_sorted[j]^2 + se_sorted[1L]^2)
      }
    }
  }
  
  out <- data.frame(
    model     = names_sorted,
    estimate  = elpd_sorted,
    se        = se_sorted,
    elpd_diff = elpd_diff,
    se_diff   = se_diff,
    rank      = seq_along(elpd_sorted),
    row.names = NULL,
    check.names = FALSE
  )
  class(out) <- c("qbrms_loo_compare", "data.frame")
  attr(out, "criterion") <- criterion
  out
}

#' @export
print.qbrms_loo_compare <- function(x, ...) {
  crit <- attr(x, "criterion") %||% "loo"
  cat(sprintf("Model comparison (%s; higher is better)\n\n", crit))
  
  # Drop special class to avoid print recursion
  df <- as.data.frame(x)
  attr(df, "criterion") <- NULL
  
  df$estimate  <- round(df$estimate, 3)
  df$se        <- ifelse(is.na(df$se), NA, round(df$se, 3))
  df$elpd_diff <- round(df$elpd_diff, 3)
  df$se_diff   <- ifelse(is.na(df$se_diff), NA, round(df$se_diff, 3))
  
  print(df, row.names = FALSE)
  invisible(x)
}

# Null-coalescing helper
`%||%` <- function(x, y) if (is.null(x)) y else x