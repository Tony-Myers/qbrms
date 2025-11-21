# =============================================================================
# R/metrics.R
# =============================================================================

# small helper
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- INTERNAL: extract pointwise CPO (if present) ---------------------------
.qb_extract_cpo <- function(fit) {
  # INLA objects are usually in fit$fit
  ff <- fit$fit %||% fit
  if (is.null(ff)) return(NULL)
  
  cp <- NULL
  # common INLA shapes:
  if (!is.null(ff$cpo)) {
    if (is.numeric(ff$cpo)) {
      cp <- ff$cpo
    } else if (is.list(ff$cpo) && !is.null(ff$cpo$cpo)) {
      cp <- ff$cpo$cpo
    }
  }
  # any other stash?
  if (is.null(cp) && !is.null(fit$cpo)) {
    if (is.numeric(fit$cpo)) cp <- fit$cpo
    if (is.list(fit$cpo) && !is.null(fit$cpo$cpo)) cp <- fit$cpo$cpo
  }
  if (!is.null(cp)) as.numeric(cp) else NULL
}

# ---- INTERNAL: LOO from CPO -------------------------------------------------
.qb_loo_from_cpo <- function(cpo_vec) {
  # Guard & stabilise zeros
  cpo_vec <- as.numeric(cpo_vec)
  if (length(cpo_vec) == 0L) return(NULL)
  cpo_vec[!is.finite(cpo_vec)] <- NA_real_
  # avoid -Inf from log(0)
  log_cpo <- log(pmax(cpo_vec, .Machine$double.xmin))
  ok <- is.finite(log_cpo)
  n_ok <- sum(ok)
  if (n_ok == 0L) return(NULL)
  
  elpd_loo <- sum(log_cpo[ok])
  looic    <- -2 * elpd_loo
  
  # SE of a sum of i.i.d. terms ~ sqrt(n) * sd(pointwise)
  se_elpd <- if (n_ok >= 2L) sqrt(n_ok) * stats::sd(log_cpo[ok]) else NA_real_
  se_loo  <- if (is.na(se_elpd)) NA_real_ else 2 * se_elpd
  
  list(
    looic        = looic,
    elpd_loo     = elpd_loo,
    se_elpd_loo  = se_elpd,
    se_looic     = se_loo,
    pointwise    = data.frame(elpd_loo_i = log_cpo)
  )
}

# ---- WAIC wrapper (kept simple/robust) --------------------------------------
# Returns $waic and, when available, $se_waic (often NA with INLA summary only)
#' @export
waic <- function(object, ...) {
  # Try to read from the INLA fit if present
  ff <- object$fit %||% object
  if (!is.null(ff$waic) && !is.null(ff$waic$waic)) {
    # INLA provides scalar WAIC; pointwise often not exposed here
    return(list(
      waic    = as.numeric(ff$waic$waic),
      se_waic = ff$waic$se %||% NA_real_
    ))
  }
  
  # Fall back: if we have pointwise elpd_loo from CPO, approximate WAIC ~ -2*elpd_loo
  cpo <- .qb_extract_cpo(object)
  if (!is.null(cpo)) {
    loo_bits <- .qb_loo_from_cpo(cpo)
    if (!is.null(loo_bits)) {
      return(list(
        waic    = -2 * loo_bits$elpd_loo,
        se_waic = if (is.null(loo_bits$se_elpd_loo)) NA_real_ else 2 * loo_bits$se_elpd_loo
      ))
    }
  }
  
  # Last resort: try DIC as a proxy (no SE)
  if (!is.null(ff$dic) && !is.null(ff$dic$dic)) {
    return(list(waic = as.numeric(ff$dic$dic), se_waic = NA_real_))
  }
  
  # Could not compute
  list(waic = NA_real_, se_waic = NA_real_)
}

# ---- LOO wrapper (with SE from CPO if available) ----------------------------
#' @export
loo <- function(object, ...) {
  # qbrms/qbrmO fits
  if (inherits(object, "qbrms_fit") || inherits(object, "qbrms_multinomial_fit") || inherits(object, "qbrmO_fit")) {
    cpo <- .qb_extract_cpo(object)
    if (!is.null(cpo)) {
      out <- .qb_loo_from_cpo(cpo)
      if (!is.null(out)) return(out)
    }
    # If no CPO, degrade gracefully: approximate from WAIC (no SE)
    w <- waic(object)
    if (!is.na(w$waic)) {
      elpd <- -0.5 * w$waic
      return(list(
        looic       = -2 * elpd,
        elpd_loo    = elpd,
        se_elpd_loo = NA_real_,
        se_looic    = NA_real_
      ))
    }
    return(list(looic = NA_real_, elpd_loo = NA_real_, se_elpd_loo = NA_real_, se_looic = NA_real_))
  }
  
  # If user passed a 'loo' object from the loo package, just return it
  if (inherits(object, "loo") || inherits(object, "psis_loo")) {
    return(object)
  }
  
  stop("Unsupported object type for loo().", call. = FALSE)
}

# ---- DIC wrapper (unchanged behavior) ---------------------------------------
#' @export
dic <- function(object, ...) {
  ff <- object$fit %||% object
  if (!is.null(ff$dic) && !is.null(ff$dic$dic)) {
    return(list(dic = as.numeric(ff$dic$dic)))
  }
  list(dic = NA_real_)
}