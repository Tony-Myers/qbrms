# =============================================================================
# R/ic_methods.R
# =============================================================================

#' Model comparison criteria for qbrms models
#'
#' Compute approximate DIC, LOO and WAIC for qbrms model fits.
#'
#' These functions provide generic interfaces (`dic()`, `loo()`, `waic()`)
#' and S3 methods for `qbrms_fit` objects that extract the corresponding
#' criteria from the underlying INLA fit where available.
#'
#' @param object A \code{qbrms_fit} object.
#' @param ... Additional arguments passed to internal methods or underlying tools.
#'
#' @return For \code{dic()}, \code{loo()} and \code{waic()} methods on
#'   \code{qbrms_fit} objects, a list containing the corresponding criterion
#'   (for example, \code{list(dic = ...)}, \code{list(looic = ..., elpd_loo = ...)},
#'   or \code{list(waic = ...)}). If the criterion cannot be computed, \code{NA_real_}
#'   is returned.
#'
#' @name qbrms-model-criteria
NULL

# -------------------------------------------------------------------------
# Public generics
# -------------------------------------------------------------------------

#' @rdname qbrms-model-criteria
#' @export
waic <- function(object, ...) UseMethod("waic")

#' @rdname qbrms-model-criteria
#' @export
loo <- function(object, ...) UseMethod("loo")

#' @rdname qbrms-model-criteria
#' @export
dic <- function(object, ...) UseMethod("dic")

# -------------------------------------------------------------------------
# S3 methods for qbrms_fit
# -------------------------------------------------------------------------

#' @rdname qbrms-model-criteria
#' @method waic qbrms_fit
#' @export
waic.qbrms_fit <- function(object, ...) {
  ir <- .qbrms_find_inla_result(object)
  # INLA style: result$waic$waic (numeric)
  if (is.list(ir) && !is.null(ir$waic)) {
    if (is.list(ir$waic) && !is.null(ir$waic$waic)) {
      return(list(waic = as.numeric(ir$waic$waic)))
    }
    if (is.numeric(ir$waic)) {
      return(list(waic = as.numeric(ir$waic[1])))
    }
  }
  NA_real_
}

#' @rdname qbrms-model-criteria
#' @method loo qbrms_fit
#' @export
loo.qbrms_fit <- function(object, ...) {
  ir <- .qbrms_find_inla_result(object)
  # CPO-based LOO: elpd_loo = sum(log(cpo_i)), LOOIC = -2 * elpd_loo
  if (is.list(ir) && !is.null(ir$cpo) && !is.null(ir$cpo$cpo)) {
    cpo <- ir$cpo$cpo
    ok  <- is.finite(cpo) & (cpo > 0)
    if (any(ok)) {
      elpd_loo <- sum(log(pmax(cpo[ok], .Machine$double.eps)))
      looic    <- -2 * elpd_loo
      return(list(
        looic    = as.numeric(looic),
        elpd_loo = as.numeric(elpd_loo)
      ))
    }
  }
  NA_real_
}

#' @rdname qbrms-model-criteria
#' @method dic qbrms_fit
#' @export
dic.qbrms_fit <- function(object, ...) {
  ir <- .qbrms_find_inla_result(object)
  if (is.list(ir) && !is.null(ir$dic)) {
    if (is.list(ir$dic) && !is.null(ir$dic$dic)) {
      return(list(dic = as.numeric(ir$dic$dic)))
    }
    if (is.numeric(ir$dic)) {
      return(list(dic = as.numeric(ir$dic[1])))
    }
  }
  NA_real_
}
