# =============================================================================
# R/compat_emmeans.R
# Compatibility wrapper so users can call emmeans(fit, ...) on qbrms_fit.
# If the object is not a qbrms_fit, we defer to the external emmeans package
# when it is available. Otherwise we give a clear error message.
# =============================================================================

#' Estimated marginal means (compatibility wrapper)
#'
#' This wrapper lets you call `emmeans()` on a `qbrms_fit` without attaching
#' the external \pkg{emmeans} package. For non-\code{qbrms_fit} objects, it
#' forwards to \pkg{emmeans} if that package is installed.
#'
#' @param object A model object; if it is a \code{qbrms_fit} we dispatch to
#'   \code{qbrms_emmeans()}.
#' @param specs Term(s) for which to compute estimated marginal means. For
#'   \code{qbrms_fit}, this is passed to \code{qbrms_emmeans()} unchanged.
#' @param ... Additional arguments forwarded either to \code{qbrms_emmeans()}
#'   or to \code{emmeans::emmeans()} as appropriate.
#'
#' @return A data frame for \code{qbrms_fit}; otherwise whatever
#'   \code{emmeans::emmeans()} returns.
#' @export
emmeans <- function(object, specs, ...) {
  # Route qbrms models to the in-package implementation
  if (inherits(object, "qbrms_fit")) {
    return(qbrms_emmeans(object, specs = specs, ...))
  }
  
  # Otherwise, forward to the external emmeans if it is available
  if (requireNamespace("emmeans", quietly = TRUE)) {
    fun <- utils::getFromNamespace("emmeans", "emmeans")
    return(fun(object, specs = specs, ...))
  }
  
  stop("The 'emmeans' package is not installed and the object is not a 'qbrms_fit'.")
}
