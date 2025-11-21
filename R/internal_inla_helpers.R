# =============================================================================
# R/internal_inla_helpers.R
# =============================================================================

# Internal helper to locate the embedded INLA result in a qbrms_fit object.
# Adjust this if the internal structure of qbrms_fit changes.
.qbrms_find_inla_result <- function(object) {
  # Common patterns for storing INLA results inside fitted objects.
  if (!is.null(object$fit) && inherits(object$fit, "inla")) {
    return(object$fit)
  }
  
  if (!is.null(object$inla_result) && inherits(object$inla_result, "inla")) {
    return(object$inla_result)
  }
  
  stop("Could not locate embedded INLA result in 'qbrms_fit' object.",
       call. = FALSE)
}