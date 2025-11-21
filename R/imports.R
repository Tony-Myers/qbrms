# Centralised imports and globals for qbrms

# Packages imported into the NAMESPACE
#' @import future
#' @import future.apply
#' @import loo
#' @import miniUI
#' @import shiny
#' @importFrom stats delete.response model.matrix terms as.formula na.omit
NULL

# Handle ggplot2 / NSE variable bindings and other globals to avoid R CMD check notes
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "draw", "y", "value", "rep_id", "x", "group", "median", "quantile",
    "obs", "lower", "upper", ".data", "d",
    "Model", "Type", "Value", "effect1__", "estimate__", "id__", "lower__", "upper__",
    "type", "ps", "Parameter", "sig_level",
    "CI_low", "CI_high", "Dist",
    "Group", "Success_Rate", "Group_Size", "Lower", "Upper", "groups",
    "is_pseudo", "group_id", "regularisation_strength", "prior_random_intercept",
    "sn", "nbinomial"
  ))
}

# TMB imports for ordinal regression (kept as a comment for clarity)

# Null-coalescing operator used throughout the package
`%||%` <- function(x, y) if (is.null(x)) y else x

#' @importFrom stats coef complete.cases dnorm model.frame qnorm reshape vcov
#' @importFrom utils capture.output modifyList
NULL