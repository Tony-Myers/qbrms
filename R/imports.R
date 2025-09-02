# R/imports.R - Final version with all imports

#' @importFrom stats as.formula binomial coef complete.cases dnorm gaussian lm median model.frame model.matrix model.response plogis qnorm quantile rbinom rnorm rpois sd vcov terms reformulate
#' @importFrom utils head modifyList
#' @importFrom graphics abline
#' @importFrom MASS ginv
#' @import ggplot2
NULL

# Handle ggplot2 variable bindings to avoid R CMD check warnings
utils::globalVariables(c(
  "draw", "y", "value", "rep_id", "x", "group", "median", "quantile", 
  "obs", "lower", "upper", ".data", "d"
))

# Define the null coalescing operator used throughout the package
`%||%` <- function(x, y) if (is.null(x)) y else x