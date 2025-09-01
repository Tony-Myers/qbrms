# R/imports.R

#' @keywords internal
#' @importFrom utils globalVariables head modifyList
#' @importFrom graphics hist lines abline plot
#' @importFrom stats density na.omit rnorm as.formula model.matrix plogis qnorm
#' @importFrom stats rbinom reformulate rpois sd terms vcov
NULL

# Suppress R CMD check warnings about global variables
utils::globalVariables(c("means", "value", "x", "y", ".data"))