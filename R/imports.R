# R/imports.R 

#' @importFrom grDevices rgb
#' @importFrom graphics hist lines abline plot
#' @importFrom stats density dnorm rnorm as.formula model.matrix plogis qnorm rbinom reformulate rpois sd terms vcov coef
#' @importFrom utils globalVariables head modifyList
NULL

# Suppress R CMD check warnings about global variables
utils::globalVariables(c("means", "value", "x", "y", ".data"))