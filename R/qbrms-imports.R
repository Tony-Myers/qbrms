# =============================================================================
# R/qbrms-imports.R - centralised import directives for NAMESPACE
# =============================================================================

#' Internal import directives for qbrms
#'
#'
#' @name qbrms-imports
#' @keywords internal
#'
#' @import ggplot2
#' @import mvtnorm
#' @importFrom cowplot plot_grid
#' @importFrom graphics abline curve hist par
#' @importFrom lme4 nobars findbars
#' @importFrom methods is
#' @importFrom patchwork wrap_plots
#' @importFrom posterior as_draws_df
#' @importFrom scales percent_format
#' @importFrom stats as.formula binomial coef complete.cases dnorm gaussian lm median
#' @importFrom stats model.frame model.matrix model.response plogis qnorm quantile
#' @importFrom stats rbinom rnorm rpois sd vcov terms reformulate rcauchy rt runif var
#' @importFrom stats fitted qqline qqnorm reshape residuals glm cor delete.response
#' @importFrom utils head modifyList capture.output getFromNamespace
#' @importFrom TMB MakeADFun compile dynlib
NULL

# Null-coalescing operator used throughout the package
`%||%` <- function(x, y) if (is.null(x)) y else x