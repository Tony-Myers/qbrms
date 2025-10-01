# R/imports.R - all imports

#' @importFrom stats as.formula binomial coef complete.cases dnorm gaussian lm median
#' @importFrom stats model.frame model.matrix model.response plogis qnorm quantile
#' @importFrom stats rbinom rnorm rpois sd vcov terms reformulate rcauchy rt runif var
#' @importFrom stats fitted qqline qqnorm reshape residuals glm cor
#' @importFrom utils head modifyList capture.output getFromNamespace
#' @importFrom graphics abline curve hist par
#' @importFrom methods is
#' @import ggplot2
#' @importFrom scales percent_format
#' @importFrom posterior as_draws_df
#' @import mvtnorm
#' @importFrom stats delete.response
#' @importFrom cowplot plot_grid
#' @importFrom patchwork wrap_plots

# If you prefer to add via usethis commands, run these:
# usethis::use_package("cowplot", "Suggests")
# usethis::use_package("patchwork", "Suggests")

# Add to your main package documentation file (typically R/qbrms-package.R):
#' @importFrom stats delete.response model.matrix terms as.formula
NULL
NULL

# Handle ggplot2 / NSE variable bindings to avoid R CMD check notes
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "draw", "y", "value", "rep_id", "x", "group", "median", "quantile",
    "obs", "lower", "upper", ".data", "d",
    "Model", "Type", "Value", "effect1__", "estimate__", "id__", "lower__", "upper__",
    "type", "ps", "Parameter", "sig_level",
    "CI_low", "CI_high", "Dist",
    "Group", "Success_Rate", "Group_Size", "Lower", "Upper", "groups",
    "is_pseudo", "group_id", "regularisation_strength", "prior_random_intercept"
  ))
}

# TMB imports for ordinal regression
if (requireNamespace("TMB", quietly = TRUE)) {
  #' @importFrom TMB MakeADFun compile dynlib
  NULL
}

# lme4 imports for formula parsing in ordinal models
if (requireNamespace("lme4", quietly = TRUE)) {
  #' @importFrom lme4 nobars findbars
  NULL
}


# Null-coalescing operator used throughout the package
`%||%` <- function(x, y) if (is.null(x)) y else x
