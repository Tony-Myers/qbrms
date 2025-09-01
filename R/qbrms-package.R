# =============================================================================
# R/qbrms-package.R
# =============================================================================

#' qbrms: Quick Bayesian Regression Models using INLA
#'
#' @description 
#' The qbrms package provides a brms-like interface for fitting Bayesian 
#' regression models using INLA (Integrated Nested Laplace Approximations). 
#' It offers faster model fitting while maintaining familiar brms syntax and 
#' output formats.
#'
#' @details
#' The main function is \code{\link{qbrms}} which fits Bayesian models using
#' INLA with brms-like syntax. The package supports:
#' \itemize{
#'   \item Fixed and mixed effects models
#'   \item Multiple probability distributions
#'   \item Conditional effects plots
#'   \item Posterior predictive checks
#'   \item Summary methods compatible with brms
#' }
#'
#' @author Tony Myers
#' @docType package
#' @name qbrms-package
#' @keywords package
"_PACKAGE"