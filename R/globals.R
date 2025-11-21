# =============================================================================
# R/globals.R - globalVariables declarations
# =============================================================================

#' Internal globals for qbrms
#'
#' Declares symbols that are used non-standardly (e.g. in NSE, Shiny) so that
#' R CMD check does not flag them as undefined global variables.
#'
#' @name qbrms-globals
#' @keywords internal
#' @importFrom utils globalVariables
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      # Shiny NSE and addin-related objects
      "shiny"
      # Add further symbols here if R CMD check complains about them
    )
  )
}

# Declare all global variables in a single call
utils::globalVariables(c(
  # General package variables
  "binary_models",
  "thresholds",
  "Estimate",
  "Q_low",
  "Q_high",
  ".data",
  "format_number",
  "format_percentage",
  "format_bf",
  "extract_model_info",
  "Group",
  "Success_Rate",
  "Group_Size",
  "Lower",
  "Upper",
  "groups",
  "is_pseudo",
  "group_id",
  "regularisation_strength",
  "prior_random_intercept",
  "Parameter",
  "mean",
  "sd",
  "severity",
  "message",
  "recommendation",
  "overall_status",
  "fitting_algorithm",
  "data_diagnostics",
  "fitting_diagnostics",
  "control.compute",
  
  # ggplot2 aesthetics and NSE variables
  "density",
  "x",
  "y",
  "ymin",
  "ymax",
  "value",
  "draw",
  "rep_id",
  "group",
  "median",
  "quantile",
  "obs",
  "lower",
  "upper",
  "d",
  
  # Model comparison variables
  "model",
  "ic",
  "se",
  "weight",
  "model1",
  "model2",
  "correlation",
  "Model",
  "Type",
  "Value",
  
  # Effect and estimate variables
  "effect1__",
  "estimate__",
  "id__",
  "lower__",
  "upper__",
  "type",
  "ps",
  "sig_level",
  "CI_low",
  "CI_high",
  "Dist",
  
  # Diagnostics variables
  "index",
  "std_residual",
  "standardised",
  "influential",
  "residuals",
  "observed",
  "fitted"
))