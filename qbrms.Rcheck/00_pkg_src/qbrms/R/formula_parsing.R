# =============================================================================
# R/formula_parsing.R
# =============================================================================

#' Parse Formula Components
#'
#' @description
#' Parse formula to identify random effects, binomial trials, etc.
#'
#' @param formula Model formula
#' @param data Data frame
#' @return List with formula components information
#'
#' @keywords internal
parse_formula_components <- function(formula, data) {
  formula_str <- deparse(formula, width.cutoff = 500L)
  has_pipe <- grepl("\\|", formula_str)
  
  if (has_pipe) {
    has_random_effects <- grepl("\\([^|]*\\|[^)]*\\)", formula_str)
    is_binomial_trials <- grepl("\\w+\\s*\\|\\s*trials\\(", formula_str)
    
    # Extract trials information if present
    trials_info <- NULL
    if (is_binomial_trials) {
      trials_match <- regmatches(formula_str, regexpr("trials\\(([^)]+)\\)", formula_str))
      if (length(trials_match) > 0) {
        trials_var <- gsub("trials\\(|\\)", "", trials_match)
        trials_info <- list(
          variable = trials_var,
          values = if (trials_var %in% names(data)) data[[trials_var]] else NULL
        )
      }
    }
  } else {
    has_random_effects <- FALSE
    is_binomial_trials <- FALSE
    trials_info <- NULL
  }
  
  return(list(
    has_random_effects = has_random_effects,
    is_binomial_trials = is_binomial_trials,
    formula_str = formula_str,
    trials_info = trials_info
  ))
}

#' Parse brms Formula Objects
#'
#' @description
#' Parse brms formula objects including bf() specifications
#'
#' @param formula Formula or brms formula object
#' @return List with parsed formula information
#'
#' @keywords internal
parse_brms_formula <- function(formula) {
  if (inherits(formula, "brmsterms")) {
    result <- list(
      response_formula = formula$formula,
      distributional_formulas = list(),
      is_distributional = length(formula$dpars) > 0 || length(formula$nlpars) > 0
    )
    
    if (length(formula$dpars) > 0) {
      for (dpar in names(formula$dpars)) {
        result$distributional_formulas[[dpar]] <- formula$dpars[[dpar]]$formula
      }
    }
    return(result)
  } else if (inherits(formula, "brmsformula")) {
    return(list(
      response_formula = formula$formula,
      distributional_formulas = list(),
      is_distributional = FALSE
    ))
  } else {
    return(list(
      response_formula = formula,
      distributional_formulas = list(),
      is_distributional = FALSE
    ))
  }
}

#' Get Predictor Variables from Formula
#'
#' @description
#' Extract predictor variables from formula, categorised by type
#'
#' @param formula Model formula
#' @param data Data frame
#' @return List with predictor variable information
#'
#' @keywords internal
get_predictor_variables <- function(formula, data) {
  # Handle different formula types safely
  if (is.character(formula)) {
    formula <- as.formula(formula)
  }
  
  # Extract variables from formula
  tryCatch({
    rhs <- formula[[3]]
    formula_str <- deparse(rhs, width.cutoff = 500L)
    # Remove random effects and f() functions
    fixed_part <- gsub("\\s*\\+\\s*\\([^)]+\\)", "", formula_str)
    fixed_part <- gsub("\\s*\\+\\s*f\\([^)]+\\)", "", fixed_part)
    
    vars <- all.vars(as.formula(paste("~", fixed_part)))
    vars <- vars[vars != "1" & vars != "group_id"]
  }, error = function(e) {
    vars <- all.vars(formula)[-1]
    vars <- vars[vars != "group_id"]
  })
  
  numeric_vars <- character(0)
  categorical_vars <- character(0)
  
  for (var in vars) {
    if (var %in% names(data)) {
      if (is.numeric(data[[var]])) {
        numeric_vars <- c(numeric_vars, var)
      } else {
        categorical_vars <- c(categorical_vars, var)
      }
    }
  }
  
  return(list(
    all_vars = vars,
    numeric_vars = numeric_vars,
    categorical_vars = categorical_vars
  ))
}
