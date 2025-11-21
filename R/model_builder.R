# =============================================================================
# R/model_builder.R
# =============================================================================

#' Interactive Model Builder for qbrms (console)
#'
#' @description
#' An interactive assistant that guides users through model specification by
#' asking questions about their data, suggesting appropriate families, helping
#' with prior selection, and building qbrms model code.
#'
#' @param data A data frame containing the variables to be modelled (optional).
#'   If not provided, the user will be prompted to specify it.
#' @param response Character string specifying the response variable name (optional).
#' @param predictors Character vector of predictor variable names (optional).
#' @param quiet Logical; if TRUE, suppresses welcome messages (default: FALSE).
#'
#' @return An object with class \code{"qbrms_model_spec"} containing:
#' \itemize{
#'   \item \code{formula}: The model formula
#'   \item \code{family}: A list with \code{name} and the constructed family \code{object}
#'   \item \code{prior}: Prior specifications (if provided)
#'   \item \code{data}: The data frame
#'   \item \code{data_name}: The symbol used for data in the emitted code
#'   \item \code{model_code}: Character string with executable qbrms code
#'   \item \code{response_info}: Summary information about the response
#' }
#'
#' @examples
#' \dontrun{
#' spec <- model_builder()
#' fit <- eval(parse(text = spec$model_code))
#' }
#'
#' @export
model_builder <- function(data = NULL, response = NULL, predictors = NULL, quiet = FALSE) {
  
  if (!quiet) {
    cat("\n")
    cat("========================================\n")
    cat("   qbrms Interactive Model Builder\n")
    cat("========================================\n\n")
    cat("I shall guide you through building your Bayesian model.\n")
    cat("Press Ctrl+C at any time to exit.\n\n")
  }
  
  # Step 1: Get or validate data
  data_info <- if (is.null(data)) {
    .mb_get_data()
  } else {
    data_name <- deparse(substitute(data))
    list(data = data, name = data_name)
  }
  
  data      <- data_info$data
  data_name <- data_info$name
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame", call. = FALSE)
  }
  
  # Step 2: Identify response variable
  if (is.null(response)) {
    response <- .mb_get_response(data)
  } else {
    if (!response %in% names(data)) {
      stop("Response variable '", response, "' not found in data", call. = FALSE)
    }
  }
  
  # Step 3: Characterise response variable and suggest family
  response_info     <- .mb_characterise_response(data[[response]], response)
  suggested_families <- .mb_suggest_families(response_info)
  
  # Step 4: Allow user to select or confirm family
  family_choice <- .mb_select_family(suggested_families, response_info)
  
  # Step 5: Get predictor variables
  if (is.null(predictors)) {
    predictors <- .mb_get_predictors(data, response)
  }
  
  # Step 6: Check for random effects
  random_effects <- .mb_get_random_effects(data, predictors)
  
  # Step 7: Build formula
  formula_obj <- .mb_build_formula(response, predictors, random_effects)
  
  # Step 8: Prior specification
  prior_spec <- .mb_get_priors(family_choice, predictors, random_effects)
  
  # Step 9: Generate model code
  model_code <- .mb_generate_code(formula_obj, family_choice, prior_spec, data_name)
  
  # Step 10: Present summary and options
  .mb_present_summary(formula_obj, family_choice, prior_spec, model_code)
  
  # Output object
  result <- structure(
    list(
      formula       = formula_obj,
      family        = family_choice,
      prior         = prior_spec,
      data          = data,
      data_name     = data_name,
      model_code    = model_code,
      response_info = response_info
    ),
    class = "qbrms_model_spec"
  )
  
  invisible(result)
}

# =============================================================================
# Small helpers
# =============================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Resolve which fitting entry point the package exports
.qbrms_entry_fn <- function() {
  ex <- getNamespaceExports("qbrms")
  if ("qbrm"  %in% ex) return("qbrm")
  if ("qbrms" %in% ex) return("qbrms")
  stop("Neither 'qbrm()' nor 'qbrms()' is exported by the qbrms package.", call. = FALSE)
}

# Translate UI/labels to actual constructor names in qbrms
# Defined here to keep the console builder self-contained.
#' @keywords internal
.qbrms_resolve_family_ctor <- function(name) {
  if (is.null(name) || !nzchar(name)) return("gaussian")
  key <- tolower(gsub("[^a-z0-9]", "", name))
  
  map <- c(
    gaussian = "gaussian",
    normal   = "gaussian",
    
    student  = "student",
    studentt = "student",
    t        = "student",
    
    skewnormal   = "sn",
    skew_normal  = "sn",
    sn           = "sn",
    
    beta      = "Beta",        # capitalised constructor in your package
    lognormal = "lognormal",
    gamma     = "gamma",
    binomial  = "binomial",
    multinomial = "multinomial",
    poisson   = "poisson",
    
    negbinomial      = "nbinomial",
    negativebinomial = "nbinomial",
    nbinomial        = "nbinomial",
    
    zeroinflatedpoisson   = "zero_inflated_poisson",
    zeroinflatedpoisson1  = "zero_inflated_poisson",
    zip                   = "zip",
    zeroinflatednegbinomial  = "zero_inflated_negbinomial",
    zeroinflatednegbinomial1 = "zero_inflated_negbinomial",
    zinb                  = "zinb",
    
    weibull        = "weibull",
    exponential    = "exponential",
    weibullsurv    = "weibull",
    exponentialsurv= "exponential",
    
    simplex        = "simplex",
    gev            = "gev",
    gumbel         = "gumbel",
    circularnormal = "circular_normal",
    circular_normal= "circular_normal",
    vonmises       = "von_mises",
    von_mises      = "von_mises",
    laplace        = "laplace",
    doubleexponential = "double_exponential",
    genstudentt    = "gen_student_t",
    generalizedt   = "gen_student_t",
    
    cumulative     = "cumulative"  # handled by routing if available
  )
  
  ctor <- unname(map[[key]])
  if (is.null(ctor)) ctor <- "gaussian"
  
  ex <- try(getNamespaceExports("qbrms"), silent = TRUE)
  if (!inherits(ex, "try-error") && !(ctor %in% ex)) {
    if ("gaussian" %in% ex) return("gaussian")
  }
  ctor
}

# =============================================================================
# HELPER FUNCTIONS FOR MODEL BUILDER
# =============================================================================

#' Get Data from User
#' @keywords internal
.mb_get_data <- function() {
  cat("Step 1: Data\n")
  cat("------------\n")
  cat("Please specify your data frame name: ")
  
  data_name <- trimws(readline())
  
  if (is.null(data_name) || length(data_name) == 0 || nchar(data_name) == 0) {
    stop("No data frame name provided", call. = FALSE)
  }
  
  if (!is.character(data_name) || length(data_name) != 1) {
    stop("Invalid data frame name", call. = FALSE)
  }
  
  if (!exists(data_name, envir = .GlobalEnv)) {
    stop("Object '", data_name, "' not found in your environment", call. = FALSE)
  }
  
  data <- get(data_name, envir = .GlobalEnv)
  
  if (!is.data.frame(data)) {
    stop("'", data_name, "' is not a data frame", call. = FALSE)
  }
  
  cat("Data loaded:", nrow(data), "observations,", ncol(data), "variables\n\n")
  return(list(data = data, name = data_name))
}

#' Get Response Variable
#' @keywords internal
.mb_get_response <- function(data) {
  cat("Step 2: Response Variable\n")
  cat("-------------------------\n")
  cat("Available variables:\n")
  cat(paste0("  ", seq_along(names(data)), ". ", names(data), "\n"))
  cat("\n")
  cat("Which variable is your response (outcome)? Enter name or number: ")
  
  choice <- trimws(readline())
  
  if (is.null(choice) || length(choice) == 0 || nchar(choice) == 0) {
    stop("No response variable selected", call. = FALSE)
  }
  
  if (grepl("^[0-9]+$", choice)) {
    idx <- as.integer(choice)
    if (is.na(idx) || idx < 1 || idx > ncol(data)) {
      stop("Invalid selection", call. = FALSE)
    }
    response <- names(data)[idx]
  } else {
    if (!choice %in% names(data)) {
      stop("Variable '", choice, "' not found in data", call. = FALSE)
    }
    response <- choice
  }
  
  cat("Response variable:", response, "\n\n")
  return(response)
}

#' Characterise Response Variable
#' @keywords internal
.mb_characterise_response <- function(y, var_name) {
  
  cat("Step 3: Understanding Your Response Variable\n")
  cat("--------------------------------------------\n")
  
  info <- list(
    name       = var_name,
    class      = class(y)[1],
    n_obs      = length(y),
    n_missing  = sum(is.na(y))
  )
  
  y_complete <- y[!is.na(y)]
  
  if (is.numeric(y)) {
    info$type            <- "numeric"
    info$range           <- range(y_complete)
    info$mean            <- mean(y_complete)
    info$sd              <- sd(y_complete)
    info$n_unique        <- length(unique(y_complete))
    info$has_zeros       <- any(y_complete == 0)
    info$all_positive    <- all(y_complete > 0)
    info$all_nonnegative <- all(y_complete >= 0)
    info$all_integers    <- all(y_complete == floor(y_complete))
    info$bounded_01      <- all(y_complete > 0 & y_complete < 1)
    info$excess_zeros    <- sum(y_complete == 0) / length(y_complete) > 0.1
    
    cat("Type: Numeric\n")
    cat("Range: [", info$range[1], ", ", info$range[2], "]\n", sep = "")
    cat("Mean (SD): ", round(info$mean, 2), " (", round(info$sd, 2), ")\n", sep = "")
    cat("Unique values:", info$n_unique, "\n")
    
  } else if (is.factor(y) || is.character(y)) {
    if (is.character(y)) y <- factor(y)
    info$type       <- "categorical"
    info$levels     <- levels(y)
    info$n_levels   <- length(info$levels)
    info$table      <- table(y)
    info$is_ordered <- is.ordered(y)
    info$is_binary  <- info$n_levels == 2
    
    cat("Type: Categorical\n")
    cat("Levels:", info$n_levels, "\n")
    cat("Categories:", paste(info$levels, collapse = ", "), "\n")
    if (info$is_ordered) cat("Ordered: Yes\n")
    
  } else {
    stop("Response variable must be numeric or categorical", call. = FALSE)
  }
  
  cat("\n")
  return(info)
}

#' Suggest Appropriate Families
#' @keywords internal
.mb_suggest_families <- function(response_info) {
  suggestions <- list()
  
  if (response_info$type == "numeric") {
    
    if (!response_info$bounded_01 && response_info$n_unique > 20) {
      suggestions[[length(suggestions) + 1]] <- list(
        family = "gaussian", reason = "Continuous data with no constraints", priority = 1
      )
      suggestions[[length(suggestions) + 1]] <- list(
        family = "student", reason = "Robust to outliers (heavier tails than normal)", priority = 2
      )
      suggestions[[length(suggestions) + 1]] <- list(
        family = "skew_normal", reason = "If your data appears skewed", priority = 3
      )
    }
    
    if (response_info$all_positive && !response_info$all_integers) {
      suggestions[[length(suggestions) + 1]] <- list(
        family = "lognormal", reason = "Positive continuous data", priority = 1
      )
      suggestions[[length(suggestions) + 1]] <- list(
        family = "gamma", reason = "Positive continuous with flexible shape", priority = 2
      )
    }
    
    if (response_info$bounded_01) {
      suggestions[[length(suggestions) + 1]] <- list(
        family = "beta", reason = "Proportions or rates in (0,1)", priority = 1
      )
    }
    
    if (response_info$all_integers && response_info$all_nonnegative) {
      if (!response_info$excess_zeros) {
        suggestions[[length(suggestions) + 1]] <- list(
          family = "poisson", reason = "Count data (non-negative integers)", priority = 1
        )
        suggestions[[length(suggestions) + 1]] <- list(
          family = "negbinomial", reason = "Overdispersed count data", priority = 2
        )
      } else {
        suggestions[[length(suggestions) + 1]] <- list(
          family = "zero_inflated_poisson", reason = "Count data with excess zeros", priority = 1
        )
        suggestions[[length(suggestions) + 1]] <- list(
          family = "zero_inflated_negbinomial", reason = "Overdispersed with excess zeros", priority = 2
        )
      }
    }
    
  } else if (response_info$type == "categorical") {
    
    if (isTRUE(response_info$is_binary)) {
      suggestions[[length(suggestions) + 1]] <- list(
        family = "binomial", reason = "Binary outcome", priority = 1
      )
    } else if (isTRUE(response_info$is_ordered)) {
      suggestions[[length(suggestions) + 1]] <- list(
        family = "cumulative", reason = "Ordered categorical outcome", priority = 1
      )
    } else {
      suggestions[[length(suggestions) + 1]] <- list(
        family = "multinomial", reason = "Unordered categorical outcome", priority = 1
      )
    }
  }
  
  suggestions <- suggestions[order(sapply(suggestions, function(x) x$priority))]
  return(suggestions)
}

#' Let User Select Family
#' @keywords internal
.mb_select_family <- function(suggestions, response_info) {
  cat("Step 4: Selecting Distribution Family\n")
  cat("--------------------------------------\n")
  cat("Based on your data characteristics, I recommend:\n\n")
  
  for (i in seq_along(suggestions)) {
    cat(i, ". ", suggestions[[i]]$family, "\n", sep = "")
    cat("   Reason: ", suggestions[[i]]$reason, "\n", sep = "")
  }
  
  cat("\nEnter number to select, or type a family name: ")
  choice <- trimws(readline())
  
  if (is.null(choice) || length(choice) == 0 || nchar(choice) == 0) {
    warning("No selection made, using first suggestion")
    choice <- "1"
  }
  
  if (grepl("^[0-9]+$", choice)) {
    idx <- as.integer(choice)
    if (is.na(idx) || idx < 1 || idx > length(suggestions)) {
      warning("Invalid selection, using first suggestion")
      idx <- 1
    }
    family_name <- suggestions[[idx]]$family
  } else {
    family_name <- tolower(choice)
  }
  
  family_obj <- .mb_create_family_object(family_name)
  
  cat("Selected family:", family_name, "\n\n")
  return(list(name = family_name, object = family_obj))
}

#' Create Family Object
#' @keywords internal
.mb_create_family_object <- function(family_name) {
  
  family_map <- list(
    gaussian    = gaussian,
    normal      = gaussian,
    
    student     = student,     # your package provides student()
    t           = student,
    
    skew_normal = sn,          # skew-normal constructor exposed as sn()
    sn          = sn,
    
    Beta        = Beta,        # capitalised
    beta        = Beta,
    
    lognormal   = lognormal,
    gamma       = gamma,
    binomial    = binomial,
    multinomial = multinomial,
    poisson     = poisson,
    
    negbinomial = nbinomial,
    nbinomial   = nbinomial,
    
    zero_inflated_poisson     = zero_inflated_poisson,
    zip                       = zip,
    zero_inflated_negbinomial = zero_inflated_negbinomial,
    zinb                      = zinb,
    
    weibull     = weibull,
    exponential = exponential,
    simplex     = simplex,
    
    cumulative  = function() list(family = "cumulative") # placeholder for routing
  )
  
  key <- tolower(gsub("[^a-z0-9_]", "", family_name))
  if (key %in% names(family_map)) {
    return(family_map[[key]]())
  } else {
    warning("Family '", family_name, "' not recognised, using gaussian()")
    return(gaussian())
  }
}

#' Get Predictor Variables
#' @keywords internal
.mb_get_predictors <- function(data, response) {
  cat("Step 5: Predictor Variables\n")
  cat("---------------------------\n")
  
  available_vars <- setdiff(names(data), response)
  
  cat("Available predictors:\n")
  for (i in seq_along(available_vars)) {
    var_class <- class(data[[available_vars[i]]])[1]
    cat(sprintf("  %2d. %-20s (%s)\n", i, available_vars[i], var_class))
  }
  
  cat("\nEnter predictor numbers separated by spaces (e.g., '1 3 5'),\n")
  cat("or press Enter to choose manually later: ")
  
  choice <- trimws(readline())
  
  if (is.null(choice) || length(choice) == 0 || nchar(choice) == 0) {
    predictors <- character(0)
  } else {
    indices <- suppressWarnings(as.integer(strsplit(choice, "\\s+")[[1]]))
    indices <- indices[!is.na(indices)]
    indices <- indices[indices >= 1 & indices <= length(available_vars)]
    predictors <- if (length(indices) == 0) character(0) else available_vars[indices]
  }
  
  if (length(predictors) == 0) {
    cat("No predictors selected now. You may add them later.\n\n")
    predictors <- "1"  # Intercept only
  } else {
    cat("Selected predictors:", paste(predictors, collapse = ", "), "\n\n")
  }
  
  return(predictors)
}

#' Get Random Effects Structure
#' @keywords internal
.mb_get_random_effects <- function(data, predictors) {
  cat("Step 6: Random Effects (Mixed Models)\n")
  cat("-------------------------------------\n")
  cat("Do you want to include random effects? (y/n): ")
  
  include_random <- tolower(trimws(readline()))
  
  if (is.null(include_random) || length(include_random) == 0 ||
      nchar(include_random) == 0 || !include_random %in% c("y", "yes")) {
    cat("No random effects included.\n\n")
    return(NULL)
  }
  
  potential_groups <- character(0)
  vars_to_scan <- if (length(predictors) == 1 && predictors[1] == "1") names(data) else predictors
  for (var in vars_to_scan) {
    if (is.factor(data[[var]]) || is.character(data[[var]])) {
      n_levels <- length(unique(data[[var]]))
      if (n_levels >= 3 && n_levels <= nrow(data) / 2) {
        potential_groups <- c(potential_groups, var)
      }
    }
  }
  
  if (length(potential_groups) == 0) {
    cat("No suitable grouping variables found.\n\n")
    return(NULL)
  }
  
  cat("Potential grouping variables:\n")
  for (i in seq_along(potential_groups)) {
    n_levels <- length(unique(data[[potential_groups[i]]]))
    cat(sprintf("  %d. %s (%d groups)\n", i, potential_groups[i], n_levels))
  }
  
  cat("\nEnter grouping variable number: ")
  choice <- trimws(readline())
  choice_int <- suppressWarnings(as.integer(choice))
  
  if (is.null(choice) || length(choice) == 0 || nchar(choice) == 0 ||
      is.na(choice_int) || choice_int < 1 || choice_int > length(potential_groups)) {
    cat("Invalid choice, no random effects included.\n\n")
    return(NULL)
  }
  
  grouping_var <- potential_groups[choice_int]
  
  cat("\nRandom intercept only, or random slopes? (intercept/slopes): ")
  re_type <- tolower(trimws(readline()))
  
  if (!is.null(re_type) && length(re_type) > 0 && nchar(re_type) > 0 && re_type == "slopes") {
    cat("Which predictor should have a random slope? Enter name: ")
    slope_var <- trimws(readline())
    
    if (is.null(slope_var) || length(slope_var) == 0 || nchar(slope_var) == 0 ||
        !slope_var %in% names(data)) {
      warning("Variable not found, using random intercept only")
      random_formula <- paste0("(1 | ", grouping_var, ")")
    } else {
      random_formula <- paste0("(", slope_var, " | ", grouping_var, ")")
    }
  } else {
    random_formula <- paste0("(1 | ", grouping_var, ")")
  }
  
  cat("Random effects:", random_formula, "\n\n")
  return(random_formula)
}

#' Build Formula
#' @keywords internal
.mb_build_formula <- function(response, predictors, random_effects = NULL) {
  
  if (length(predictors) == 1 && predictors[1] == "1") {
    fixed_part <- "1"
  } else {
    fixed_part <- paste(predictors, collapse = " + ")
  }
  
  if (!is.null(random_effects)) {
    formula_str <- paste(response, "~", fixed_part, "+", random_effects)
  } else {
    formula_str <- paste(response, "~", fixed_part)
  }
  
  as.formula(formula_str)
}

#' Get Prior Specifications
#' @keywords internal
.mb_get_priors <- function(family_choice, predictors, random_effects) {
  cat("Step 7: Prior Specification\n")
  cat("---------------------------\n")
  cat("Would you like to specify custom priors? (y/n): ")
  cat("[Default priors will be used if you select 'n']\n")
  
  specify_priors <- tolower(trimws(readline()))
  
  if (is.null(specify_priors) || length(specify_priors) == 0 ||
      nchar(specify_priors) == 0 || !specify_priors %in% c("y", "yes")) {
    cat("Using default priors.\n\n")
    return(NULL)
  }
  
  cat("\nPrior options:\n")
  cat("1. Weakly informative (recommended for most cases)\n")
  cat("2. Informative (if you have strong prior knowledge)\n")
  cat("3. Flat/vague (minimal prior information)\n")
  cat("Note: Specific prior values must be added manually to the generated code\n")
  cat("\nSelect option (1-3): ")
  
  prior_choice     <- trimws(readline())
  prior_choice_int <- suppressWarnings(as.integer(prior_choice))
  
  if (is.null(prior_choice) || length(prior_choice) == 0 ||
      nchar(prior_choice) == 0 || is.na(prior_choice_int) || !prior_choice_int %in% 1:3) {
    cat("Invalid choice, using default priors.\n\n")
    return(NULL)
  }
  
  prior_type <- c("weakly_informative", "informative", "flat")[prior_choice_int]
  
  cat("Prior type selected:", prior_type, "\n")
  cat("You will need to add specific prior values to the generated code.\n\n")
  
  list(type = prior_type, specifications = NULL)
}

#' Generate Model Code
#' @keywords internal
.mb_generate_code <- function(formula, family, prior, data_name) {
  
  entry <- .qbrms_entry_fn()
  
  fam_label   <- (family$name %||% family)
  ctor_sym    <- .qbrms_resolve_family_ctor(fam_label)
  family_call <- paste0("qbrms::", ctor_sym, "()")
  
  args <- c(
    paste0("  formula = ", deparse(formula)),
    paste0("  data = ", data_name),
    paste0("  family = ", family_call)
  )
  
  if (!is.null(prior) && !is.null(prior$specifications)) {
    args <- c(args, paste0("  prior = ", prior$specifications))
  }
  
  code <- paste0(
    "fit <- qbrms::", entry, "(\n",
    paste(args, collapse = ",\n"),
    "\n)"
  )
  
  code
}

#' Present Summary
#' @keywords internal
.mb_present_summary <- function(formula, family, prior, model_code) {
  
  cat("\n")
  cat("========================================\n")
  cat("   Model Specification Complete\n")
  cat("========================================\n\n")
  
  cat("Formula:", deparse(formula), "\n")
  cat("Family:", family$name, "\n")
  cat("Priors:", if (is.null(prior)) "Default" else prior$type, "\n\n")
  
  cat("Generated Code:\n")
  cat("---------------\n")
  cat(model_code, "\n\n")
  
  cat("Copy the code above to fit your model.\n")
  cat("To run it now, assign the result to a variable:\n")
  cat("  spec <- model_builder(...)\n")
  cat("  fit <- eval(parse(text = spec$model_code))\n\n")
}

#' Print Method for qbrms_model_spec
#' @param x A qbrms_model_spec object
#' @param ... Additional arguments (unused)
#' @export
print.qbrms_model_spec <- function(x, ...) {
  cat("qbrms Model Specification\n")
  cat("=========================\n\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Family:", x$family$name, "\n")
  cat("Data:",   x$data_name, "\n\n")
  cat("Code to execute:\n")
  cat(x$model_code, "\n")
  invisible(x)
}
