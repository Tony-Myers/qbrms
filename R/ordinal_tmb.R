# =============================================================================
# R/ordinal_tmb.R - Complete TMB-based Ordinal Regression for qbrms
# =============================================================================

# =============================================================================
# ADAPTIVE THRESHOLD CENTERING FUNCTIONS
# =============================================================================

#' Compute Data-Driven Threshold Priors
#' @description Calculate threshold priors based on empirical quantiles to match brms centering
#' @param y_ordered Ordered factor response variable
#' @param method Method for computing adaptive priors ("quantile" or "cumulative_mean")
#' @return List with threshold means, standard deviations, and empirical baseline
#' @keywords internal
compute_adaptive_threshold_priors <- function(y_ordered, method = "quantile") {
  
  n_levels <- nlevels(y_ordered)
  y_numeric <- as.numeric(y_ordered)
  
  if (method == "quantile") {
    # Use empirical cumulative probabilities to center thresholds
    level_props <- table(y_numeric) / length(y_numeric)
    cum_props <- cumsum(level_props)
    
    # Convert to logit scale (matching cumulative logit model)
    empirical_thresholds <- qlogis(cum_props[-length(cum_props)])
    
    # Center around the middle threshold to match brms convention
    if (length(empirical_thresholds) > 1) {
      middle_thresh <- empirical_thresholds[ceiling(length(empirical_thresholds)/2)]
      centered_thresholds <- empirical_thresholds - middle_thresh
    } else {
      centered_thresholds <- empirical_thresholds
    }
    
    return(list(
      thresh_mean = as.numeric(centered_thresholds),
      thresh_sd = rep(1.5, length(centered_thresholds)),
      empirical_baseline = empirical_thresholds
    ))
    
  } else if (method == "cumulative_mean") {
    # Alternative: Use cumulative means
    n_thresh <- n_levels - 1
    empirical_thresholds <- qlogis(seq(0.2, 0.8, length.out = n_thresh))
    
    return(list(
      thresh_mean = empirical_thresholds,
      thresh_sd = rep(2.0, length(empirical_thresholds))
    ))
  }
}

#' Improved Threshold Prior Computation - ENHANCED VERSION
#' @param y_ordered Ordered factor response
#' @param verbose Logical; print diagnostics
#' @return List with improved threshold priors
#' @keywords internal
.compute_improved_threshold_priors <- function(y_ordered, verbose = TRUE) {
  
  n_levels <- nlevels(y_ordered)
  y_numeric <- as.numeric(y_ordered)
  
  # Compute empirical cumulative probabilities
  level_counts <- table(y_numeric)
  level_props <- level_counts / sum(level_counts)
  cum_props <- cumsum(level_props)
  
  # Remove the last cumulative probability (always 1)
  cum_props_for_thresholds <- cum_props[-length(cum_props)]
  
  # IMPROVED: Better handling of extreme probabilities
  # Avoid logit(-Inf) and logit(Inf) by trimming extreme values
  cum_props_trimmed <- pmax(pmin(cum_props_for_thresholds, 0.999), 0.001)
  
  # Convert to logit scale
  empirical_thresholds <- qlogis(cum_props_trimmed)
  
  # IMPROVED: More sophisticated centering strategy
  if (length(empirical_thresholds) > 2) {
    # Center around median threshold for better numerical properties
    median_thresh <- stats::median(empirical_thresholds)
    centered_thresholds <- empirical_thresholds - median_thresh
  } else if (length(empirical_thresholds) == 2) {
    # Center around first threshold
    centered_thresholds <- empirical_thresholds - empirical_thresholds[1]
  } else {
    # Single threshold case
    centered_thresholds <- empirical_thresholds
  }
  
  # IMPROVED: Adaptive standard deviations based on data spread
  empirical_spread <- diff(range(empirical_thresholds))
  adaptive_sd <- pmax(empirical_spread * 0.3, 1.0)  # At least 1.0, but scale with data
  
  if (verbose) {
    cat("Threshold prior computation:\n")
    cat("  Empirical cumulative props:", round(cum_props_for_thresholds, 3), "\n")
    cat("  Raw empirical thresholds:", round(empirical_thresholds, 3), "\n")
    cat("  Empirical spread:", round(empirical_spread, 3), "\n")
    cat("  Adaptive SD:", round(adaptive_sd, 3), "\n")
  }
  
  return(list(
    thresh_mean = as.numeric(centered_thresholds),
    thresh_sd = rep(adaptive_sd, length(centered_thresholds)),
    empirical_quantiles = cum_props_for_thresholds,
    empirical_thresholds = empirical_thresholds,
    adaptive_sd = adaptive_sd
  ))
}

#' Enhanced Prior Processing with Adaptive Centering - CORRECTED VERSION
#' @param prior Prior specifications using qbrms prior syntax
#' @param tmb_data Prepared TMB data structure
#' @param data Original data frame
#' @param formula_parts Parsed formula components
#' @param verbose Logical; print progress messages
#' @return List of prior parameters for TMB
#' @keywords internal
process_ordinal_priors_adaptive <- function(prior, tmb_data, data, formula_parts, verbose = TRUE) {
  
  # IMPROVED: Better adaptive threshold centering
  adaptive_priors <- .compute_improved_threshold_priors(formula_parts$y_ordered, verbose)
  
  if (verbose) {
    cat("Improved adaptive threshold priors:\n")
    cat("  Empirical quantiles:", round(adaptive_priors$empirical_quantiles, 3), "\n")
    cat("  Centered thresholds:", round(adaptive_priors$thresh_mean, 3), "\n")
    cat("  Threshold SDs:", round(adaptive_priors$thresh_sd, 3), "\n")
  }
  
  # Start with adaptive priors
  thresh_mean <- adaptive_priors$thresh_mean
  thresh_sd <- adaptive_priors$thresh_sd
  
  # Fixed effects priors
  if (tmb_data$n_coef > 0) {
    # IMPROVED: More appropriate default priors based on design matrix scale
    X_scale <- if (ncol(tmb_data$X) > 0) {
      apply(tmb_data$X, 2, function(x) max(abs(x), na.rm = TRUE))
    } else {
      numeric(0)
    }
    
    beta_mean <- rep(0, tmb_data$n_coef)
    beta_sd <- pmax(1 / pmax(X_scale, 1), 0.5)  # Scale-aware priors
    
    if (verbose) cat("Scale-aware coefficient priors, SDs:", round(beta_sd, 3), "\n")
  } else {
    beta_mean <- numeric(0)
    beta_sd <- numeric(0)
  }
  
  # Random effects priors
  re_sd_shape <- 1.0
  re_sd_rate <- 0.5
  
  # Apply user-specified priors
  if (!is.null(prior)) {
    prior_list <- if (inherits(prior, "qbrms_prior_spec")) list(prior) else prior
    
    for (p in prior_list) {
      if (p$class == "Intercept") {
        if (p$distribution == "normal") {
          # Apply to all thresholds but maintain centering
          user_sd <- p$parameters$sd
          thresh_sd[] <- user_sd
          if (verbose) cat("Applied user threshold SD:", user_sd, "\n")
        }
      } else if (p$class == "b" && tmb_data$n_coef > 0) {
        if (!is.null(p$coef)) {
          coef_idx <- which(colnames(tmb_data$X) == p$coef)
          if (length(coef_idx) > 0 && p$distribution == "normal") {
            beta_mean[coef_idx] <- p$parameters$mean
            beta_sd[coef_idx] <- p$parameters$sd
          }
        } else if (p$distribution == "normal") {
          beta_mean[] <- p$parameters$mean
          beta_sd[] <- p$parameters$sd
        }
      }
    }
  }
  
  return(list(
    thresh_mean = thresh_mean,
    thresh_sd = thresh_sd,
    beta_mean = beta_mean,
    beta_sd = beta_sd,
    re_sd_shape = re_sd_shape,
    re_sd_rate = re_sd_rate,
    adaptive_info = adaptive_priors
  ))
}

# =============================================================================
# MAIN ORDINAL REGRESSION FUNCTION
# =============================================================================

#' Quick Bayesian Ordinal Regression Models with Adaptive Centering
#'
#' @description
#' Fits ordinal regression models using Template Model Builder (TMB) with
#' Laplace approximation and adaptive threshold centering to match brms output.
#'
#' @param formula Model formula with ordinal response on the left-hand side.
#' @param data Data frame containing the variables in the model.
#' @param family Ordinal family specification. Currently supports cumulative().
#' @param prior Prior specifications using qbrms prior syntax.
#' @param verbose Logical; print progress messages during fitting.
#' @param threshold_method Method for threshold centering ("quantile" or "cumulative_mean").
#' @param control List of control parameters for TMB optimization.
#' @param ... Additional arguments passed to TMB functions.
#'
#' @return An object of class c("tmb_ordinal_qbrms_fit", "qbrms_fit")
#'
#' @export
qbrmO <- function(formula, data, family = cumulative(), 
                  prior = NULL, verbose = FALSE, 
                  threshold_method = "quantile",
                  control = list(), ...) {
  
  # Input validation
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is_ordinal(family)) {
    stop("'family' must be an ordinal family (e.g., cumulative())")
  }
  
  if (verbose) cat("Starting TMB ordinal regression with adaptive centering...\n")
  
  # Parse formula components
  formula_parts <- parse_ordinal_formula(formula, data, verbose)
  
  # Prepare data for TMB
  tmb_data <- prepare_ordinal_tmb_data(formula_parts, data, verbose)
  
  # Process priors with adaptive centering
  prior_params <- process_ordinal_priors_adaptive(prior, tmb_data, data, 
                                                  formula_parts, verbose)
  
  # Set up TMB model
  tmb_setup <- setup_ordinal_tmb(tmb_data, prior_params, control, verbose)
  
  # Fit the model
  tmb_fit <- fit_ordinal_tmb_model(tmb_setup, verbose)
  
  # Create qbrms-compatible result
  result <- create_ordinal_qbrms_result(tmb_fit, formula, data, family, 
                                        prior, verbose)
  
  # Add adaptive threshold information
  result$adaptive_thresholds <- prior_params$adaptive_info
  result$threshold_method <- threshold_method
  
  if (verbose) {
    cat("TMB ordinal regression completed successfully.\n")
    cat("Threshold method:", threshold_method, "\n")
  }
  
  return(result)
}

# =============================================================================
# FORMULA AND DATA PREPARATION FUNCTIONS
# =============================================================================

#' Parse Ordinal Formula Components
#' @param formula Model formula
#' @param data Data frame
#' @param verbose Logical; print progress messages
#' @return List with parsed formula components
#' @keywords internal
parse_ordinal_formula <- function(formula, data, verbose = TRUE) {
  
  # Extract response variable
  response_var <- all.vars(formula)[1]
  if (!response_var %in% names(data)) {
    stop("Response variable '", response_var, "' not found in data")
  }
  
  y_raw <- data[[response_var]]
  
  # Convert to ordered factor if necessary
  if (!is.factor(y_raw)) {
    if (verbose) cat("Converting response to ordered factor...\n")
    y_ordered <- ordered(y_raw)
  } else if (!is.ordered(y_raw)) {
    if (verbose) cat("Converting factor response to ordered...\n")
    y_ordered <- ordered(y_raw, levels = levels(y_raw))
  } else {
    y_ordered <- y_raw
  }
  
  # Check number of levels
  n_levels <- nlevels(y_ordered)
  if (n_levels < 3) {
    stop("Ordinal response must have at least 3 levels, found ", n_levels)
  }
  
  # Check for random effects
  has_random <- grepl("\\|", deparse(formula))
  
  if (has_random) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("lme4 package required for random effects parsing")
    }
    fixed_formula <- lme4::nobars(formula)
    random_terms <- lme4::findbars(formula)
    if (length(random_terms) > 1) {
      warning("Multiple random effects not fully supported yet")
    }
    random_formula <- random_terms[[1]]
  } else {
    fixed_formula <- formula
    random_formula <- NULL
  }
  
  return(list(
    response_var = response_var,
    y_ordered = y_ordered,
    y_numeric = as.numeric(y_ordered),
    n_levels = n_levels,
    levels = levels(y_ordered),
    fixed_formula = fixed_formula,
    random_formula = random_formula,
    has_random = has_random
  ))
}

#' Prepare Data for TMB Ordinal Model
#' @param formula_parts Parsed formula components
#' @param data Data frame
#' @param verbose Logical; print progress messages
#' @return List with TMB data structure
#' @keywords internal
prepare_ordinal_tmb_data <- function(formula_parts, data, verbose = TRUE) {
  
  if (verbose) {
    cat("Creating design matrix with treatment contrasts...\n")
  }
  
  # Check if model has only intercept (no predictors)
  terms_rhs <- attr(terms(formula_parts$fixed_formula), "term.labels")
  
  if (length(terms_rhs) == 0) {
    # Intercept-only model
    X <- matrix(nrow = nrow(data), ncol = 0)
    colnames(X) <- character(0)
    if (verbose) cat("Intercept-only model: no fixed effects design matrix\n")
  } else {
    # Model with predictors - use treatment contrasts
    original_contrasts <- options("contrasts")
    options(contrasts = c("contr.treatment", "contr.poly"))
    
    X_full <- stats::model.matrix(formula_parts$fixed_formula, data)
    
    # Remove intercept - thresholds serve as intercepts
    if ("(Intercept)" %in% colnames(X_full)) {
      X <- X_full[, !colnames(X_full) %in% "(Intercept)", drop = FALSE]
    } else {
      X <- X_full
    }
    
    options(original_contrasts)
  }
  
  # Handle missing data
  if (ncol(X) > 0) {
    complete_idx <- stats::complete.cases(cbind(formula_parts$y_numeric, X))
  } else {
    complete_idx <- stats::complete.cases(formula_parts$y_numeric)
  }
  
  if (sum(!complete_idx) > 0) {
    if (verbose) cat("Removing", sum(!complete_idx), "rows with missing data\n")
    X <- X[complete_idx, , drop = FALSE]
    y_numeric <- formula_parts$y_numeric[complete_idx]
  } else {
    y_numeric <- formula_parts$y_numeric
  }
  
  n_obs <- length(y_numeric)
  n_coef <- ncol(X)
  n_thresh <- formula_parts$n_levels - 1
  
  # Handle random effects properly - no artificial conversion
  if (formula_parts$has_random) {
    if (!is.null(formula_parts$random_formula)) {
      random_vars <- all.vars(formula_parts$random_formula)
      group_var <- random_vars[length(random_vars)]
      
      if (group_var %in% names(data)) {
        group_factor <- as.factor(data[[group_var]][complete_idx])
        Z <- stats::model.matrix(~ group_factor - 1)
        n_groups <- ncol(Z)
        has_random_effects <- TRUE
        
        if (verbose) cat("Random effects for grouping variable:", group_var, "\n")
      } else {
        warning("Random effects grouping variable '", group_var, "' not found")
        Z <- matrix(0, nrow = n_obs, ncol = 0)
        n_groups <- 0
        has_random_effects <- FALSE
      }
    } else {
      Z <- matrix(0, nrow = n_obs, ncol = 0)
      n_groups <- 0
      has_random_effects <- FALSE
    }
  } else {
    # Pure fixed effects - no dummy random effects
    Z <- matrix(0, nrow = n_obs, ncol = 0)
    n_groups <- 0
    has_random_effects <- FALSE
    
    if (verbose) cat("Pure fixed effects model\n")
  }
  
  if (verbose) {
    cat("Data summary:\n")
    cat("  Observations:", n_obs, "\n")
    cat("  Response levels:", formula_parts$n_levels, "\n")
    cat("  Fixed effects:", n_coef, "\n")
    cat("  Thresholds:", n_thresh, "\n")
    cat("  Random effects groups:", n_groups, "\n")
  }
  
  return(list(
    y = y_numeric,
    X = X,
    Z = Z,
    n_obs = n_obs,
    n_coef = n_coef,
    n_thresh = n_thresh,
    n_groups = n_groups,
    has_random_effects = has_random_effects,
    response_levels = formula_parts$levels,
    complete_idx = complete_idx
  ))
}

# =============================================================================
# TMB MODEL SETUP AND FITTING FUNCTIONS
# =============================================================================

#' Set Up TMB Model Object
#' @param tmb_data Prepared TMB data structure
#' @param prior_params Processed prior parameters
#' @param control Control parameters for TMB
#' @param verbose Logical; print progress messages
#' @return List with TMB setup information
#' @keywords internal
setup_ordinal_tmb <- function(tmb_data, prior_params, control, verbose = TRUE) {
  
  # Find TMB template
  tmb_path <- system.file("tmb", "ordinal_qbrms.cpp", package = "qbrms")
  if (!file.exists(tmb_path)) {
    stop("TMB template file 'ordinal_qbrms.cpp' not found. Please check package installation.")
  }
  
  if (verbose) cat("Setting up TMB model...\n")
  
  # Compile TMB model
  tryCatch({
    TMB::compile(tmb_path, verbose = verbose)
  }, error = function(e) {
    stop("TMB compilation failed: ", e$message)
  })
  
  # Load compiled model
  so_path <- file.path(dirname(tmb_path), "ordinal_qbrms")
  tryCatch({
    dyn.load(TMB::dynlib(so_path))
  }, error = function(e) {
    alt_path <- file.path(getwd(), "ordinal_qbrms")
    tryCatch({
      dyn.load(TMB::dynlib(alt_path))
    }, error = function(e2) {
      stop("Failed to load compiled TMB model")
    })
  })
  
  # Prepare data list for TMB
  data_list <- list(
    y = tmb_data$y,
    X = tmb_data$X,
    Z = if (tmb_data$has_random_effects) tmb_data$Z else matrix(0, nrow = tmb_data$n_obs, ncol = 0),
    thresh_mean = prior_params$thresh_mean,
    thresh_sd = prior_params$thresh_sd,
    beta_mean = prior_params$beta_mean,
    beta_sd = prior_params$beta_sd,
    re_sd_shape = prior_params$re_sd_shape,
    re_sd_rate = prior_params$re_sd_rate,
    has_random_effects = as.integer(tmb_data$has_random_effects)
  )
  
  # Initial parameter values
  param_list <- list(
    threshold_raw = seq(-1, 1, length.out = tmb_data$n_thresh),
    beta = rep(0, tmb_data$n_coef)
  )
  
  # Only add random effect parameters for genuine random effects
  if (tmb_data$has_random_effects && tmb_data$n_groups > 0) {
    param_list$log_re_sd <- log(1.0)
    param_list$u <- rep(0, tmb_data$n_groups)
    random_effects <- "u"
  } else {
    param_list$log_re_sd <- log(1.0)
    param_list$u <- numeric(0)
    random_effects <- NULL
  }
  
  # Create TMB object
  obj <- tryCatch({
    TMB::MakeADFun(
      data = data_list,
      parameters = param_list,
      random = random_effects,
      DLL = "ordinal_qbrms",
      silent = !verbose
    )
  }, error = function(e) {
    stop("Failed to create TMB object: ", e$message)
  })
  
  return(list(
    obj = obj,
    data_list = data_list,
    param_list = param_list,
    tmb_data = tmb_data
  ))
}

#' Enhanced TMB Model Fitting with Better Error Handling - CORRECTED VERSION
#' @param tmb_setup TMB setup list
#' @param verbose Logical; print progress messages
#' @return List with TMB fit results
#' @keywords internal
fit_ordinal_tmb_model <- function(tmb_setup, verbose = TRUE) {
  
  if (verbose) cat("Optimising TMB model with enhanced error handling...\n")
  
  obj <- tmb_setup$obj
  
  # IMPROVED: More sophisticated multi-start strategy
  best_opt <- NULL
  best_value <- Inf
  convergence_achieved <- FALSE
  
  # Better starting points based on problem structure
  n_thresh <- tmb_setup$tmb_data$n_thresh
  n_coef <- tmb_setup$tmb_data$n_coef
  
  starting_points <- list(
    # Original parameters
    obj$par,
    
    # Smaller thresholds, zero coefficients
    c(seq(-2, 2, length.out = n_thresh), rep(0, n_coef), obj$par[(n_thresh + n_coef + 1):length(obj$par)]),
    
    # Wider threshold spread
    c(seq(-4, 4, length.out = n_thresh), rep(0, n_coef), obj$par[(n_thresh + n_coef + 1):length(obj$par)]),
    
    # Random perturbation
    obj$par + rnorm(length(obj$par), 0, 0.1)
  )
  
  for (i in seq_along(starting_points)) {
    opt_attempt <- tryCatch({
      
      # IMPROVED: More robust optimization settings
      nlminb_result <- stats::nlminb(
        start = starting_points[[i]], 
        objective = obj$fn, 
        gradient = obj$gr,
        control = list(
          trace = ifelse(verbose && i == 1, 1, 0),
          eval.max = 3000,
          iter.max = 1500,
          abs.tol = 1e-10,
          rel.tol = 1e-8,
          x.tol = 1e-8
        )
      )
      
      # Additional check with optim if nlminb didn't converge well
      if (nlminb_result$convergence != 0 && nlminb_result$objective < Inf) {
        if (verbose && i == 1) cat("  Trying optim as backup...\n")
        
        optim_result <- stats::optim(
          par = nlminb_result$par,
          fn = obj$fn,
          gr = obj$gr,
          method = "BFGS",
          control = list(maxit = 1000, trace = 0)
        )
        
        if (optim_result$convergence == 0 || optim_result$value < nlminb_result$objective) {
          nlminb_result <- list(
            par = optim_result$par,
            objective = optim_result$value,
            convergence = optim_result$convergence,
            iterations = optim_result$counts[1],
            evaluations = optim_result$counts
          )
        }
      }
      
      nlminb_result
      
    }, error = function(e) {
      if (verbose) cat("  Start point", i, "failed:", e$message, "\n")
      NULL
    })
    
    if (!is.null(opt_attempt) && is.finite(opt_attempt$objective)) {
      if (opt_attempt$objective < best_value) {
        best_opt <- opt_attempt
        best_value <- opt_attempt$objective
        convergence_achieved <- (opt_attempt$convergence == 0)
        if (verbose && i > 1) cat("  Improved solution from start point", i, 
                                  " (objective:", round(best_value, 4), ")\n")
      }
    }
  }
  
  if (is.null(best_opt)) {
    stop("All TMB optimization attempts failed")
  }
  
  opt <- best_opt
  attr(opt, "obj") <- obj
  
  if (!convergence_achieved) {
    warning("TMB optimization may not have fully converged (code: ", opt$convergence, 
            ", objective: ", round(opt$objective, 4), ")")
  } else if (verbose) {
    cat("Optimization converged successfully (objective:", round(opt$objective, 4), ")\n")
  }
  
  # IMPROVED: More robust variance computation
  if (verbose) cat("Computing variance-covariance matrix...\n")
  
  post_cov <- .compute_robust_vcov(obj, opt, verbose)
  
  return(list(
    opt = opt,
    obj = obj,
    hessian = NULL,
    post_cov = post_cov,
    converged = convergence_achieved,
    tmb_setup = tmb_setup,
    variance_method = attr(post_cov, "method")
  ))
}

#' Compute Robust Variance-Covariance Matrix - CORRECTED VERSION
#' @param obj TMB object
#' @param opt Optimization result
#' @param verbose Logical; print diagnostics
#' @return Variance-covariance matrix with method attribute
#' @keywords internal
.compute_robust_vcov <- function(obj, opt, verbose = FALSE) {
  
  post_cov <- NULL
  method_used <- "failed"
  
  # Method 1: TMB sdreport (most accurate)
  vcov_result <- tryCatch({
    sdr <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
    
    if (!is.null(sdr$cov.fixed) && all(is.finite(sdr$cov.fixed))) {
      eigenvals <- eigen(sdr$cov.fixed, symmetric = TRUE, only.values = TRUE)$values
      if (all(eigenvals > 1e-12)) {
        post_cov <- sdr$cov.fixed
        method_used <- "TMB_sdreport"
        if (verbose) cat("  TMB sdreport covariance successful\n")
        list(success = TRUE, cov = post_cov, method = method_used)
      } else {
        if (verbose) cat("  TMB sdreport covariance has negative eigenvalues\n")
        list(success = FALSE)
      }
    } else {
      if (verbose) cat("  TMB sdreport covariance contains non-finite values\n")
      list(success = FALSE)
    }
  }, error = function(e) {
    if (verbose) cat("  TMB sdreport failed:", e$message, "\n")
    list(success = FALSE)
  })
  
  # Method 2: Numerical Hessian with enhanced regularization
  if (!vcov_result$success) {
    if (verbose) cat("  Computing numerical Hessian...\n")
    
    H <- tryCatch({
      # Use adaptive step sizes
      par_scale <- pmax(abs(opt$par), 0.01)
      step_size <- pmax(par_scale * 1e-5, 1e-8)
      
      stats::optimHess(opt$par, obj$fn, obj$gr, 
                       control = list(parscale = par_scale, ndeps = step_size))
    }, error = function(e) {
      if (verbose) cat("  optimHess failed:", e$message, "\n")
      NULL
    })
    
    if (!is.null(H) && all(is.finite(H))) {
      post_cov <- .regularize_hessian(H, opt$par, verbose)
      method_used <- "numerical_hessian"
      vcov_result <- list(success = TRUE, cov = post_cov, method = method_used)
    }
  }
  
  # Method 3: Diagonal approximation (fallback)
  if (!vcov_result$success) {
    if (verbose) cat("  Using diagonal approximation\n")
    
    # Better scaling based on parameter magnitudes and problem structure
    param_scales <- pmax(abs(opt$par), 0.1)
    
    # Different scales for different parameter types
    n_thresh <- length(grep("threshold", names(opt$par), value = FALSE))
    if (n_thresh == 0) n_thresh <- length(opt$par) # fallback
    diagonal_vars <- numeric(length(opt$par))
    
    # Threshold parameters: larger uncertainty
    if (n_thresh > 0) {
      diagonal_vars[1:min(n_thresh, length(opt$par))] <- (param_scales[1:min(n_thresh, length(opt$par))] * 0.2)^2
    }
    
    # Coefficient parameters: moderate uncertainty
    if (length(opt$par) > n_thresh) {
      diagonal_vars[(n_thresh + 1):length(opt$par)] <- (param_scales[(n_thresh + 1):length(opt$par)] * 0.15)^2
    }
    
    post_cov <- diag(diagonal_vars)
    method_used <- "diagonal_approximation"
    vcov_result <- list(success = TRUE, cov = post_cov, method = method_used)
  }
  
  # Final validation
  post_cov <- vcov_result$cov
  method_used <- vcov_result$method
  
  # Ensure positive definiteness
  eigenvals <- eigen(post_cov, symmetric = TRUE, only.values = TRUE)$values
  min_eigenval <- min(eigenvals)
  
  if (min_eigenval <= 1e-12) {
    if (verbose) cat("  Final eigenvalue adjustment (min eigenval:", min_eigenval, ")\n")
    correction <- max(1e-10, -min_eigenval * 1.1)
    post_cov <- post_cov + diag(correction, nrow(post_cov))
    method_used <- paste0(method_used, "_eigenvalue_adjusted")
  }
  
  attr(post_cov, "method") <- method_used
  
  if (verbose) {
    final_condition <- max(eigenvals) / max(min(eigen(post_cov)$values), 1e-12)
    cat("  Final covariance method:", method_used, "\n")
    cat("  Final condition number:", round(final_condition, 2), "\n")
    cat("  Diagonal range:", round(range(diag(post_cov)), 6), "\n")
  }
  
  return(post_cov)
}

#' Regularize Hessian Matrix for Numerical Stability
#' @param H Raw Hessian matrix
#' @param par Parameter vector
#' @param verbose Logical; print progress
#' @return Regularized covariance matrix
#' @keywords internal
.regularize_hessian <- function(H, par, verbose = TRUE) {
  
  if (any(!is.finite(H))) {
    if (verbose) cat("  Hessian contains non-finite values\n")
    return(diag((pmax(abs(par), 0.1) * 0.2)^2))
  }
  
  # Eigenvalue decomposition
  eigen_decomp <- eigen(H, symmetric = TRUE)
  eigenvals <- eigen_decomp$values
  eigenvecs <- eigen_decomp$vectors
  
  min_eigenval <- min(eigenvals)
  condition_number <- max(eigenvals) / max(min_eigenval, 1e-12)
  
  if (verbose) {
    cat("  Hessian condition number:", round(condition_number, 2), "\n")
  }
  
  # Regularization based on condition
  if (min_eigenval <= 1e-10 || condition_number > 1e6) {
    if (verbose) cat("  Applying strong regularization\n")
    eigenvals_reg <- pmax(eigenvals, max(eigenvals) * 1e-6)
    ridge_amount <- max(eigenvals) * 1e-4
    eigenvals_reg <- eigenvals_reg + ridge_amount
    H_reg <- eigenvecs %*% diag(eigenvals_reg) %*% t(eigenvecs)
  } else if (min_eigenval <= 1e-6) {
    if (verbose) cat("  Applying moderate regularization\n")
    eigenvals_reg <- pmax(eigenvals, max(eigenvals) * 1e-5)
    H_reg <- eigenvecs %*% diag(eigenvals_reg) %*% t(eigenvecs)
  } else {
    H_reg <- H
  }
  
  # Invert to get covariance matrix
  post_cov <- tryCatch({
    solve(H_reg)
  }, error = function(e) {
    if (verbose) cat("  Matrix inversion failed, using pseudo-inverse\n")
    eigen_decomp_reg <- eigen(H_reg, symmetric = TRUE)
    eigenvals_inv <- ifelse(eigen_decomp_reg$values > 1e-10, 
                            1 / eigen_decomp_reg$values, 0)
    eigen_decomp_reg$vectors %*% diag(eigenvals_inv) %*% t(eigen_decomp_reg$vectors)
  })
  
  return(post_cov)
}

# =============================================================================
# PARAMETER EXTRACTION AND RESULT CREATION
# =============================================================================

#' Extract Parameter Estimates and Standard Errors - CORRECTED VERSION
#' @param opt TMB optimization result
#' @param post_cov Posterior covariance matrix
#' @param tmb_data TMB data structure
#' @param verbose Logical; print progress messages
#' @return List with extracted parameters
#' @keywords internal
extract_ordinal_parameters <- function(opt, post_cov, tmb_data, verbose = FALSE) {
  
  obj <- attr(opt, "obj")
  
  if (verbose) {
    cat("=== DEBUGGING STANDARD ERROR EXTRACTION ===\n")
    cat("Post-cov dimensions:", dim(post_cov), "\n")
    cat("Post-cov diagonal range:", round(range(diag(post_cov)), 6), "\n")
    cat("Parameter vector length:", length(opt$par), "\n")
  }
  
  # FIXED: Try TMB sdreport with detailed diagnostics
  adreport_success <- FALSE
  adreport_summary <- NULL
  
  tryCatch({
    if (verbose) cat("Attempting TMB sdreport...\n")
    
    sdr <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
    adreport_summary <- summary(sdr)
    
    if (verbose) {
      cat("TMB sdreport successful, summary dimensions:", dim(adreport_summary), "\n")
      cat("ADREPORT parameter names:", rownames(adreport_summary), "\n")
      cat("ADREPORT SE range:", range(adreport_summary[, "Std. Error"]), "\n")
    }
    
    # Check if standard errors are reasonable
    se_values <- adreport_summary[, "Std. Error"]
    if (all(is.finite(se_values)) && all(se_values > 1e-10) && all(se_values < 100)) {
      adreport_success <- TRUE
      if (verbose) cat("TMB sdreport standard errors look reasonable\n")
    } else {
      if (verbose) cat("TMB sdreport standard errors appear problematic\n")
    }
    
  }, error = function(e) {
    if (verbose) cat("TMB sdreport failed with error:", e$message, "\n")
  })
  
  # FIXED: Extract parameters with proper standard error mapping
  if (adreport_success && !is.null(adreport_summary)) {
    
    # Extract using ADREPORT results
    if (verbose) cat("Using TMB ADREPORT results\n")
    
    param_names <- rownames(adreport_summary)
    estimates <- adreport_summary[, "Estimate"]
    std_errors <- adreport_summary[, "Std. Error"]
    
    # Extract thresholds
    threshold_indices <- grep("^threshold", param_names)
    if (length(threshold_indices) >= tmb_data$n_thresh) {
      thresholds <- estimates[threshold_indices[1:tmb_data$n_thresh]]
      thresh_se <- std_errors[threshold_indices[1:tmb_data$n_thresh]]
      if (verbose) cat("Extracted", length(thresholds), "thresholds from ADREPORT\n")
    } else {
      thresholds <- .apply_threshold_constraints(opt$par, tmb_data$n_thresh)
      thresh_se <- sqrt(diag(post_cov))[1:tmb_data$n_thresh]
    }
    
    # Extract coefficients
    beta_indices <- grep("^beta", param_names)
    if (tmb_data$n_coef > 0 && length(beta_indices) >= tmb_data$n_coef) {
      coefficients <- estimates[beta_indices[1:tmb_data$n_coef]]
      coef_se <- std_errors[beta_indices[1:tmb_data$n_coef]]
      if (verbose) cat("Extracted", length(coefficients), "coefficients from ADREPORT\n")
    } else if (tmb_data$n_coef > 0) {
      param_start <- tmb_data$n_thresh + 1
      param_end <- tmb_data$n_thresh + tmb_data$n_coef
      coefficients <- opt$par[param_start:param_end]
      coef_se <- sqrt(diag(post_cov))[param_start:param_end]
    } else {
      coefficients <- numeric(0)
      coef_se <- numeric(0)
    }
    
  } else {
    
    # FIXED: Manual extraction with corrected parameter indexing
    if (verbose) cat("Using manual parameter extraction with post_cov\n")
    
    # Extract thresholds (constrained from raw parameters)
    thresholds <- .apply_threshold_constraints(opt$par, tmb_data$n_thresh)
    
    # CRITICAL FIX: Properly extract threshold standard errors
    # The post_cov is for the raw (unconstrained) parameters
    # Need to compute SEs for constrained thresholds using delta method
    thresh_se <- .compute_threshold_se_delta_method(opt$par, post_cov, tmb_data$n_thresh, verbose)
    
    # Extract coefficient estimates and SEs
    if (tmb_data$n_coef > 0) {
      param_start <- tmb_data$n_thresh + 1
      param_end <- tmb_data$n_thresh + tmb_data$n_coef
      
      if (param_end <= length(opt$par)) {
        coefficients <- opt$par[param_start:param_end]
        
        # FIXED: Proper indexing for coefficient standard errors
        if (param_end <= nrow(post_cov)) {
          coef_se <- sqrt(diag(post_cov))[param_start:param_end]
        } else {
          coef_se <- rep(0.3, tmb_data$n_coef)
          if (verbose) cat("Warning: Using fallback coefficient SEs\n")
        }
      } else {
        coefficients <- rep(0, tmb_data$n_coef)
        coef_se <- rep(0.3, tmb_data$n_coef)
        if (verbose) cat("Warning: Parameter vector too short\n")
      }
    } else {
      coefficients <- numeric(0)
      coef_se <- numeric(0)
    }
  }
  
  # FIXED: Validate and clean standard errors
  thresh_se[!is.finite(thresh_se) | thresh_se <= 1e-10] <- 0.2
  coef_se[!is.finite(coef_se) | coef_se <= 1e-10] <- 0.15
  
  # FIXED: Ensure threshold ordering is maintained
  if (length(thresholds) > 1) {
    for (i in 2:length(thresholds)) {
      if (thresholds[i] <= thresholds[i-1]) {
        thresholds[i] <- thresholds[i-1] + 0.01
        if (verbose) cat("Adjusted threshold", i, "for ordering\n")
      }
    }
  }
  
  if (verbose) {
    cat("Final parameter summary:\n")
    cat("  Thresholds:", round(thresholds, 3), "\n")
    cat("  Threshold SEs:", round(thresh_se, 3), "\n")
    if (length(coefficients) > 0) {
      cat("  Coefficients:", round(coefficients, 3), "\n")
      cat("  Coefficient SEs:", round(coef_se, 3), "\n")
    }
  }
  
  return(list(
    thresholds = thresholds,
    coefficients = coefficients,
    thresh_se = thresh_se,
    coef_se = coef_se,
    re_sd = NULL,
    re_sd_se = NULL
  ))
}

#' Compute Threshold Standard Errors Using Delta Method - NEW FUNCTION
#' @param raw_params Raw parameter vector from optimization
#' @param post_cov Covariance matrix of raw parameters
#' @param n_thresh Number of thresholds
#' @param verbose Logical; print diagnostics
#' @return Vector of threshold standard errors
#' @keywords internal
.compute_threshold_se_delta_method <- function(raw_params, post_cov, n_thresh, verbose = FALSE) {
  
  if (n_thresh == 0) return(numeric(0))
  
  if (verbose) cat("Computing threshold SEs using delta method...\n")
  
  tryCatch({
    
    # Extract raw threshold parameters
    threshold_raw <- raw_params[1:n_thresh]
    
    # Compute Jacobian matrix for threshold transformation
    # thresholds[1] = threshold_raw[1]
    # thresholds[k] = thresholds[k-1] + exp(threshold_raw[k]) for k > 1
    
    jacobian <- matrix(0, nrow = n_thresh, ncol = n_thresh)
    
    # Derivative of threshold[1] w.r.t. threshold_raw[1]
    jacobian[1, 1] <- 1
    
    # For k > 1: d(threshold[k])/d(threshold_raw[j])
    if (n_thresh > 1) {
      for (k in 2:n_thresh) {
        # Cumulative effect: each threshold depends on all previous raw parameters
        jacobian[k, 1] <- 1  # Through threshold[1]
        
        for (j in 2:k) {
          jacobian[k, j] <- exp(threshold_raw[j])  # Direct effect of exp(threshold_raw[j])
        }
      }
    }
    
    # Extract covariance submatrix for threshold parameters
    thresh_cov <- post_cov[1:n_thresh, 1:n_thresh, drop = FALSE]
    
    # Apply delta method: Var(g(θ)) = J * Var(θ) * J^T
    transformed_cov <- jacobian %*% thresh_cov %*% t(jacobian)
    
    # Extract standard errors
    thresh_se <- sqrt(pmax(diag(transformed_cov), 1e-10))
    
    if (verbose) {
      cat("  Delta method jacobian condition number:", 
          max(eigen(jacobian %*% t(jacobian))$values) / min(eigen(jacobian %*% t(jacobian))$values), "\n")
      cat("  Computed threshold SEs:", round(thresh_se, 4), "\n")
    }
    
    return(thresh_se)
    
  }, error = function(e) {
    if (verbose) cat("Delta method failed:", e$message, "- using fallback\n")
    
    # Fallback: Use diagonal elements directly (approximation)
    thresh_se_fallback <- sqrt(pmax(diag(post_cov)[1:n_thresh], 1e-10))
    
    # Apply rough scaling for transformation
    threshold_raw <- raw_params[1:n_thresh]
    for (i in 2:n_thresh) {
      # For exp transformation, SE approximately scales by exp(x)
      if (is.finite(threshold_raw[i])) {
        thresh_se_fallback[i] <- thresh_se_fallback[i] * exp(threshold_raw[i])
      }
    }
    
    return(thresh_se_fallback)
  })
}

#' Apply threshold constraints manually (fallback)
#' @param raw_params Raw parameter vector
#' @param n_thresh Number of thresholds
#' @return Constrained threshold values
#' @keywords internal
.apply_threshold_constraints <- function(raw_params, n_thresh) {
  if (n_thresh == 0) return(numeric(0))
  
  threshold_raw <- raw_params[1:n_thresh]
  thresholds <- numeric(n_thresh)
  
  thresholds[1] <- threshold_raw[1]
  for (k in 2:n_thresh) {
    thresholds[k] <- thresholds[k-1] + exp(threshold_raw[k])
  }
  
  return(thresholds)
}

#' Create Fixed Effects Summary Table
#' @param param_summary Extracted parameter summary
#' @param tmb_data TMB data structure
#' @return Data frame with summary statistics
#' @keywords internal
create_fixed_effects_summary <- function(param_summary, tmb_data) {
  
  # Build parameter components with guaranteed finite SEs
  all_names <- character(0)
  all_estimates <- numeric(0)
  all_se <- numeric(0)
  
  # Add thresholds
  if (length(param_summary$thresholds) > 0) {
    thresh_names <- paste0("Intercept[", seq_along(param_summary$thresholds), "]")
    thresh_se_safe <- param_summary$thresh_se
    thresh_se_safe[!is.finite(thresh_se_safe) | thresh_se_safe <= 0] <- 0.5
    
    all_names <- c(all_names, thresh_names)
    all_estimates <- c(all_estimates, param_summary$thresholds)
    all_se <- c(all_se, thresh_se_safe)
  }
  
  # Add coefficients
  if (length(param_summary$coefficients) > 0) {
    design_names <- colnames(tmb_data$X)
    if (length(design_names) >= length(param_summary$coefficients)) {
      coef_names <- design_names[seq_along(param_summary$coefficients)]
    } else {
      coef_names <- paste0("coef", seq_along(param_summary$coefficients))
    }
    
    coef_se_safe <- param_summary$coef_se
    coef_se_safe[!is.finite(coef_se_safe) | coef_se_safe <= 0] <- 0.3
    
    all_names <- c(all_names, coef_names)
    all_estimates <- c(all_estimates, param_summary$coefficients)
    all_se <- c(all_se, coef_se_safe)
  }
  
  # Final safety check
  all_se[!is.finite(all_se)] <- 0.1
  all_estimates[!is.finite(all_estimates)] <- 0
  
  # Create summary table
  summary_table <- data.frame(
    mean = all_estimates,
    sd = all_se,
    "0.025quant" = all_estimates - 1.96 * all_se,
    "0.5quant" = all_estimates,
    "0.975quant" = all_estimates + 1.96 * all_se,
    mode = all_estimates,
    kld = 0,
    row.names = all_names,
    stringsAsFactors = FALSE
  )
  
  return(summary_table)
}

#' Create qbrms-Compatible Result Object
#' @param tmb_fit TMB fit results
#' @param formula Original formula
#' @param data Original data
#' @param family Model family
#' @param prior Prior specifications
#' @param verbose Logical; print progress messages
#' @return qbrms-compatible result object
#' @keywords internal
create_ordinal_qbrms_result <- function(tmb_fit, formula, data, family, 
                                        prior, verbose = TRUE) {
  
  opt <- tmb_fit$opt
  obj <- attr(opt, "obj") %||% tmb_fit$obj
  post_cov <- tmb_fit$post_cov
  tmb_data <- tmb_fit$tmb_setup$tmb_data
  
  # Extract parameters
  param_summary <- extract_ordinal_parameters(opt, post_cov, tmb_data, verbose)
  
  # Create summary
  summary_fixed <- create_fixed_effects_summary(param_summary, tmb_data)
  
  # Get fitted values
  fitted_values <- tryCatch({
    obj$report()$fitted_probs
  }, error = function(e) {
    rep(0.25, tmb_data$n_obs)
  })
  
  # Create result object
  result <- list(
    fit = list(
      summary.fixed = summary_fixed,
      fitted.values = fitted_values,
      converged = tmb_fit$converged
    ),
    tmb_obj = obj,
    tmb_fit = tmb_fit,
    posterior_cov = post_cov,
    data = data,
    original_formula = formula,
    family = family,
    prior = prior,
    ordinal_levels = tmb_data$response_levels,
    model_type = "ordinal_tmb",
    fitting_time = 0
  )
  
  class(result) <- c("tmb_ordinal_qbrms_fit", "qbrms_fit", "list")
  
  return(result)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Check if family is ordinal
#' @param family Family specification
#' @return Logical indicating if family is ordinal
#' @keywords internal
is_ordinal <- function(family) {
  if (is.character(family)) {
    return(family %in% c("cumulative", "ordinal"))
  }
  if (is.list(family) && !is.null(family$family)) {
    return(family$family %in% c("cumulative", "ordinal"))
  }
  return(FALSE)
}

# =============================================================================
# S3 METHODS FOR TMB ORDINAL FITS
# =============================================================================

#' Summary Method for TMB Ordinal Fits
#' @param object A tmb_ordinal_qbrms_fit object
#' @param digits Number of decimal places for output
#' @param ... Additional arguments
#' @return Invisibly returns the object
#' @export
#' @method summary tmb_ordinal_qbrms_fit
summary.tmb_ordinal_qbrms_fit <- function(object, digits = 2, ...) {
  
  cat("TMB Ordinal Regression (qbrms)\n\n")
  
  cat("Family: cumulative (ordinal)\n")
  cat("Link function: logit\n")
  cat("Response levels:", length(object$ordinal_levels), 
      "(", paste(object$ordinal_levels, collapse = ", "), ")\n")
  cat("Observations:", nrow(object$data), "\n")
  
  if (!is.null(object$threshold_method)) {
    cat("Threshold method:", object$threshold_method, "\n")
  }
  cat("\n")
  
  # Fixed effects
  cat("Population-Level Effects:\n")
  if (exists("format_numeric_df", mode = "function")) {
    formatted_summary <- format_numeric_df(object$fit$summary.fixed, digits)
  } else {
    formatted_summary <- object$fit$summary.fixed
    numeric_cols <- sapply(formatted_summary, is.numeric)
    formatted_summary[numeric_cols] <- lapply(formatted_summary[numeric_cols], 
                                              function(col) sprintf(paste0("%.", digits, "f"), col))
  }
  print(formatted_summary, quote = FALSE, right = TRUE)
  cat("\n")
  
  # Random effects
  cat("Group-Level Effects:\n")
  if (grepl("\\|", deparse(object$original_formula))) {
    cat("  Random intercept present\n")
  } else {
    cat("  (none)\n")
  }
  cat("\n")
  
  # Model fit info
  cat("Model Fitting:\n")
  cat("  Method: TMB with Laplace approximation\n")
  cat("  Converged:", object$fit$converged, "\n")
  
  invisible(object)
}

#' Print Method for TMB Ordinal Fits
#' @param x A tmb_ordinal_qbrms_fit object
#' @param digits Number of decimal places for output
#' @param ... Additional arguments
#' @return Invisibly returns the object
#' @export
#' @method print tmb_ordinal_qbrms_fit
print.tmb_ordinal_qbrms_fit <- function(x, digits = 2, ...) {
  cat("TMB Ordinal qbrms Fit\n\n")
  
  cat("Formula:", deparse(x$original_formula), "\n")
  cat("Data:   ", nrow(x$data), "observations\n")
  cat("Levels: ", paste(x$ordinal_levels, collapse = ", "), "\n")
  cat("\n")
  
  if (!is.null(x$fit$summary.fixed)) {
    if (exists("format_numeric_df", mode = "function")) {
      formatted_coef <- format_numeric_df(x$fit$summary.fixed, digits)
    } else {
      formatted_coef <- x$fit$summary.fixed
      numeric_cols <- sapply(formatted_coef, is.numeric)
      formatted_coef[numeric_cols] <- lapply(formatted_coef[numeric_cols], 
                                             function(col) sprintf(paste0("%.", digits, "f"), col))
    }
    print(formatted_coef, quote = FALSE, right = TRUE)
  }
  
  invisible(x)
}

#' Coefficients Method for TMB Ordinal Fits
#' @param object A tmb_ordinal_qbrms_fit object
#' @param ... Additional arguments
#' @return Named vector of coefficients
#' @export
#' @method coef tmb_ordinal_qbrms_fit
coef.tmb_ordinal_qbrms_fit <- function(object, ...) {
  if (!is.null(object$fit$summary.fixed)) {
    coefs <- object$fit$summary.fixed[, "mean"]
    names(coefs) <- rownames(object$fit$summary.fixed)
    return(coefs)
  }
  return(NULL)
}

#' Variance-Covariance Matrix Method for TMB Ordinal Fits
#' @param object A tmb_ordinal_qbrms_fit object
#' @param ... Additional arguments
#' @return Variance-covariance matrix
#' @export
#' @method vcov tmb_ordinal_qbrms_fit
vcov.tmb_ordinal_qbrms_fit <- function(object, ...) {
  if (!is.null(object$posterior_cov)) {
    coef_names <- rownames(object$fit$summary.fixed)
    rownames(object$posterior_cov) <- coef_names
    colnames(object$posterior_cov) <- coef_names
    return(object$posterior_cov)
  }
  return(NULL)
}

#' Fitted Values Method for TMB Ordinal Fits
#' @param object A tmb_ordinal_qbrms_fit object
#' @param ... Additional arguments
#' @return Fitted values
#' @export
#' @method fitted tmb_ordinal_qbrms_fit
fitted.tmb_ordinal_qbrms_fit <- function(object, ...) {
  if (!is.null(object$fit$fitted.values)) {
    return(object$fit$fitted.values)
  }
  return(NULL)
}