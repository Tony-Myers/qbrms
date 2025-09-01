# =============================================================================
# R/conditional_effects.R - SIMPLIFIED ROBUST VERSION
# =============================================================================

#' Conditional Effects for qbrms Models
#'
#' @description
#' Generic function for computing conditional effects
#'
#' @param object Model object
#' @param ... Additional arguments
#'
#' @return Object of class "brms_conditional_effects"
#' @export
conditional_effects <- function(object, ...) {
  UseMethod("conditional_effects")
}

#' Conditional Effects for qbrms Models
#'
#' @description
#' Compute conditional effects for qbrms model objects, similar to brms.
#'
#' @param object A qbrms model object
#' @param effects Character vector specifying effects to plot. If NULL, plots all numeric predictors.
#' @param conditions Named list of values for conditional variables
#' @param resolution Integer specifying number of prediction points
#' @param prob Probability for credible intervals (default 0.95)
#' @param method Method for generating predictions ("fitted" is default)
#' @param re_formula Formula for random effects (currently ignored)
#' @param spaghetti Logical, whether to add spaghetti plots
#' @param ndraws Number of posterior draws for spaghetti plots
#' @param ribbon Logical, whether to show ribbon (confidence band)
#' @param ... Additional arguments
#'
#' @return Object of class "brms_conditional_effects"
#'
#' @export
conditional_effects.qbrms_fit <- function(object, effects = NULL,
                                          conditions = NULL,
                                          resolution = 100,
                                          prob = 0.95,
                                          method = "fitted",
                                          re_formula = NULL,
                                          spaghetti = FALSE,
                                          ndraws = 100,
                                          ribbon = TRUE,
                                          ...) {
  
  predictors <- get_predictor_variables_simple(object)
  
  if (is.null(effects)) {
    effects <- predictors$numeric_vars
    if (length(effects) == 0) {
      effects <- predictors$all_vars[1]
    }
  }
  
  results <- list()
  
  for (effect in effects) {
    cat("Processing effect:", effect, "\n")
    
    tryCatch({
      pred_data <- create_prediction_grid_simple(effect, object$data, predictors,
                                                 conditions, resolution)
      
      if (spaghetti) {
        predictions <- generate_spaghetti_predictions_simple(object, pred_data, prob, ndraws)
      } else {
        predictions <- generate_predictions_simple(object, pred_data, prob)
      }
      
      effect_data <- format_brms_conditional_effects(predictions, effect, prob, pred_data, spaghetti)
      results[[effect]] <- effect_data
      
    }, error = function(e) {
      cat("Error processing effect", effect, ":", e$message, "\n")
      # Create fallback data
      n_pred <- if (effect %in% predictors$numeric_vars) resolution else length(unique(object$data[[effect]]))
      fallback_pred <- list(
        estimate__ = rep(0, n_pred),
        se__ = rep(1, n_pred),
        lower__ = rep(-2, n_pred),
        upper__ = rep(2, n_pred)
      )
      pred_data_fallback <- data.frame(x = seq_len(n_pred))
      names(pred_data_fallback) <- effect
      effect_data <- format_brms_conditional_effects(fallback_pred, effect, prob, pred_data_fallback, FALSE)
      results[[effect]] <- effect_data
    })
  }
  
  class(results) <- c("brms_conditional_effects", "list")
  attr(results, "effects") <- effects
  attr(results, "prob") <- prob
  attr(results, "method") <- method
  attr(results, "spaghetti") <- spaghetti
  attr(results, "ndraws") <- if(spaghetti) ndraws else NULL
  attr(results, "ribbon") <- ribbon
  
  return(results)
}

#' @export
conditional_effects.bru_brms_fit <- conditional_effects.qbrms_fit

# Simplified helper function to get predictor variables
get_predictor_variables_simple <- function(object) {
  # Simple extraction that avoids complex conditionals
  data <- object$data
  
  # Try to get variables from the original formula
  vars <- c()
  
  tryCatch({
    if (is.character(object$original_formula)) {
      formula_obj <- as.formula(object$original_formula)
    } else {
      formula_obj <- object$original_formula
    }
    
    # Get all variables except response
    all_vars <- all.vars(formula_obj)
    vars <- all_vars[-1]  # Remove response variable
    
    # Remove common grouping variables
    vars <- vars[!vars %in% c("group_id", "athlete_id", "site_id", "trials")]
    
  }, error = function(e) {
    # Fallback: use numeric columns from data
    vars <- names(data)[sapply(data, is.numeric)]
  })
  
  # Categorise variables
  numeric_vars <- c()
  categorical_vars <- c()
  
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

# Simplified helper function to create prediction grid
create_prediction_grid_simple <- function(effect, data, predictors, conditions, resolution) {
  
  # Create effect variable values
  if (effect %in% predictors$numeric_vars) {
    effect_range <- range(data[[effect]], na.rm = TRUE)
    effect_values <- seq(effect_range[1], effect_range[2], length.out = resolution)
  } else {
    # FIXED: For categorical, get ALL unique levels, properly sorted
    effect_values <- sort(unique(data[[effect]]))
    cat("Categorical levels for", effect, ":", paste(effect_values, collapse = ", "), "\n")
  }
  
  # Create prediction data
  pred_data <- data.frame(effect_values, stringsAsFactors = FALSE)
  names(pred_data) <- effect
  
  # FIXED: For categorical effects, ensure proper factor structure
  if (!effect %in% predictors$numeric_vars && is.factor(data[[effect]])) {
    pred_data[[effect]] <- factor(pred_data[[effect]], levels = levels(data[[effect]]))
    cat("Made", effect, "a factor with levels:", paste(levels(pred_data[[effect]]), collapse = ", "), "\n")
  }
  
  # Add other predictors at their means/modes
  for (var in predictors$all_vars[predictors$all_vars != effect]) {
    if (var %in% names(data)) {
      if (var %in% predictors$numeric_vars) {
        constant_value <- mean(data[[var]], na.rm = TRUE)
        pred_data[[var]] <- constant_value
        cat("Setting", var, "to constant value:", constant_value, "\n")
      } else {
        # Use most frequent level for factors
        mode_value <- names(sort(table(data[[var]]), decreasing = TRUE))[1]
        pred_data[[var]] <- mode_value
        cat("Setting", var, "to mode value:", mode_value, "\n")
      }
    }
  }
  
  # Apply conditions if specified (overwrites defaults)
  if (!is.null(conditions)) {
    for (var in names(conditions)) {
      if (var %in% names(pred_data)) {
        pred_data[[var]] <- conditions[[var]]
        cat("Applied condition:", var, "=", conditions[[var]], "\n")
      }
    }
  }
  
  # Ensure factors have proper levels
  for (var in names(pred_data)) {
    if (var %in% names(data)) {
      if (is.factor(data[[var]])) {
        pred_data[[var]] <- factor(pred_data[[var]], levels = levels(data[[var]]))
      }
    }
  }
  
  cat("Final prediction data (first 6 rows):\n")
  print(head(pred_data))
  
  return(pred_data)
}

# Simplified prediction generation function with FIXED categorical handling
generate_predictions_simple <- function(object, pred_data, prob) {
  
  n_pred <- nrow(pred_data)
  alpha <- 1 - prob
  upper_q <- 1 - alpha / 2
  
  # Get family info safely
  family_name <- "gaussian"
  if (!is.null(object$family)) {
    if (is.list(object$family)) {
      family_name <- object$family$family
    } else {
      family_name <- as.character(object$family)[1]
    }
  }
  
  tryCatch({
    # Get fixed effects
    fixed_effects <- object$fit$summary.fixed
    coef_names <- rownames(fixed_effects)
    
    cat("Available coefficients:", paste(coef_names, collapse = ", "), "\n")
    cat("Prediction data columns:", paste(names(pred_data), collapse = ", "), "\n")
    
    # FIXED: Create proper design matrix handling factors correctly
    # Remove any response variables first
    clean_pred_data <- pred_data
    
    # For binomial models, remove response variables  
    if (family_name == "binomial") {
      # Remove common response variable names
      response_vars <- c("successes", "failures", "y", "response")
      for (rv in response_vars) {
        if (rv %in% names(clean_pred_data)) {
          clean_pred_data[[rv]] <- NULL
        }
      }
    }
    
    # Create model matrix using the ORIGINAL data factor levels
    # This is crucial for categorical variables to work properly
    for (var in names(clean_pred_data)) {
      if (var %in% names(object$data) && is.factor(object$data[[var]])) {
        # Ensure the prediction data factor has the same levels as original data
        clean_pred_data[[var]] <- factor(clean_pred_data[[var]], 
                                         levels = levels(object$data[[var]]))
        cat("Set factor levels for", var, ":", paste(levels(clean_pred_data[[var]]), collapse = ", "), "\n")
      }
    }
    
    # Create design matrix - this will automatically handle factor encoding
    if (ncol(clean_pred_data) > 0) {
      design_formula <- as.formula(paste("~", paste(names(clean_pred_data), collapse = " + ")))
      cat("Design matrix formula:", deparse(design_formula), "\n")
      
      X <- model.matrix(design_formula, data = clean_pred_data)
      cat("Design matrix dimensions:", dim(X), "\n")
      cat("Design matrix columns:", paste(colnames(X), collapse = ", "), "\n")
    } else {
      # Intercept only
      X <- matrix(1, nrow = n_pred, ncol = 1)
      colnames(X) <- "(Intercept)"
    }
    
    # FIXED: Match coefficients more carefully
    available_coefs <- intersect(coef_names, colnames(X))
    cat("Matched coefficients:", paste(available_coefs, collapse = ", "), "\n")
    
    if (length(available_coefs) > 0) {
      beta <- as.numeric(fixed_effects[available_coefs, "mean"])
      beta_sd <- as.numeric(fixed_effects[available_coefs, "sd"])
      
      cat("Using beta values:", paste(round(beta, 4), collapse = ", "), "\n")
      
      # Calculate predictions - this should now work correctly for categorical variables
      X_matched <- X[, available_coefs, drop = FALSE]
      linear_pred <- as.numeric(X_matched %*% beta)
      se_pred <- sqrt(rowSums((X_matched^2) * rep(beta_sd^2, each = n_pred)))
      
      cat("Linear predictions range:", paste(round(range(linear_pred), 4), collapse = " to "), "\n")
      
      # Apply link function
      if (family_name == "binomial") {
        estimate <- plogis(linear_pred)
        lower <- plogis(linear_pred - qnorm(upper_q) * se_pred)
        upper <- plogis(linear_pred + qnorm(upper_q) * se_pred)
        se_estimate <- (upper - lower) / (2 * qnorm(upper_q))
      } else if (family_name == "poisson") {
        estimate <- exp(linear_pred)
        lower <- exp(linear_pred - qnorm(upper_q) * se_pred)
        upper <- exp(linear_pred + qnorm(upper_q) * se_pred)
        se_estimate <- (upper - lower) / (2 * qnorm(upper_q))
      } else {
        estimate <- linear_pred
        se_estimate <- se_pred
        lower <- estimate - qnorm(upper_q) * se_estimate
        upper <- estimate + qnorm(upper_q) * se_estimate
      }
      
      cat("Final estimates range:", paste(round(range(estimate), 4), collapse = " to "), "\n")
      
    } else {
      # No matching coefficients - use fallbacks
      cat("No matching coefficients found, using fallbacks\n")
      if (family_name == "binomial") {
        estimate <- rep(0.5, n_pred)
        lower <- rep(0.3, n_pred)
        upper <- rep(0.7, n_pred)
      } else {
        estimate <- rep(0, n_pred)
        lower <- rep(-2, n_pred) 
        upper <- rep(2, n_pred)
      }
      se_estimate <- rep(1, n_pred)
    }
    
    return(list(
      estimate__ = as.numeric(estimate),
      se__ = as.numeric(se_estimate),
      lower__ = as.numeric(lower),
      upper__ = as.numeric(upper)
    ))
    
  }, error = function(e) {
    cat("Error in generate_predictions_simple:", e$message, "\n")
    
    # Ultimate fallback
    if (family_name == "binomial") {
      estimate <- rep(0.5, n_pred)
      lower <- rep(0.3, n_pred)
      upper <- rep(0.7, n_pred)
    } else {
      estimate <- rep(0, n_pred)
      lower <- rep(-2, n_pred)
      upper <- rep(2, n_pred)
    }
    
    return(list(
      estimate__ = estimate,
      se__ = rep(1, n_pred),
      lower__ = lower,
      upper__ = upper
    ))
  })
}

# RESTORED: Proper spaghetti predictions with original functionality
generate_spaghetti_predictions_simple <- function(object, pred_data, prob, ndraws) {
  # Check for MASS package
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("MASS package required for multivariate normal sampling in spaghetti plots")
  }
  
  # Get the base predictions first
  base_predictions <- generate_predictions_simple(object, pred_data, prob)
  
  # Get INLA fit and extract variance-covariance matrix
  inla_fit <- object$fit
  fixed_effects <- inla_fit$summary.fixed$mean
  fixed_names <- rownames(inla_fit$summary.fixed)
  
  cat("Generating spaghetti predictions with", ndraws, "draws\n")
  
  tryCatch({
    # Create design matrix for predictions (matching the simple prediction approach)
    predictor_vars <- intersect(names(pred_data), 
                                gsub("(Intercept)", "", fixed_names, fixed = TRUE))
    
    if (length(predictor_vars) > 0) {
      design_formula <- as.formula(paste("~", paste(predictor_vars, collapse = " + ")))
      X <- model.matrix(design_formula, data = pred_data)
    } else {
      X <- matrix(1, nrow = nrow(pred_data), ncol = 1)
      colnames(X) <- "(Intercept)"
    }
    
    # Ensure all fixed effects are represented
    missing_cols <- setdiff(fixed_names, colnames(X))
    for (col in missing_cols) {
      if (col != "(Intercept)") {
        X <- cbind(X, 0)
        colnames(X)[ncol(X)] <- col
      }
    }
    X <- X[, fixed_names, drop = FALSE]
    
    # Get variance-covariance matrix
    vcov_matrix <- tryCatch({
      if (inherits(inla_fit, "quantile_inla")) {
        inla_fit$vcov
      } else {
        vcov(inla_fit)
      }
    }, error = function(e) {
      # Fallback to diagonal approximation
      diag(pmax(inla_fit$summary.fixed$sd^2, 1e-8))
    })
    
    # Validate vcov matrix
    if (any(is.na(vcov_matrix)) || any(diag(vcov_matrix) <= 0)) {
      warning("Issues with variance-covariance matrix, using diagonal approximation")
      vcov_matrix <- diag(pmax(inla_fit$summary.fixed$sd^2, 1e-8))
    }
    
    # Generate proper posterior draws using MASS::mvrnorm
    set.seed(123)  # For reproducibility
    posterior_draws <- MASS::mvrnorm(n = ndraws,
                                     mu = fixed_effects,
                                     Sigma = vcov_matrix)
    
    # Calculate predictions for each draw
    all_predictions <- matrix(NA, nrow = nrow(pred_data), ncol = ndraws)
    
    # Get family info for link function
    family_name <- "gaussian"
    if (!is.null(object$family)) {
      if (is.list(object$family)) {
        family_name <- object$family$family
      } else {
        family_name <- as.character(object$family)[1]
      }
    }
    
    for (i in 1:ndraws) {
      beta_draw <- posterior_draws[i, ]
      linear_pred <- X %*% beta_draw
      
      # Apply link function if needed
      if (family_name == "binomial") {
        fitted_pred <- plogis(linear_pred)
      } else if (family_name == "poisson") {
        fitted_pred <- exp(linear_pred)
      } else {
        fitted_pred <- linear_pred
      }
      
      all_predictions[, i] <- as.numeric(fitted_pred)
    }
    
    # Create spaghetti data
    effect_var <- names(pred_data)[1]  # First column is effect variable
    spaghetti_data <- data.frame()
    
    for (i in 1:ndraws) {
      draw_data <- data.frame(
        pred_data[effect_var],
        estimate__ = all_predictions[, i],
        sample__ = paste0("draw_", i)
      )
      names(draw_data)[1] <- effect_var
      spaghetti_data <- rbind(spaghetti_data, draw_data)
    }
    
    # Add spaghetti data as attribute
    attr(base_predictions, "spaghetti") <- spaghetti_data
    cat("Spaghetti data created with", nrow(spaghetti_data), "rows\n")
    
    return(base_predictions)
    
  }, error = function(e) {
    warning("Error generating spaghetti predictions: ", e$message)
    cat("Falling back to base predictions only\n")
    return(base_predictions)
  })
}

# Keep the original helper functions
format_brms_conditional_effects <- function(predictions, effect_var, prob, pred_data, spaghetti = FALSE) {
  result <- data.frame(
    effect__ = pred_data[[effect_var]], 
    estimate__ = predictions$estimate__,
    se__ = predictions$se__,
    lower__ = predictions$lower__,
    upper__ = predictions$upper__
  )
  
  names(result)[1] <- effect_var
  
  attr(result, "effect") <- effect_var
  attr(result, "prob") <- prob
  attr(result, "spaghetti") <- spaghetti
  
  if (spaghetti && !is.null(attr(predictions, "spaghetti"))) {
    attr(result, "spaghetti") <- attr(predictions, "spaghetti")
  }
  
  class(result) <- c("data.frame")
  return(result)
}

# =============================================================================
# R/plotting_methods.R
# =============================================================================

#' Plot Conditional Effects
#'
#' @description
#' Plot method for conditional effects from qbrms models.
#'
#' @param x Object of class "brms_conditional_effects"
#' @param ask Logical, whether to ask before showing each plot
#' @param spaghetti_args List of arguments for spaghetti plot appearance
#' @param line_args List of arguments for main line appearance
#' @param ribbon Logical, whether to show confidence ribbon
#' @param ribbon_args List of arguments for ribbon appearance
#' @param ... Additional arguments passed to ggplot2
#'
#' @return List of ggplot objects (returned invisibly)
#'
#' @export
plot.brms_conditional_effects <- function(x, ask = TRUE,
                                          spaghetti_args = list(),
                                          line_args = list(),
                                          ribbon = TRUE,
                                          ribbon_args = list(),
                                          ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting conditional effects")
  }
  
  plots <- list()
  
  for (i in seq_along(x)) {
    effect_name <- names(x)[i]
    effect_data <- x[[i]]
    
    # Check if this specific effect has spaghetti data
    spaghetti_data <- attr(effect_data, "spaghetti")
    has_spaghetti_this_effect <- !is.null(spaghetti_data) &&
      is.data.frame(spaghetti_data) &&
      nrow(spaghetti_data) > 0
    
    # Debug output
    cat("Effect:", effect_name, "Has spaghetti:", has_spaghetti_this_effect, "\n")
    
    # Create base plot
    p <- ggplot2::ggplot(effect_data, ggplot2::aes(x = .data[[effect_name]], y = .data[["estimate__"]]))
    
    # FIXED: Only add ribbon if NO spaghetti data AND ribbon is TRUE
    if (ribbon && !has_spaghetti_this_effect) {
      cat("Adding ribbon for effect:", effect_name, "\n")
      default_ribbon_args <- list(
        ggplot2::aes(ymin = .data[["lower__"]], ymax = .data[["upper__"]]),
        alpha = 0.3,
        fill = "blue"
      )
      
      ribbon_args_final <- modifyList(default_ribbon_args, ribbon_args)
      p <- p + do.call(ggplot2::geom_ribbon, ribbon_args_final)
    } else if (has_spaghetti_this_effect) {
      cat("Skipping ribbon for spaghetti effect:", effect_name, "\n")
    }
    
    # Add spaghetti lines if present
    if (has_spaghetti_this_effect) {
      cat("Adding spaghetti lines for effect:", effect_name, "\n")
      default_spaghetti_args <- list(
        data = spaghetti_data,
        ggplot2::aes(x = .data[[effect_name]], y = .data[["estimate__"]],
                     group = .data[["sample__"]]),
        alpha = 0.2,
        linewidth = 0.3,
        color = "steelblue"
      )
      
      spaghetti_args_final <- modifyList(default_spaghetti_args, spaghetti_args)
      p <- p + do.call(ggplot2::geom_line, spaghetti_args_final)
    }
    
    # Add main prediction line with categorical handling
    if (is.factor(effect_data[[effect_name]]) || is.character(effect_data[[effect_name]])) {
      # For categorical variables: use points with error bars
      
      # Add error bars first (so they appear behind points)
      if (!has_spaghetti_this_effect) {
        default_errorbar_args <- list(
          ggplot2::aes(ymin = .data[["lower__"]], ymax = .data[["upper__"]]),
          width = 0.2,
          alpha = 0.8,
          color = "blue"
        )
        
        # Allow customization of error bars via ribbon_args for consistency
        errorbar_args_final <- modifyList(default_errorbar_args, ribbon_args)
        p <- p + do.call(ggplot2::geom_errorbar, errorbar_args_final)
        cat("Adding error bars for categorical effect:", effect_name, "\n")
      }
      
      # Add points
      default_point_args <- list(
        color = if(has_spaghetti_this_effect) "darkblue" else "blue",
        size = if(has_spaghetti_this_effect) 1.5 else 2,
        alpha = if(has_spaghetti_this_effect) 0.8 else 1
      )
      
      point_args_final <- modifyList(default_point_args, line_args)
      p <- p + do.call(ggplot2::geom_point, point_args_final)
      
    } else {
      # For continuous variables: use lines with your existing logic
      default_line_args <- list(
        color = if(has_spaghetti_this_effect) "darkblue" else "blue",
        linewidth = if(has_spaghetti_this_effect) 0 else 1,  # Keep your linewidth = 0 for spaghetti
        alpha = if(has_spaghetti_this_effect) 0.8 else 1
      )
      
      line_args_final <- modifyList(default_line_args, line_args)
      p <- p + do.call(ggplot2::geom_line, line_args_final)
    }
    
    # Add labels and theme
    p <- p +
      ggplot2::labs(
        title = paste("Conditional Effect of", effect_name),
        x = effect_name,
        y = "Predicted Response"
      ) +
      ggplot2::theme_minimal()
    
    plots[[effect_name]] <- p
    print(p)
    
    if (ask && i < length(x)) {
      readline("Press Enter for next plot...")
    }
  }
  
  invisible(plots)
}