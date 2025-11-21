# =============================================================================
# R/ordinal_plots.R
# Ordinal posterior predictive checks + conditional effects for TMB ordinal fits
# =============================================================================

#' Ordinal Plots and Posterior Predictive Checks
#'
#' @description
#' File-level imports and utilities for ordinal model plotting functions.
#' 
#' @name ordinal_plots
#' @keywords internal
#' @import ggplot2
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "category","proportion","predicted_prob","observed_freq","bin_size",
    "effect_value","probability","mean_prob","lower_prob","upper_prob",
    "observation","type","sample__","bin_center","mean_observed","n_in_bin",
    "pred_mean", "pred_lower", "pred_upper", "observed", "predicted"
  ))
}

# ---- helper: build RHS design matrix using training contrasts ----------------
#' @keywords internal
#' @noRd
.qbrms_build_X <- function(object, newdata) {
  dat <- object$data
  if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0L) {
    if (!is.null(object$fit) && is.data.frame(object$fit$frame) && nrow(object$fit$frame) > 0L) {
      dat <- object$fit$frame
    } else {
      stop("Training data not found in object$data or object$fit$frame.")
    }
  }
  
  # Build fixed-effects terms (strips any "(...|...)" random-effect chunks)
  # Note: assumes .qbrms_terms_fixed is available in the package namespace
  tr_terms <- .qbrms_terms_fixed(object, dat)
  
  X_train  <- stats::model.matrix(tr_terms, dat)
  contrs   <- attr(X_train, "contrasts")
  
  stats::model.matrix(tr_terms, newdata, contrasts.arg = contrs)
}

# ---- parameter extraction helpers --------------------------------------------
#' @keywords internal
#' @noRd
.qbrms_ordinal_params <- function(object) {
  stopifnot(!is.null(object$fit$summary.fixed))
  S  <- object$fit$summary.fixed
  nm <- rownames(S)
  th_idx <- grepl("^Intercept\\[", nm)
  thresholds <- setNames(S$mean[th_idx], nm[th_idx])
  betas <- setNames(S$mean[!th_idx], nm[!th_idx])
  ord_thr <- order(as.integer(gsub("\\D", "", names(thresholds))))
  list(thresholds = unname(thresholds[ord_thr]), betas = betas)
}

#' @keywords internal
#' @noRd
.qbrms_make_pred_grid_ord <- function(object, effect, resolution = 100, conditions = NULL) {
  dat <- object$data
  if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0L) {
    if (!is.null(object$fit) && is.data.frame(object$fit$frame) && nrow(object$fit$frame) > 0L) {
      dat <- object$fit$frame
    } else {
      stop("Training data not found in object$data or object$fit$frame.")
    }
  }
  if (!effect %in% names(dat)) {
    stop("Effect '", effect, "' not found in training data.")
  }
  
  if (is.numeric(dat[[effect]])) {
    r <- range(dat[[effect]], na.rm = TRUE)
    if (!is.finite(r[1]) || !is.finite(r[2])) stop("Cannot determine range for numeric effect '", effect, "'.")
    grid <- data.frame(setNames(list(seq(r[1], r[2], length.out = resolution)), effect))
  } else {
    levs <- levels(as.factor(dat[[effect]]))
    grid <- data.frame(setNames(list(factor(levs, levels = levs)), effect))
  }
  
  # Hold other RHS variables fixed (mean/mode) unless supplied in `conditions`
  # Note: assumes .qbrms__remove_random_effects is available in package namespace
  rhs_all <- all.vars(.qbrms__remove_random_effects(object$original_formula))
  rhs_vars <- setdiff(rhs_all[-1], effect)
  for (v in rhs_vars) {
    if (!is.null(conditions) && v %in% names(conditions)) {
      grid[[v]] <- conditions[[v]]
    } else if (v %in% names(dat)) {
      if (is.numeric(dat[[v]])) {
        grid[[v]] <- mean(dat[[v]], na.rm = TRUE)
      } else {
        levs <- levels(as.factor(dat[[v]]))
        mode_val <- names(sort(table(dat[[v]]), decreasing = TRUE))[1]
        grid[[v]] <- factor(mode_val, levels = levs)
      }
    }
  }
  
  for (v in names(grid)) {
    if (v %in% names(dat) && is.factor(dat[[v]])) {
      grid[[v]] <- factor(grid[[v]], levels = levels(dat[[v]]))
    }
  }
  grid
}

#' @keywords internal
#' @noRd
.qbrms_eta_ord <- function(object, grid, betas) {
  X <- .qbrms_build_X(object, grid)
  if ("(Intercept)" %in% colnames(X))
    X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  b <- betas[colnames(X)]; b[is.na(b)] <- 0
  as.numeric(X %*% b)
}

#' @keywords internal
#' @noRd
.qbrms_conditional_probs_ord <- function(object, effect, resolution = 100, conditions = NULL) {
  pars <- .qbrms_ordinal_params(object)
  grid <- .qbrms_make_pred_grid_ord(object, effect, resolution, conditions)
  eta  <- .qbrms_eta_ord(object, grid, pars$betas)
  
  if (any(!is.finite(eta)))
    stop("Linear predictor contains non-finite values; check factor levels in prediction grid.")
  
  resp <- all.vars(object$original_formula)[1]
  y    <- object$data[[resp]]
  levs <- if (is.factor(y) || is.ordered(y)) levels(y) else sort(unique(y))
  K    <- length(levs)
  if (length(pars$thresholds) != (K - 1L))
    stop("Mismatch: number of thresholds != K-1")
  
  P  <- matrix(NA_real_, nrow = nrow(grid), ncol = K)
  for (i in seq_len(nrow(grid))) {
    cum <- numeric(K)
    if (K > 1L) cum[1:(K-1)] <- stats::plogis(pars$thresholds - eta[i])
    cum[K] <- 1
    P[i, 1] <- cum[1]
    for (k in 2:K) P[i, k] <- cum[k] - cum[k-1]
    P[i, ] <- pmax(1e-12, P[i, ]); P[i, ] <- P[i, ] / sum(P[i, ])
  }
  
  list(prob = P, grid = grid, levels = levs)
}

# ---- safe tabulation of predicted categories per draw ------------------------
#' @keywords internal
#' @noRd
.qbrms_tabulate_pred_counts <- function(pred_categories, n_cats) {
  nd <- nrow(pred_categories)
  out <- matrix(0, nrow = n_cats, ncol = nd)
  for (d in seq_len(nd)) out[, d] <- tabulate(pred_categories[d, ], nbins = n_cats)
  out # rows = categories, cols = draws
}

# =============================================================================
# pp_check() for ordinal TMB fits
# =============================================================================

#' Posterior predictive checks for TMB ordinal models
#'
#' @param object A fitted TMB ordinal qbrms model object
#' @param type Character; type of posterior predictive check
#' @param ndraws Integer; number of posterior draws to use
#' @param seed Random seed for reproducibility.
#' @param newdata Optional data frame for predictions. If NULL, uses original data.
#' @param prob Probability mass for credible intervals (default 0.95).
#' @param ... Additional arguments passed to methods.
#'
#' @return A ggplot object showing the posterior predictive check
#' @export
#' @method pp_check tmb_ordinal_qbrms_fit
pp_check.tmb_ordinal_qbrms_fit <- function(object, type = "bars", ndraws = 100,
                                           seed = NULL, newdata = NULL, prob = 0.9, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for ordinal posterior predictive checks")
  }
  if (!is.null(seed)) set.seed(seed)
  
  allowed_types <- c("bars", "rootogram", "ribbon", "intervals", "calibration")
  if (!is.character(type) || length(type) != 1L || !(type %in% allowed_types)) {
    stop("For ordinal models, type must be one of: ",
         paste(allowed_types, collapse = ", "))
  }
  
  if (is.null(newdata)) newdata <- object$data
  pred_results <- .generate_ordinal_posterior_predictions(object, newdata, ndraws)
  
  plot_result <- switch(
    type,
    "bars"       = .pp_check_ordinal_bars(pred_results),
    "rootogram"  = .pp_check_ordinal_rootogram(pred_results),
    "ribbon"     = .pp_check_ordinal_ribbon(pred_results, prob),
    "intervals"  = .pp_check_ordinal_intervals(pred_results, prob),
    "calibration"= .pp_check_ordinal_calibration(pred_results)
  )
  
  class(plot_result) <- c("qbrms_ordinal_pp_check", class(plot_result))
  plot_result
}

# ---- draw predictions for pp_check -------------------------------------------

#' @keywords internal
#' @noRd
.generate_ordinal_parameter_draws <- function(object, ndraws) {
  if (is.null(object$fit$summary.fixed)) stop("Cannot extract parameter estimates")
  sf <- object$fit$summary.fixed
  nm <- rownames(sf); m <- sf[, "mean"]; s <- sf[, "sd"]
  
  th_idx <- grepl("^Intercept\\[", nm)
  co_idx <- !th_idx
  
  th_m <- as.numeric(m[th_idx]); th_s <- as.numeric(s[th_idx])
  co_m <- as.numeric(m[co_idx]);  co_s <- as.numeric(s[co_idx])
  names(co_m) <- nm[co_idx]
  
  if (length(th_m)) {
    th_draws <- matrix(NA_real_, nrow = ndraws, ncol = length(th_m))
    for (j in seq_along(th_m)) th_draws[, j] <- stats::rnorm(ndraws, th_m[j], th_s[j])
  } else th_draws <- matrix(0, nrow = ndraws, ncol = 0)
  
  if (length(co_m)) {
    co_draws <- matrix(NA_real_, nrow = ndraws, ncol = length(co_m))
    for (j in seq_along(co_m)) co_draws[, j] <- stats::rnorm(ndraws, co_m[j], co_s[j])
    colnames(co_draws) <- names(co_m)
  } else co_draws <- matrix(0, nrow = ndraws, ncol = 0)
  
  list(thresholds = th_draws, coefs = co_draws)
}

#' @keywords internal
#' @noRd
.generate_ordinal_p <- function(thr, eta_vec) {
  K <- length(thr) + 1L
  n <- length(eta_vec)
  P <- matrix(NA_real_, nrow = n, ncol = K)
  for (i in seq_len(n)) {
    cum <- numeric(K)
    if (K > 1L) cum[1:(K-1)] <- stats::plogis(thr - eta_vec[i])
    cum[K] <- 1
    P[i, 1] <- cum[1]
    for (k in 2:K) P[i, k] <- cum[k] - cum[k-1]
    P[i, ] <- pmax(1e-12, P[i, ]); P[i, ] <- P[i, ] / sum(P[i, ])
  }
  P
}

#' @keywords internal
#' @noRd
.generate_ordinal_posterior_predictions <- function(object, newdata, ndraws) {
  response_var <- all.vars(object$original_formula)[1]
  y_obs <- newdata[[response_var]]
  
  if (is.factor(y_obs) || is.ordered(y_obs)) {
    level_labels <- levels(y_obs)
    y_obs_num    <- as.integer(y_obs)
  } else {
    y_obs_num    <- as.integer(as.factor(y_obs))
    level_labels <- levels(factor(y_obs))
  }
  
  n_obs  <- length(y_obs_num)
  n_cats <- length(level_labels)
  
  draws <- .generate_ordinal_parameter_draws(object, ndraws)
  
  X_full <- .qbrms_build_X(object, newdata)
  if ("(Intercept)" %in% colnames(X_full))
    X_full <- X_full[, setdiff(colnames(X_full), "(Intercept)"), drop = FALSE]
  
  if (ncol(X_full) > 0 && ncol(draws$coefs) > 0) {
    common <- intersect(colnames(X_full), colnames(draws$coefs))
    B_all  <- matrix(0, nrow = ndraws, ncol = ncol(X_full)); colnames(B_all) <- colnames(X_full)
    if (length(common)) B_all[, common] <- draws$coefs[, common, drop = FALSE]
    eta_all <- B_all %*% t(X_full) # ndraws x n_obs
  } else {
    eta_all <- matrix(0, nrow = ndraws, ncol = n_obs)
  }
  
  if (ncol(draws$thresholds) != (n_cats - 1L))
    stop("Mismatch between number of thresholds and response levels (K-1).")
  
  pred_probs_array <- array(NA_real_, dim = c(ndraws, n_obs, n_cats))
  pred_categories  <- matrix(NA_integer_, nrow = ndraws, ncol = n_obs)
  
  for (d in seq_len(ndraws)) {
    P <- .generate_ordinal_p(draws$thresholds[d, ], eta_all[d, ])
    pred_probs_array[d, , ] <- P
    for (i in seq_len(n_obs)) {
      pred_categories[d, i] <- sample.int(n_cats, size = 1, prob = P[i, ])
    }
  }
  
  list(
    y_obs           = y_obs_num,
    pred_probs      = pred_probs_array,
    pred_categories = pred_categories,
    level_labels    = level_labels,
    n_cats          = n_cats,
    ndraws          = ndraws,
    n_obs           = n_obs
  )
}

# ---- plotting helpers for ordinal pp_check -----------------------------------

#' @keywords internal
#' @noRd
.pp_check_ordinal_bars <- function(pred_results) {
  n_cats <- pred_results$n_cats
  levs   <- pred_results$level_labels
  
  obs_counts <- as.numeric(table(factor(pred_results$y_obs, levels = seq_len(n_cats))))
  obs_props  <- obs_counts / sum(obs_counts)
  
  pred_counts_draws <- .qbrms_tabulate_pred_counts(pred_results$pred_categories, n_cats) # n_cats x ndraws
  pred_props_mean   <- rowMeans(pred_counts_draws) / pred_results$n_obs
  pred_props_q025   <- apply(pred_counts_draws, 1, stats::quantile, 0.025) / pred_results$n_obs
  pred_props_q975   <- apply(pred_counts_draws, 1, stats::quantile, 0.975) / pred_results$n_obs
  
  df_obs <- data.frame(
    category   = factor(levs, levels = levs),
    proportion = obs_props,
    type       = "Observed",
    lower      = obs_props,
    upper      = obs_props
  )
  df_pred <- data.frame(
    category   = factor(levs, levels = levs),
    proportion = pred_props_mean,
    type       = "Predicted",
    lower      = pred_props_q025,
    upper      = pred_props_q975
  )
  plot_data <- rbind(df_obs, df_pred)
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = category, y = proportion, fill = type)) +
    ggplot2::geom_col(position = "dodge", alpha = 0.85) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.22
    ) +
    ggplot2::scale_fill_manual(values = c("Observed" = "#2C3E50", "Predicted" = "#3498DB"), name = NULL) +
    ggplot2::labs(
      title = "Posterior Predictive Check: Category Proportions",
      x = "Category", y = "Proportion"
    ) +
    ggplot2::theme_minimal()
}

#' @keywords internal
#' @noRd
.pp_check_ordinal_rootogram <- function(pred_results) {
  n_cats <- pred_results$n_cats
  levs   <- pred_results$level_labels
  
  obs_counts <- as.numeric(table(factor(pred_results$y_obs, levels = seq_len(n_cats))))
  pred_counts_draws <- .qbrms_tabulate_pred_counts(pred_results$pred_categories, n_cats)
  pred_counts_mean <- rowMeans(pred_counts_draws)
  
  plot_data <- data.frame(
    category = factor(seq_len(n_cats), labels = levs),
    observed = sqrt(obs_counts),
    predicted = sqrt(pred_counts_mean)
  )
  plot_data$residual <- plot_data$observed - plot_data$predicted
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = category)) +
    ggplot2::geom_col(ggplot2::aes(y = predicted), alpha = 0.6, fill = "#3498DB") +
    ggplot2::geom_point(ggplot2::aes(y = observed),
                        size = 3, shape = 21, fill = "#2C3E50", colour = "white") +
    ggplot2::geom_segment(ggplot2::aes(xend = category, y = predicted, yend = observed),
                          colour = "#2C3E50", linewidth = 0.8) +
    ggplot2::labs(title = "Rootogram (sqrt counts): Observed vs Predicted",
                  x = "Category", y = "sqrt(Count)") +
    ggplot2::theme_minimal()
}

#' @keywords internal
#' @noRd
.pp_check_ordinal_ribbon <- function(pred_results, prob) {
  alpha <- 1 - prob
  mean_probs   <- apply(pred_results$pred_probs, c(2, 3), mean)  # n_obs x K
  lower_probs  <- apply(pred_results$pred_probs, c(2, 3), stats::quantile, alpha/2)
  upper_probs  <- apply(pred_results$pred_probs, c(2, 3), stats::quantile, 1 - alpha/2)
  
  df <- do.call(rbind, lapply(seq_len(pred_results$n_cats), function(k) {
    data.frame(
      observation = seq_len(pred_results$n_obs),
      category    = factor(pred_results$level_labels[k], levels = pred_results$level_labels),
      mean_prob   = mean_probs[, k],
      lower_prob  = lower_probs[, k],
      upper_prob  = upper_probs[, k],
      observed    = as.numeric(pred_results$y_obs == k)
    )
  }))
  
  ggplot2::ggplot(df, ggplot2::aes(x = observation)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower_prob, ymax = upper_prob, fill = category),
      alpha = 0.3, na.rm = TRUE
    ) +
    ggplot2::geom_line(ggplot2::aes(y = mean_prob, color = category), linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(y = observed, color = category),
                        size = 0.8, alpha = 0.7, na.rm = TRUE) +
    ggplot2::facet_wrap(~ category, scales = "free_y") +
    ggplot2::labs(
      title = paste0("Prediction Uncertainty by Category (", prob * 100, "%)"),
      x = "Observation", y = "Probability"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

#' @keywords internal
#' @noRd
.pp_check_ordinal_intervals <- function(pred_results, prob) {
  alpha <- 1 - prob
  n_cats <- pred_results$n_cats
  levs   <- pred_results$level_labels
  
  obs_counts <- as.numeric(table(factor(pred_results$y_obs, levels = seq_len(n_cats))))
  
  pred_counts_draws <- .qbrms_tabulate_pred_counts(pred_results$pred_categories, n_cats)
  pred_counts_mean  <- rowMeans(pred_counts_draws)
  pred_counts_lower <- apply(pred_counts_draws, 1, stats::quantile, alpha/2)
  pred_counts_upper <- apply(pred_counts_draws, 1, stats::quantile, 1 - alpha/2)
  
  plot_data <- data.frame(
    category   = factor(seq_len(n_cats), labels = levs),
    observed   = obs_counts,
    pred_mean  = pred_counts_mean,
    pred_lower = pred_counts_lower,
    pred_upper = pred_counts_upper
  )
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = category)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(y = pred_mean, ymin = pred_lower, ymax = pred_upper),
      color = "#3498DB", size = 1, linewidth = 1.2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = observed),
      color = "#E74C3C", size = 3, shape = 18
    ) +
    ggplot2::labs(
      title = paste0("Prediction Intervals vs Observed Counts (", prob * 100, "%)"),
      x = "Category", y = "Count"
    ) +
    ggplot2::theme_minimal()
}

#' @keywords internal
#' @noRd
.pp_check_ordinal_calibration <- function(pred_results) {
  breaks <- seq(0, 1, by = 0.1)
  bin_centers_all <- (head(breaks, -1) + tail(breaks, -1)) / 2
  
  rows <- list()
  for (cat in seq_len(pred_results$n_cats)) {
    # FIX: Use colMeans for drawing average probability per observation
    # pred_probs is [ndraws, n_obs, n_cats] -> we want [n_obs]
    p <- colMeans(pred_results$pred_probs[, , cat]) 
    o <- as.numeric(pred_results$y_obs == cat)
    
    bins <- cut(p, breaks = breaks, include.lowest = TRUE, right = TRUE)
    # For each bin: mean observed & count
    means <- tapply(o, bins, function(x) if (length(x)) mean(x) else NA_real_)
    counts <- tapply(o, bins, length)
    
    # Keep only bins that actually occurred
    idx <- which(!is.na(means))
    if (length(idx)) {
      rows[[length(rows) + 1L]] <- data.frame(
        category      = factor(pred_results$level_labels[cat], levels = pred_results$level_labels),
        bin_center    = bin_centers_all[idx],
        mean_observed = as.numeric(means[idx]),
        n_in_bin      = as.numeric(counts[idx])
      )
    }
  }
  
  df <- if (length(rows)) do.call(rbind, rows) else {
    stop("No calibration bins contained data; try fewer categories or different data.")
  }
  
  ggplot2::ggplot(df, ggplot2::aes(x = bin_center, y = mean_observed, color = category, group = category)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_point() +
    ggplot2::geom_line() +  # will quietly drop for categories with single point
    ggplot2::labs(title = "Calibration Plot by Category",
                  x = "Predicted probability (bin center)", y = "Observed frequency") +
    ggplot2::theme_minimal()
}

# =============================================================================
# Conditional effects for TMB ordinal fits
# =============================================================================

#' Conditional Effects for TMB Ordinal Models
#'
#' @param object A tmb_ordinal_qbrms_fit object
#' @param effects Character vector of effect names (defaults to auto-detected)
#' @param prob Confidence level
#' @param ndraws Number of draws
#' @param spaghetti Logical
#' @param n_points Number of points for continuous predictors
#' @param plot Logical, whether to return plots
#' @param at Named list of conditioning values
#' @param seed Random seed
#' @param conditions Ordinal-specific conditions (for backwards compatibility)
#' @param categorical Whether to show categorical plot (for backwards compatibility)
#' @param resolution Grid resolution (for backwards compatibility)
#' @param ... Additional arguments
#'
#' @return List of conditional effects
#' @export
#' @method conditional_effects tmb_ordinal_qbrms_fit
conditional_effects.tmb_ordinal_qbrms_fit <- function(object, 
                                                      effects = NULL,
                                                      prob = 0.95,
                                                      ndraws = 100,
                                                      spaghetti = FALSE,
                                                      n_points = 100,
                                                      plot = TRUE,
                                                      at = list(),
                                                      seed = NULL,
                                                      # Backwards compatibility parameters
                                                      conditions = NULL,
                                                      categorical = TRUE,
                                                      resolution = NULL,
                                                      ...) {
  
  # Handle backwards compatibility - map old parameters to new ones
  if (!is.null(conditions)) at <- conditions
  if (!is.null(resolution)) n_points <- resolution
  
  # Input validation
  if (!inherits(object, "tmb_ordinal_qbrms_fit")) {
    stop("object must be a tmb_ordinal_qbrms_fit object")
  }
  
  if (is.null(object$data) || is.null(object$original_formula)) {
    stop("Model object missing required data or formula components")
  }
  
  # Auto-detect effects if not provided
  if (is.null(effects)) {
    # Choose data source
    dat <- object$data
    if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0L) {
      if (!is.null(object$fit) && is.data.frame(object$fit$frame) && nrow(object$fit$frame) > 0L) {
        dat <- object$fit$frame
      } else {
        stop("Training data not found in object$data or object$fit$frame.")
      }
    }
    # Strip random effects and find RHS predictors
    f_fixed <- .qbrms__remove_random_effects(object$original_formula)
    allv    <- all.vars(f_fixed)
    if (length(allv) < 2L) stop("Formula appears to have no predictors.")
    resp    <- allv[1]
    cand    <- setdiff(allv, resp)
    cand    <- cand[cand %in% names(dat)]
    num     <- cand[vapply(dat[cand], is.numeric, TRUE)]
    effects <- if (length(num)) num[1] else if (length(cand)) cand[1] else stop("No usable predictors found to plot.")
  }
  
  
  out <- list()
  for (eff in effects) {
    cp <- .qbrms_conditional_probs_ord(object, eff, n_points, at)
    df <- data.frame(
      effect_value = rep(cp$grid[[eff]], times = length(cp$levels)),
      category     = factor(rep(cp$levels, each = nrow(cp$grid)), levels = cp$levels),
      probability  = as.vector(cp$prob)
    )
    df <- df[order(df$category, df$effect_value), ]
    
    if (categorical && plot) {
      p <- ggplot2::ggplot(df, ggplot2::aes(effect_value, probability, fill = category, group = category)) +
        ggplot2::geom_area(alpha = 0.7, position = "stack") +
        ggplot2::labs(title = paste("Conditional Effects:", eff),
                      subtitle = "Predicted category probabilities",
                      x = eff, y = "Probability", fill = "Category") +
        ggplot2::theme_minimal()
    } else if (plot) {
      p <- ggplot2::ggplot(df, ggplot2::aes(effect_value, probability, color = category)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::facet_wrap(~ category, scales = "free_y") +
        ggplot2::labs(title = paste("Conditional Effects:", eff),
                      subtitle = "Predicted category probabilities by category",
                      x = eff, y = "Probability") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
    } else {
      p <- NULL
    }
    
    result <- if (plot) {
      structure(list(plot = p, data = df, effect_var = eff),
                class = c("qbrms_ordinal_conditional_effect", "list"))
    } else {
      df
    }
    
    out[[eff]] <- result
  }
  
  structure(out,
            class = c("qbrms_ordinal_conditional_effects", "list"),
            effects = effects, categorical = categorical)
}

#' @export
#' @method plot qbrms_ordinal_conditional_effects
plot.qbrms_ordinal_conditional_effects <- function(x, ask = FALSE, ...) {
  if (length(x) == 1) return(x[[1]]$plot)
  lapply(x, `[[`, "plot")
}

# Tiny alias to match earlier references in older code
#' @keywords internal
.create_ordinal_prediction_grid_safe <- function(object, effect, conditions, resolution) {
  .qbrms_make_pred_grid_ord(object, effect, resolution, conditions)
}