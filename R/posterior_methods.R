# R/posterior_methods.R

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("kind", "idx", "value", "y", "ybar", "d"))
}

# -----------------------------------------------------------------------------
# Generic
# -----------------------------------------------------------------------------

#' Posterior and prior predictive checks
#'
#' Create posterior or prior predictive diagnostic plots for fitted qbrms models.
#'
#' @param object A model object.
#' @param ... Additional arguments passed to methods.
#' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}

# -----------------------------------------------------------------------------
# Internal helpers (shared)
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.extract_family_name_simple <- function(family) {
  if (is.null(family)) return("unknown")
  if (is.character(family)) return(as.character(family)[1L])
  if (is.list(family) && !is.null(family$family)) return(as.character(family$family)[1L])
  "unknown"
}

#' Core pp_check implementation (continuous families)
#' @keywords internal
#' @noRd
.pp_check_core <- function(
    object,
    type = "dens_overlay",
    ndraws = 5000,
    seed = NULL,
    is_prior = inherits(object, "qbrms_prior_fit"),
    show_observed = FALSE,
    ...
) {
  allowed <- c("dens_overlay", "hist", "scatter", "scatter_avg")
  if (is.null(type) || !is.character(type) || length(type) != 1L || !(type %in% allowed)) {
    stop("Unsupported pp_check type: ", type, call. = FALSE)
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # observed response (robust)
  response_var <- tryCatch(all.vars(object$original_formula)[1], error = function(e) NULL)
  y_raw <- if (!is.null(response_var) && !is.null(object$data)) object$data[[response_var]] else NULL
  y <- NULL
  idx <- NULL
  if (!is.null(y_raw)) {
    y_num <- if (is.ordered(y_raw)) {
      as.numeric(y_raw)
    } else if (is.factor(y_raw)) {
      suppressWarnings(as.numeric(as.character(y_raw)))
    } else {
      as.numeric(y_raw)
    }
    idx <- is.finite(y_num)
    if (any(idx)) y <- y_num[idx]
  }
  
  # replicated draws
  if (is_prior) {
    yrep <- object$prior_samples
    if (!is.null(yrep) && is.matrix(yrep) && nrow(yrep) > ndraws) {
      yrep <- yrep[seq_len(ndraws), , drop = FALSE]
    }
  } else {
    yrep <- tryCatch(generate_posterior_predictive_samples(object, ndraws), error = function(e) NULL)
  }
  
  # align to observed mask
  if (!is.null(idx) && is.logical(idx) && length(idx) > 0L &&
      !is.null(yrep) && is.matrix(yrep) && ncol(yrep) == length(idx)) {
    yrep <- yrep[, idx, drop = FALSE]
  }
  
  # fallback if empty
  if (is.null(yrep) || !is.matrix(yrep) || nrow(yrep) < 1L || ncol(yrep) < 1L) {
    n_obs <- if (!is.null(y) && length(y) > 0L) length(y) else 100L
    n_draw <- max(1L, ndraws)
    yrep <- matrix(stats::rnorm(n_draw * n_obs, 0, 1), nrow = n_draw, ncol = n_obs)
  }
  if (is.null(y) || length(y) < 1L) y <- as.numeric(yrep[1, ])
  
  title_text <- if (is_prior) {
    t <- "Prior Predictive Check"
    if (show_observed) paste(t, "(with observed data)") else t
  } else {
    "Posterior Predictive Check"
  }
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    col_obs <- "#000000"
    col_rep <- "#b3cde0"
    col_obs_prior <- "#E31A1C"
    rep_size <- 0.45
    rep_alpha <- 0.35
    p <- switch(
      type,
      "dens_overlay" = {
        k <- min(50L, nrow(yrep))
        rep_list <- lapply(seq_len(k), function(i) {
          di <- stats::density(as.numeric(yrep[i, ]))
          data.frame(x = di$x, d = di$y, draw = i)
        })
        df_rep <- do.call(rbind, rep_list)
        show_obs_data <- if (is_prior) isTRUE(show_observed) else TRUE
        plot <- ggplot2::ggplot() +
          ggplot2::geom_line(
            data = df_rep,
            ggplot2::aes(x = x, y = d, group = draw),
            linewidth = rep_size, alpha = rep_alpha, colour = col_rep
          )
        if (show_obs_data) {
          d_obs <- stats::density(as.numeric(y))
          df_obs <- data.frame(x = d_obs$x, d = d_obs$y)
          obs_colour <- if (is_prior && show_observed) col_obs_prior else col_obs
          plot <- plot + ggplot2::geom_line(
            data = df_obs,
            ggplot2::aes(x = x, y = d),
            linewidth = 1.15, colour = obs_colour
          )
        }
        plot + 
          ggplot2::labs(title = title_text, x = NULL, y = NULL) +
          ggplot2::theme_classic()
      },
      "hist" = {
        n_panels <- 10L
        k <- max(1L, min(nrow(yrep), n_panels - 1L))
        df_y <- data.frame(
          value = as.numeric(y),
          panel = factor("y", levels = c("y", paste0("yrep[", seq_len(k), "]"))),
          source = factor("y", levels = c("y", "yrep"))
        )
        rep_list <- lapply(seq_len(k), function(i) {
          data.frame(
            value = as.numeric(yrep[i, ]),
            panel = paste0("yrep[", i, "]"),
            source = "yrep"
          )
        })
        df_rep <- do.call(rbind, rep_list)
        df_rep$panel <- factor(df_rep$panel, levels = levels(df_y$panel))
        df_rep$source <- factor(df_rep$source, levels = levels(df_y$source))
        df_all <- rbind(df_y, df_rep)
        col_rep_line <- "#6f93ad"
        col_y_fill <- col_obs
        col_y_line <- "#1a1a1a"
        ggplot2::ggplot(df_all, ggplot2::aes(x = value)) +
          ggplot2::geom_histogram(
            ggplot2::aes(fill = source, colour = source),
            bins = 30,
            alpha = 0.9,
            linewidth = 0.25
          ) +
          ggplot2::facet_wrap(~ panel, ncol = 5, scales = "free_y") +
          ggplot2::scale_fill_manual(values = c("y" = col_y_fill, "yrep" = col_rep), name = NULL) +
          ggplot2::scale_colour_manual(values = c("y" = col_y_line, "yrep" = col_rep_line), guide = "none") +
          ggplot2::labs(title = title_text, x = NULL, y = NULL) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(
            strip.text = ggplot2::element_text(size = 8),
            legend.position = "right",
            panel.grid.minor = ggplot2::element_blank()
          )
      },
      "scatter" = {
        n_panels <- 8L
        y_num <- as.numeric(y)
        stopifnot(length(y_num) >= 1L, is.matrix(yrep), ncol(yrep) == length(y_num))
        total_draws <- nrow(yrep)
        k <- max(1L, min(total_draws, n_panels))
        draw_ids <- if (k < total_draws) sample.int(total_draws, k) else seq_len(k)
        panel_levels <- paste0("yrep[", seq_len(k), "]")
        rep_data <- do.call(rbind, lapply(seq_len(k), function(j) {
          i <- draw_ids[j]
          data.frame(
            x = y_num,
            y = as.numeric(yrep[i, ]),
            panel = factor(paste0("yrep[", j, "]"), levels = panel_levels)
          )
        }))
        xlim <- range(rep_data$x, na.rm = TRUE)
        ylim <- range(rep_data$y, na.rm = TRUE)
        p_sc <- ggplot2::ggplot(rep_data, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_point(
            shape = 21,
            size = 3.8,
            stroke = 0.6,
            alpha = 0.9,
            fill = col_rep,
            colour = "#1a1a1a"
          ) +
          ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, colour = "#1a1a1a") +
          ggplot2::facet_wrap(~ panel, ncol = 3) +
          ggplot2::scale_x_continuous(limits = xlim) +
          ggplot2::scale_y_continuous(limits = ylim) +
          ggplot2::labs(title = title_text, x = NULL, y = NULL) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(
            strip.text = ggplot2::element_text(size = 8),
            axis.text = ggplot2::element_text(size = 7),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "none"
          )
        if (isTRUE(show_observed)) {
          obs_df <- data.frame(x = y_num, y = y_num)
          p_sc <- p_sc + ggplot2::geom_point(
            data = obs_df,
            shape = 21, size = 3.8, stroke = 0.6, alpha = 0.95,
            fill = col_obs, colour = "#1a1a1a", inherit.aes = FALSE
          )
        }
        p_sc
      },
      "scatter_avg" = {
        dots <- list(...)
        point_size <- if (!is.null(dots$point_size)) dots$point_size else 3.8
        point_fill <- if (!is.null(dots$point_fill)) dots$point_fill else "#b3cde0"
        point_edge <- if (!is.null(dots$point_edge)) dots$point_edge else "#1a4c78"
        point_alpha <- if (!is.null(dots$point_alpha)) dots$point_alpha else 0.95
        line_colour <- if (!is.null(dots$line_colour)) dots$line_colour else "#000000"
        line_width <- if (!is.null(dots$line_width)) dots$line_width else 0.9
        line_linetype <- if (!is.null(dots$line_linetype)) dots$line_linetype else "dashed"
        ybar <- as.numeric(colMeans(yrep))
        len <- min(length(y), length(ybar))
        df <- data.frame(ybar = ybar[seq_len(len)], y = as.numeric(y[seq_len(len)]))
        ggplot2::ggplot(df, ggplot2::aes(x = ybar, y = y)) +
          ggplot2::geom_point(
            shape = 21,
            size = point_size,
            stroke = 0.6,
            alpha = point_alpha,
            fill = point_fill,
            colour = point_edge
          ) +
          ggplot2::geom_smooth(
            method = "lm",
            se = FALSE,
            linetype = line_linetype,
            linewidth = line_width,
            colour = line_colour
          ) +
          ggplot2::labs(title = title_text, x = "Average yrep", y = "y") +
          ggplot2::theme_minimal(base_size = 11)
      }
    ) # end of switch
    
    if (is_prior) {
      class(p) <- c("qbrms_prior_ggplot", class(p))
      attr(p, "qbrms_meta") <- list(
        type = type, y = y, yrep = yrep, title = title_text, show_observed = show_observed
      )
    }
    return(p)
  }
  
  # base graphics fallback
  cls <- if (is_prior) "qbrms_prior_plot" else "qbrms_plot"
  structure(list(type = type, y = y, yrep = yrep, title = title_text), class = cls)
}

# -----------------------------------------------------------------------------
# Gaussian / generic methods
# -----------------------------------------------------------------------------

#' Posterior predictive checks for qbrms fits
#'
#' @rdname pp_check
#' @param type Character string indicating the check type: one of
#'   \code{"dens_overlay"}, \code{"hist"}, \code{"scatter"}, \code{"scatter_avg"}.
#' @param ndraws Integer number of draws to use.
#' @param seed Optional RNG seed.
#' @param show_observed Logical; show observed data where applicable.
#' @export
pp_check.qbrms_fit <- function(
    object,
    type = "dens_overlay",
    ndraws = 5000,
    seed = NULL,
    show_observed = FALSE,
    ...
) {
  .pp_check_core(
    object, type = type, ndraws = ndraws, seed = seed,
    is_prior = inherits(object, "qbrms_prior_fit"),
    show_observed = show_observed, ...
  )
}

#' Prior predictive checks for qbrms prior objects
#'
#' @rdname pp_check
#' @export
pp_check.qbrms_prior <- function(
    object,
    type = "dens_overlay",
    ndraws = 5000,
    seed = NULL,
    show_observed = FALSE,
    ...
) {
  .pp_check_core(
    object, type = type, ndraws = ndraws, seed = seed,
    is_prior = TRUE, show_observed = show_observed, ...
  )
}

#' @export
pp_check.default <- function(object, ...) {
  stop("No pp_check() method for objects of class: ", paste(class(object), collapse = "/"))
}

# example calibration plot
calibration_plot <- function(plot_data) {
  if (nrow(plot_data) > 0) {
    ggplot2::ggplot(plot_data, ggplot2::aes(x = predicted_prob, y = observed_freq)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_point(ggplot2::aes(size = bin_size, color = category), alpha = 0.7) +
      ggplot2::geom_smooth(
        ggplot2::aes(color = category),
        method = "loess", se = FALSE, linewidth = 1
      ) +
      ggplot2::facet_wrap(~ category) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(
        title = "Calibration Plot by Category",
        x = "Predicted Probability", y = "Observed Frequency"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  } else {
    ggplot2::ggplot() + ggplot2::labs(title = "Calibration plot: insufficient data")
  }
}

# -----------------------------------------------------------------------------
# base-graphics printing (fallback objects)
# -----------------------------------------------------------------------------

#' @export
print.qbrms_plot <- function(x, ...) {
  y <- x$y
  yrep <- x$yrep
  title <- x$title
  if (x$type == "dens_overlay") {
    graphics::hist(y, breaks = 30, probability = TRUE, main = title, xlab = "Value")
    k <- min(5L, nrow(yrep))
    for (i in seq_len(k)) {
      d <- try(stats::density(as.numeric(yrep[i, ])), silent = TRUE)
      if (!inherits(d, "try-error")) graphics::lines(d)
    }
  } else if (x$type == "hist") {
    graphics::hist(y, breaks = 30, main = title, xlab = "Value")
    m <- mean(y, na.rm = TRUE); if (is.finite(m)) graphics::abline(v = m, lwd = 2)
  } else if (x$type == "scatter") {
    graphics::plot(seq_along(y), y, pch = 16, main = title, xlab = "Index", ylab = "Value")
  } else if (x$type == "scatter_avg") {
    ybar <- as.numeric(colMeans(yrep))
    graphics::plot(ybar, y, pch = 16, main = title, xlab = "Average yrep", ylab = "y")
    graphics::abline(stats::lm(y ~ ybar), lty = 2)
  } else {
    graphics::plot.new(); graphics::title(main = title)
  }
  invisible(x)
}

#' @export
print.qbrms_prior_plot <- function(x, ...) {
  print(structure(x, class = "qbrms_plot"))
  cat("[Prior Predictive Plot]\n")
  invisible(x)
}

#' @export
print.qbrms_prior_ggplot <- function(x, ...) {
  NextMethod("print")
  cat("[Prior Predictive Check]\n")
  invisible(x)
}

# -----------------------------------------------------------------------------
# Posterior predictive generator (continuous fallback)
# -----------------------------------------------------------------------------

#' Generate posterior predictive samples
#' @keywords internal
generate_posterior_predictive_samples <- function(object, ndraws = 100) {
  if (!is.null(object$yrep) && is.matrix(object$yrep)) {
    if (nrow(object$yrep) >= ndraws) return(object$yrep[seq_len(ndraws), , drop = FALSE])
    return(object$yrep)
  }
  
  fam <- tryCatch(.extract_family_name_simple(object$family), error = function(e) "unknown")
  if (!is.null(object$data) && fam %in% c("gaussian", "Gaussian", "normal")) {
    X <- tryCatch(.qbrms_model_matrix_fixed(object, object$data),
                  error = function(e) NULL)
    beta <- tryCatch({
      if (!is.null(object$fit$summary.fixed)) {
        as.numeric(object$fit$summary.fixed[, "mean"])
      } else if (!is.null(coef(object))) {
        as.numeric(coef(object))
      } else NULL
    }, error = function(e) NULL)
    beta_names <- tryCatch({
      if (!is.null(object$fit$summary.fixed)) rownames(object$fit$summary.fixed)
      else if (!is.null(coef(object))) names(coef(object)) else NULL
    }, error = function(e) NULL)
    if (!is.null(X) && !is.null(beta) && !is.null(beta_names)) {
      b <- rep(0, ncol(X)); names(b) <- colnames(X)
      match_ok <- intersect(names(b), beta_names)
      b[match_ok] <- beta[match(names(b), beta_names, nomatch = 0L)]
      mu <- as.vector(X %*% b)
      sigma <- tryCatch({
        if (!is.null(object$sigma)) {
          as.numeric(object$sigma)
        } else if (!is.null(object$fit$summary.hyperpar)) {
          hp <- object$fit$summary.hyperpar
          prec_row <- grep("Precision|prec", rownames(hp), ignore.case = TRUE, value = TRUE)
          if (length(prec_row) >= 1) 1 / sqrt(as.numeric(hp[prec_row[1], "mean"])) else NA_real_
        } else NA_real_
      }, error = function(e) NA_real_)
      if (!is.finite(sigma)) {
        y <- object$data[[all.vars(object$original_formula)[1]]]
        y <- suppressWarnings(as.numeric(y))
        sigma <- stats::sd(stats::na.omit(y - mu))
        if (!is.finite(sigma) || sigma <= 0) sigma <- stats::sd(stats::na.omit(y))
        if (!is.finite(sigma) || sigma <= 0) sigma <- 1
      }
      n <- nrow(X)
      out <- matrix(NA_real_, nrow = ndraws, ncol = n)
      for (i in seq_len(ndraws)) out[i, ] <- stats::rnorm(n, mean = mu, sd = sigma)
      return(out)
    }
  }
  n_obs <- if (!is.null(object$data)) nrow(object$data) else 100L
  matrix(stats::rnorm(ndraws * n_obs, mean = 0, sd = 1), nrow = ndraws, ncol = n_obs)
}
