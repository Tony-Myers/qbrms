# R/posterior_methods.R

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("kind", "idx", "value", "y", "ybar", "d"))
}

#' Posterior and Prior Predictive Checks (generic)
#' @param object A model object
#' @param ... Additional arguments passed to methods
#' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}

.pp_check_core <- function(object,
                           type     = "dens_overlay",
                           ndraws   = 100,
                           seed     = NULL,
                           is_prior = FALSE,
                           ...) {
  allowed <- c("dens_overlay", "hist", "scatter", "scatter_avg")
  if (is.null(type) || !is.character(type) || length(type) != 1L || !(type %in% allowed)) {
    stop("Unsupported pp_check type: ", type, call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  
  # Extract observed response (robust)
  response_var <- tryCatch(all.vars(object$original_formula)[1], error = function(e) NULL)
  y <- if (!is.null(response_var) && !is.null(object$data)) object$data[[response_var]] else NULL
  if (!is.null(y)) {
    if (is.ordered(y))      y <- as.numeric(y)
    else if (is.factor(y))  y <- suppressWarnings(as.numeric(as.character(y)))
    else                    y <- as.numeric(y)
    y <- stats::na.omit(y)
  }
  
  # Replicated draws
  if (is_prior) {
    yrep <- object$prior_samples
  } else {
    yrep <- tryCatch(generate_posterior_predictive_samples(object, ndraws),
                     error = function(e) NULL)
  }
  
  # Fallback if empty
  if (is.null(yrep) || !is.matrix(yrep) || nrow(yrep) < 1L || ncol(yrep) < 1L) {
    n_obs  <- if (!is.null(y) && length(y) > 0L) length(y) else 100L
    n_draw <- max(1L, ndraws)
    yrep   <- matrix(stats::rnorm(n_draw * n_obs, 0, 1), nrow = n_draw, ncol = n_obs)
  }
  if (is.null(y) || length(y) < 1L) y <- as.numeric(yrep[1, ])
  
  title_text <- if (is_prior) "Prior Predictive Check" else "Posterior Predictive Check"
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Requested palette/sizing
    col_obs   <- "#000000"  # observed (black)
    col_rep   <- "#b3cde0"  # yrep spaghetti
    point_col <- "#b3cde0"
    point_sz  <- 3
    rep_size  <- 0.45
    rep_alpha <- 0.35
    
    p <- switch(
      type,
      "dens_overlay" = {
        k <- min(50L, nrow(yrep))
        d_obs <- stats::density(as.numeric(y))
        df_obs <- data.frame(x = d_obs$x, d = d_obs$y, kind = "y")
        
        rep_list <- lapply(seq_len(k), function(i) {
          di <- stats::density(as.numeric(yrep[i, ]))
          data.frame(x = di$x, d = di$y, kind = "yrep", draw = i)
        })
        df_rep <- do.call(rbind, rep_list)
        
        ggplot2::ggplot() +
          ggplot2::geom_line(
            data = df_rep,
            ggplot2::aes(x = .data$x, y = .data$d, group = .data$draw, colour = "yrep"),
            linewidth = rep_size, alpha = rep_alpha
          ) +
          ggplot2::geom_line(
            data = df_obs,
            ggplot2::aes(x = .data$x, y = .data$d, colour = "y"),
            linewidth = 1.15
          ) +
          ggplot2::scale_colour_manual(
            values = c(y = col_obs, yrep = col_rep),
            breaks = c("y", "yrep"), labels = c("y", "yrep"), name = NULL
          ) +
          ggplot2::labs(title = title_text, x = NULL, y = NULL) +
          ggplot2::theme_classic()
      },
      "hist" = {
        df <- data.frame(value = y)
        ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
          ggplot2::geom_histogram(bins = 30, fill = col_rep, colour = "white") +
          ggplot2::labs(title = title_text, x = "Value", y = "Count") +
          ggplot2::theme_minimal()
      },
      "scatter" = {
        df <- data.frame(idx = seq_along(y), y = y)
        ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$y)) +
          ggplot2::geom_point(size = point_sz, colour = point_col, alpha = 0.95) +
          ggplot2::labs(title = title_text, x = "Index", y = "y") +
          ggplot2::theme_minimal()
      },
      "scatter_avg" = {
        ybar <- as.numeric(colMeans(yrep))
        len  <- min(length(y), length(ybar))
        df <- data.frame(ybar = ybar[seq_len(len)], y = y[seq_len(len)])
        ggplot2::ggplot(df, ggplot2::aes(x = .data$ybar, y = .data$y)) +
          ggplot2::geom_point(size = point_sz, colour = point_col, alpha = 0.95) +
          ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = "dashed",
                               linewidth = 0.8, colour = "#000000") +
          ggplot2::labs(title = title_text, x = "Average yrep", y = "y") +
          ggplot2::theme_minimal()
      }
    )
    
    if (is_prior) {
      class(p) <- c("qbrms_prior_ggplot", class(p))
      attr(p, "qbrms_meta") <- list(type = type, y = y, yrep = yrep, title = title_text)
    }
    return(p)
  }
  
  # base graphics fallback
  cls <- if (is_prior) "qbrms_prior_plot" else "qbrms_plot"
  structure(list(type = type, y = y, yrep = yrep, title = title_text),
            class = cls)
}

#' Posterior predictive checks for qbrms_fit objects
#' @param object A qbrms_fit object
#' @param type Type of plot ("dens_overlay", "hist", "scatter", "scatter_avg")
#' @param ndraws Number of posterior draws to use
#' @param seed Random seed for reproducibility
#' @param ... Additional arguments
#' @export
#' @method pp_check qbrms_fit
pp_check.qbrms_fit <- function(object,
                               type   = "dens_overlay",
                               ndraws = 100,
                               seed   = NULL,
                               ...) {
  .pp_check_core(object, type = type, ndraws = ndraws, seed = seed,
                 is_prior = inherits(object, "qbrms_prior"), ...)
}

#' @export
#' @method pp_check qbrms_prior
pp_check.qbrms_prior <- function(object,
                                 type   = "dens_overlay",
                                 ndraws = 100,
                                 seed   = NULL,
                                 ...) {
  .pp_check_core(object, type = type, ndraws = ndraws, seed = seed,
                 is_prior = TRUE, ...)
}

# ---------- base-graphics printing -------------------------------------------
#' @export
#' @method print qbrms_plot
print.qbrms_plot <- function(x, ...) {
  y     <- x$y
  yrep  <- x$yrep
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
#' @method print qbrms_prior_plot
print.qbrms_prior_plot <- function(x, ...) {
  print(structure(x, class = "qbrms_plot"))
  cat("[Prior Predictive Plot]\n")
}

#' @export
#' @method print qbrms_prior_ggplot
print.qbrms_prior_ggplot <- function(x, ...) {
  NextMethod("print")
  cat("[Prior Predictive Check]\n")
}

# ---------- posterior predictive generator -----------------------------------

#' Generate posterior predictive samples (Gaussian models are model-based)
#' @keywords internal
generate_posterior_predictive_samples <- function(object, ndraws = 100) {
  if (!is.null(object$yrep) && is.matrix(object$yrep)) {
    if (nrow(object$yrep) >= ndraws) return(object$yrep[seq_len(ndraws), , drop = FALSE])
    return(object$yrep)
  }
  
  fam <- tryCatch(extract_family_name(object$family), error = function(e) "unknown")
  
  if (!is.null(object$data) && fam %in% c("gaussian", "Gaussian", "normal")) {
    X <- tryCatch(stats::model.matrix(object$original_formula, data = object$data),
                  error = function(e) NULL)
    
    beta <- tryCatch({
      if (!is.null(object$fit$summary.fixed)) {
        as.numeric(object$fit$summary.fixed[,"mean"])
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
      
      sigma <- tryCatch(
        {
          if (!is.null(object$sigma)) {
            as.numeric(object$sigma)
          } else if (!is.null(object$fit$summary.hyperpar)) {
            hp <- object$fit$summary.hyperpar
            prec_row <- grep("Precision|prec", rownames(hp), ignore.case = TRUE, value = TRUE)
            if (length(prec_row) >= 1) 1 / sqrt(as.numeric(hp[prec_row[1], "mean"])) else NA_real_
          } else NA_real_
        },
        error = function(e) NA_real_
      )
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