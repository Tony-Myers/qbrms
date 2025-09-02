# =============================================================================
# Conditional effects (qbrms) with spaghetti and ggplot-friendly output
# =============================================================================

#' Conditional effects (generic)
#'
#' @description
#' Generic for computing one-dimensional conditional effects / marginal fitted
#' values as a predictor varies (others held fixed). Methods should return a
#' \pkg{ggplot2} object (or a named list of them) that users can further modify.
#'
#' @param object A model object.
#' @param ... Passed to methods.
#' @export
conditional_effects <- function(object, ...) {
  UseMethod("conditional_effects")
}

#' Conditional effects for qbrms models
#'
#' @description
#' Compute one-dimensional conditional effects for \code{qbrms_fit} objects.
#' Returns a single \code{ggplot} (if one effect) or a named list of
#' \code{ggplot}s (if multiple effects). Users can modify the plots with
#' ggplot syntax (e.g., \code{+ theme_minimal()}).
#'
#' @param object A \code{qbrms_fit} object.
#' @param effects Character vector of predictor names to vary. If \code{NULL}
#'   (default) we try to infer main-effect variables from the model formula.
#' @param prob Credible interval mass (default \code{0.95}).
#' @param spaghetti Logical; if \code{TRUE}, overlay lines from coefficient
#'   draws. Default \code{FALSE}.
#' @param ndraws Number of coefficient draws for spaghetti (default \code{200}).
#' @param grid_points Number of x-points in the grid (default \code{100}).
#' @param seed Optional RNG seed.
#' @param ... Unused.
#'
#' @return A \code{ggplot} or a named list of \code{ggplot}s with class
#'   \code{"qbrms_conditional_effects"}.
#'
#' @export
#' @method conditional_effects qbrms_fit
#' @importFrom stats coef vcov model.matrix qnorm terms
conditional_effects.qbrms_fit <- function(object,
                                          effects = NULL,
                                          prob = 0.95,
                                          spaghetti = FALSE,
                                          ndraws = 200,
                                          grid_points = 100,
                                          seed = NULL,
                                          ...) {
  stopifnot(inherits(object, "qbrms_fit"))
  
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  typical_value <- function(v) {
    if (is.numeric(v)) return(mean(v, na.rm = TRUE))
    if (is.logical(v)) return(FALSE)
    if (is.factor(v) && !is.ordered(v)) return(levels(v)[1])
    if (is.ordered(v)) {
      lev <- levels(v); return(lev[ceiling(length(lev) / 2)])
    }
    if (is.character(v)) return(v[which.max(table(v))][1])
    stats::na.omit(v)[1] %||% NA
  }
  
  build_reference_row <- function(df) {
    as.data.frame(lapply(df, typical_value), stringsAsFactors = FALSE)
  }
  
  safe_vcov <- function(obj) {
    v <- try(stats::vcov(obj), silent = TRUE)
    if (inherits(v, "try-error") || is.null(v) || anyNA(v)) {
      cf <- stats::coef(obj)
      return(diag(rep(1e-6, length(cf))))
    }
    v
  }
  
  draw_betas <- function(beta, V, ndraws) {
    if (!is.matrix(V)) V <- as.matrix(V)
    p <- length(beta)
    if (nrow(V) != p || ncol(V) != p) V <- diag(rep(1e-6, p))
    L <- try(chol(V), silent = TRUE)
    if (inherits(L, "try-error")) {
      eps <- 1e-7
      for (k in 1:6) {
        L <- try(chol(V + diag(eps, p)), silent = TRUE)
        if (!inherits(L, "try-error")) break
        eps <- eps * 10
      }
    }
    if (inherits(L, "try-error")) L <- diag(sqrt(pmax(1e-6, diag(V))), p, p)
    Z <- matrix(stats::rnorm(p * ndraws), nrow = p, ncol = ndraws)
    sweep(L %*% Z, 1, beta, `+`)
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  data <- object$data %||% stop("object$data is required.", call. = FALSE)
  form <- object$original_formula %||% stop("object$original_formula missing.", call. = FALSE)
  
  beta <- stats::coef(object)
  if (is.null(beta)) stop("coef(object) returned NULL.", call. = FALSE)
  V <- safe_vcov(object)
  
  mm_all <- stats::model.matrix(form, data)
  if (is.null(effects)) {
    trm <- stats::terms(form)
    tl  <- attr(trm, "term.labels") %||% character()
    # try to keep only simple main effects
    cand <- unique(unlist(strsplit(tl, "[:*^+]")))
    cand <- cand[cand %in% names(data)]
    if (length(cand) == 0L) cand <- setdiff(names(data), "(Intercept)")
    effects <- cand
  }
  
  ref <- build_reference_row(data)
  out_plots <- list()
  
  for (ef in effects) {
    if (!ef %in% names(data)) next
    
    xseq <- if (is.numeric(data[[ef]])) {
      rng <- range(data[[ef]], na.rm = TRUE)
      seq(rng[1], rng[2], length.out = grid_points)
    } else if (is.factor(data[[ef]])) {
      levels(data[[ef]])
    } else {
      rng <- range(as.numeric(data[[ef]]), na.rm = TRUE)
      seq(rng[1], rng[2], length.out = grid_points)
    }
    
    pred <- ref[rep(1, length(xseq)), , drop = FALSE]
    pred[[ef]] <- xseq
    if (is.factor(data[[ef]])) {
      pred[[ef]] <- factor(pred[[ef]], levels = levels(data[[ef]]),
                           ordered = is.ordered(data[[ef]]))
    }
    
    X <- stats::model.matrix(form, data = pred)
    
    # align beta with X
    if (!is.null(names(beta))) {
      if ("(Intercept)" %in% colnames(X) && !("(Intercept)" %in% names(beta))) {
        if ("Intercept" %in% names(beta)) {
          names(beta)[names(beta) == "Intercept"] <- "(Intercept)"
        }
      }
      beta_use <- beta[colnames(X)]
      nas <- is.na(beta_use)
      if (any(nas)) beta_use[nas] <- 0
    } else {
      beta_use <- beta
    }
    
    mu <- as.numeric(X %*% beta_use)
    
    row_se <- vapply(seq_len(nrow(X)), function(i) {
      xi <- X[i, , drop = FALSE]
      sqrt(drop(xi %*% V %*% t(xi)))
    }, numeric(1))
    
    z <- stats::qnorm((1 + prob) / 2)
    lower <- mu - z * row_se
    upper <- mu + z * row_se
    
    df_sp <- NULL
    if (isTRUE(spaghetti) && ndraws > 0) {
      B <- draw_betas(beta_use, V, ndraws)       # p x ndraws
      eta <- X %*% B                             # n x ndraws
      df_sp <- data.frame(
        effect1__  = rep(pred[[ef]], ndraws),
        estimate__ = as.numeric(eta),
        id__       = rep(seq_len(ndraws), each = nrow(pred)),
        stringsAsFactors = FALSE
      )
    }
    
    df_main <- data.frame(
      effect1__  = pred[[ef]],
      estimate__ = mu,
      lower__    = lower,
      upper__    = upper,
      stringsAsFactors = FALSE
    )
    
    # FIXED: Only show ribbon when spaghetti=FALSE
    p <- ggplot2::ggplot(df_main, ggplot2::aes(x = .data$effect1__, y = .data$estimate__)) +
      { if (!isTRUE(spaghetti))
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower__, ymax = .data$upper__),
                             fill = "grey70", alpha = 0.5)
        else
          ggplot2::geom_blank()
      } +
      { if (!is.null(df_sp))
        ggplot2::geom_line(
          data = df_sp,
          ggplot2::aes(x = .data$effect1__, y = .data$estimate__, group = .data$id__),
          linewidth = 0.3, alpha = 0.25, colour = "#6497b1"
        )
        else
          ggplot2::geom_blank()
      } +
      { if (!isTRUE(spaghetti))
        ggplot2::geom_line(linewidth = 0.6, colour = "blue")
        else
          ggplot2::geom_blank()
      } +
      ggplot2::labs(
        title = paste("Conditional effect:", ef),
        subtitle = paste0(round(prob * 100), "% credible interval"),
        x = ef, y = "Expected value"
      ) +
      ggplot2::theme_minimal()
    
    out_plots[[ef]] <- p
  }
  
  class(out_plots) <- c("qbrms_conditional_effects", "list")
  if (length(out_plots) == 1L) return(out_plots[[1]])
  out_plots
}

#' @export
#' @method print qbrms_conditional_effects
print.qbrms_conditional_effects <- function(x, ...) {
  if (inherits(x, "list")) {
    print(x[[1]])
    invisible(x)
  } else {
    NextMethod()
  }
}

#' @export
#' @method plot qbrms_conditional_effects
plot.qbrms_conditional_effects <- function(x, ...) {
  print.qbrms_conditional_effects(x, ...)
}