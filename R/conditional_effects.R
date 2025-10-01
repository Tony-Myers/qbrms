# =============================================================================
# R/conditional_effects.R - Conditional effects (Gaussian) for qbrms
# =============================================================================

#' Conditional effects (generic)
#'
#' Compute one-dimensional conditional effects / marginal fitted values as a
#' predictor varies while other covariates are held fixed (typically at means /
#' modes). Methods should return an object that `plot()` can visualise.
#'
#' @importFrom stats contrasts<-
#' @param object A model object.
#' @param ... Passed to methods.
#' @export
conditional_effects <- function(object, ...) {
  UseMethod("conditional_effects")
}

# ---- Internal utilities ------------------------------------------------------

.qbrms__ref_row <- function(dat) {
  if (nrow(dat) == 0L) stop("No data found in model object.")
  take_mode <- function(z) {
    tz <- table(z)
    names(tz)[which.max(tz)]
  }
  ref <- lapply(dat[1, , drop = FALSE], function(z) {
    if (is.numeric(z)) return(mean(z, na.rm = TRUE))
    if (is.logical(z)) return(as.logical(take_mode(z)))
    if (is.factor(z))  return(factor(take_mode(z), levels = levels(z), ordered = is.ordered(z)))
    if (inherits(z, "Date"))   return(as.Date(mean(as.numeric(z), na.rm = TRUE), origin = "1970-01-01"))
    if (inherits(z, "POSIXt")) return(as.POSIXct(mean(as.numeric(z), na.rm = TRUE), origin = "1970-01-01", tz = attr(z, "tzone")))
    z[which.max(!is.na(z))]
  })
  as.data.frame(ref, stringsAsFactors = FALSE, check.names = FALSE)
}

.qbrms__apply_at <- function(ref, at) {
  if (!length(at)) return(ref)
  for (nm in names(at)) if (nm %in% names(ref)) ref[[nm]] <- at[[nm]]
  ref
}

.qbrms__build_grid <- function(dat, xvar, n_points = 100L) {
  x <- dat[[xvar]]
  if (!is.numeric(x)) stop("`effects` must name a numeric predictor in the model data.")
  rng <- range(x, na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n_points)
}

.qbrms__cov_source_from <- function(object, V_used) {
  src <- attr(V_used, "cov_source")
  if (is.null(src)) src <- attr(V_used, "source")
  if (is.null(src)) src <- "unknown"
  src
}

.qbrms__stamp_cov_attrs <- function(target, object, V_used) {
  src <- .qbrms__cov_source_from(object, V_used)
  attr(target, "cov_source") <- src
  offdiag <- V_used[upper.tri(V_used, diag = FALSE)]
  attr(target, "used_full_cov") <- isTRUE(any(abs(offdiag) > .Machine$double.eps^0.5))
  target
}

.qbrms__package_blue <- function() "#3F7AE0"

# --- recursively replace any formula that still contains a random-effect bar ---
.qbrms__sanitize_formulas_in_object <- function(x, f_fixed) {
  # Only replace formulas that actually contain "(...|...)" to avoid breaking anything else
  has_bar <- function(frm) {
    if (!inherits(frm, "formula")) return(FALSE)
    grepl("\\([^()]*\\|[^()]*\\)", paste(deparse(frm, width.cutoff = 500L), collapse = ""), perl = TRUE)
  }
  rec <- function(obj) {
    if (is.list(obj)) {
      for (nm in names(obj)) {
        val <- obj[[nm]]
        if (inherits(val, "formula")) {
          if (has_bar(val)) obj[[nm]] <- f_fixed
        } else if (is.list(val)) {
          obj[[nm]] <- rec(val)
        } else if (is.environment(val)) {
          # environments can hold calls/formulas; be conservative
          # don't modify environments to avoid side effects
        }
      }
    }
    obj
  }
  rec(x)
}

# --- helpers ---------------------------------------------------------------

.qb_nobars <- function(form) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required for conditional_effects.")
  lme4::nobars(form)
}

.qb_ref_values <- function(dat) {
  stopifnot(is.data.frame(dat))
  take_mode <- function(z) {
    tz <- table(z)
    names(tz)[which.max(tz)]
  }
  out <- lapply(dat, function(z) {
    if (is.numeric(z)) {
      mean(z, na.rm = TRUE)
    } else if (is.factor(z)) {
      factor(take_mode(z), levels = levels(z))
    } else {
      z[1]
    }
  })
  as.data.frame(out, stringsAsFactors = FALSE)
}

.qb_invlink <- function(fam, eta) {
  fam <- tolower(fam)
  if (grepl("poisson", fam))      return(exp(eta))
  if (grepl("gaussian|normal", fam)) return(eta)
  if (grepl("binomial|bernoulli|logit", fam)) return(1 / (1 + exp(-eta)))
  stop("Family not yet supported in conditional_effects().")
}

# ---- helper: inverse-link from the object -----------------------------------
.qbrms__invlink <- function(object, eta) {
  # Try to discover a link; default to identity if unknown
  link <- NULL
  fam  <- NULL
  if (!is.null(object$family)) fam <- object$family
  if (is.null(fam) && !is.null(object$fit) && !is.null(object$fit$family))
    fam <- object$fit$family
  
  if (is.character(fam)) {
    f <- tolower(fam)
    if (grepl("poisson", f)) link <- "log"
    else if (grepl("binom|bernoulli", f)) link <- "logit"
    else if (grepl("gauss|normal", f)) link <- "identity"
  } else if (is.list(fam) && !is.null(fam$link)) {
    link <- fam$link
  }
  
  if (is.null(link) && is.list(fam) && !is.null(fam$linkfun)) {
    # Some objects carry linkfun/linkinv like stats::family
    return(fam$linkinv(eta))
  }
  
  switch(tolower(link %||% "identity"),
         "log"     = exp(eta),
         "logit"   = 1 / (1 + exp(-eta)),
         "probit"  = stats::pnorm(eta),
         "cloglog" = 1 - exp(-exp(eta)),
         "identity"= eta,
         exp(eta))  # safe default
}

#' Discrete-slice conditional effects (brms-style) for qbrms
#' Build point/interval summaries at a few values of a numeric moderator,
#' plotted against the factor on the x-axis.
#' @param object A qbrms_fit object.
#' @param effects Character vector specifying effects to plot. If NULL, all numeric predictors are used.
#' @param slices Named list of variables and values at which to slice the data.
#' @param nslices Number of slices to use for each slicing variable.
#' @param prob Probability mass to include in uncertainty intervals (default 0.95).
#' @param ndraws Number of posterior draws to use for predictions.
#' @param at Named list of values at which to fix other predictors.
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to prediction functions.
#' @export
conditional_effects_slices <- function(
    object,
    effects,                     # e.g. "Trt:zBase" (factor:numeric OR numeric:factor)
    slices    = NULL,            # numeric values for the moderator; NULL = auto
    nslices   = 3L,              # only used if slices is NULL
    prob      = 0.95,
    ndraws    = 200L,
    at        = list(),
    seed      = NULL,
    ...
) {
  dat <- object$data
  if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0L)
    stop("Model data not found in `object$data`.")
  
  parts <- strsplit(effects, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2 || !all(parts %in% names(dat)))
    stop("`effects` must be like 'factor:numeric' (or 'numeric:factor') present in the model data.")
  
  is_num <- vapply(dat[parts], is.numeric, TRUE)
  if (sum(is_num) != 1L)
    stop("Interaction must contain exactly one numeric and one factor variable.")
  
  numvar <- parts[is_num]
  facvar <- parts[!is_num]
  
  if (!is.factor(dat[[facvar]]))
    dat[[facvar]] <- factor(dat[[facvar]])
  levs <- levels(dat[[facvar]])
  
  if (is.null(slices)) {
    m  <- mean(dat[[numvar]], na.rm = TRUE)
    sd <- stats::sd(dat[[numvar]], na.rm = TRUE)
    if (nslices == 3L && is.finite(m) && is.finite(sd) && abs(m) < 0.15 && abs(sd - 1) < 0.15) {
      slices <- c(-1, 0, 1)             # standardised default
    } else {
      p <- if (nslices == 3L) c(0.1, 0.5, 0.9) else stats::ppoints(nslices)
      slices <- as.numeric(stats::quantile(dat[[numvar]], probs = p, na.rm = TRUE))
    }
  }
  slice_lab <- paste0(numvar, "=", format(slices, digits = 3))
  
  if (!is.null(seed)) set.seed(as.integer(seed))
  
  ref  <- .qbrms__apply_at(.qbrms__ref_row(dat), at)
  
  newd <- do.call(rbind, lapply(seq_along(slices), function(i) {
    v <- slices[i]
    tmp <- ref[rep(1, length(levs)), , drop = FALSE]
    tmp[[facvar]] <- factor(levs, levels = levs)
    tmp[[numvar]] <- v
    tmp[[".slice"]] <- factor(rep(slice_lab[i], length(levs)), levels = slice_lab)
    tmp
  }))
  rownames(newd) <- NULL
  
  # fixed-effects design with training contrasts
  # fixed-effects design with training contrasts (bar-proof)
  f_raw   <- object$original_formula
  f_fixed <- .qbrms__remove_random_effects(f_raw)
  
  # Belt-and-braces: if anything with '|' survived, strip "(...|...)" by regex
  f_txt <- paste(deparse(f_fixed), collapse = "")
  if (grepl("\\|", f_txt, perl = TRUE)) {
    f_fixed <- stats::as.formula(
      gsub("\\([^()]*\\|[^()]*\\)", "", paste(deparse(f_raw), collapse = ""))
    )
  }
  environment(f_fixed) <- environment(f_raw)
  
  tt     <- stats::delete.response(stats::terms(f_fixed, data = dat))
  Xtrain <- stats::model.matrix(tt, dat)
  ctr    <- attr(Xtrain, "contrasts")
  
  for (nm in names(dat)) if (is.factor(dat[[nm]])) {
    newd[[nm]] <- factor(newd[[nm]], levels = levels(dat[[nm]]))
    if (!is.null(ctr) && !is.null(ctr[[nm]])) contrasts(newd[[nm]]) <- ctr[[nm]]
  }
  Xgrid <- stats::model.matrix(tt, newd, contrasts.arg = ctr)
  
  # ---- Align to fitted coefficient order
  summ   <- object$fit$summary.fixed
  if (is.null(summ)) stop("No summary.fixed available in fit.")
  bnames <- rownames(summ)
  bmean  <- setNames(as.numeric(summ[, "mean"]), bnames)
  
  if (!all(bnames %in% colnames(Xgrid))) {
    stop("Design matrix is missing columns: ", paste(setdiff(bnames, colnames(Xgrid)), collapse = ", "))
  }
  Xgrid <- Xgrid[, bnames, drop = FALSE]
  
  # ---- Use a sanitized copy of the object for vcov retrieval (prevents '|' warnings)
  obj2 <- object
  obj2$original_formula <- f_fixed
  if (!is.null(obj2$fit) && (is.list(obj2$fit) || is.environment(obj2$fit))) {
    if (!is.null(obj2$fit$formula) && inherits(obj2$fit$formula, "formula")) {
      obj2$fit$formula <- f_fixed
    }
    if (!is.null(obj2$fit$call)) {
      cl <- as.list(obj2$fit$call)
      if ("formula" %in% names(cl)) obj2$fit$call[["formula"]] <- f_fixed
    }
  }
  obj2 <- .qbrms__sanitize_formulas_in_object(obj2, f_fixed)
  
  V <- try(.qbrms_get_vcov_fixed(obj2), silent = TRUE)
  if (inherits(V, "try-error") || is.null(V)) {
    sds <- as.numeric(summ[, "sd"])
    V   <- diag(sds^2, nrow = length(bnames))
    rownames(V) <- colnames(V) <- bnames
  } else {
    V <- V[bnames, bnames, drop = FALSE]
  }
  
  # ---- coefficient draws
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required for conditional_effects_slices().")
  beta_draws <- MASS::mvrnorm(n = ndraws, mu = bmean, Sigma = V)
  
  
  eta_draws <- Xgrid %*% t(beta_draws)
  mu_draws  <- .qbrms__invlink(object, eta_draws)
  
  alpha <- (1 - prob) / 2
  
  df <- data.frame(
    x          = newd[[facvar]],
    estimate__ = rowMeans(mu_draws),
    lower__    = apply(mu_draws, 1, stats::quantile, probs = alpha),
    upper__    = apply(mu_draws, 1, stats::quantile, probs = 1 - alpha),
    slice      = newd[[".slice"]],
    check.names = FALSE
  )
  names(df)[names(df) == "x"] <- facvar
  attr(df, "kind") <- "summary"
  attr(df, "slice_name") <- numvar
  df <- .qbrms__stamp_cov_attrs(df, object, V)
  
  out <- list()
  out[[effects]] <- df
  class(out) <- c("qbrms_conditional_effects", "list")
  out
}


# ---- Method for qbrms fits (Gaussian) ---------------------------------------

#' Conditional effects for qbrms Gaussian models
#'
#' @param object A qbrms fit object (Gaussian).
#' @param effects Character vector: names of predictors to vary. Supports simple
#'   two-way interactions "num:fac" or "fac:num" where one is numeric and the other factor.
#' @param spaghetti Logical; if TRUE draw per-draw "spaghetti" lines. If FALSE,
#'   draw a mean line with a credible-interval ribbon.
#' @param ndraws Number of joint coefficient draws for uncertainty (default 200).
#' @param n_points Size of the x-grid across the observed range (default 100).
#' @param at Optional named list of covariate values to hold constant.
#' @param seed Optional integer seed for reproducibility.
#' @param prob Interval probability for ribbons (default 0.95).
#' @param ... Ignored.
#' @export
#' @export
conditional_effects.qbrms_fit <- function(
    object,
    effects  = NULL,
    spaghetti = FALSE,
    ndraws   = 200L,
    n_points = 100L,
    at       = list(),
    seed     = NULL,
    prob     = 0.95,
    ...
) {
  dat <- object$data
  if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0L) {
    if (!is.null(object$fit) && is.data.frame(object$fit$frame) && nrow(object$fit$frame) > 0L) {
      dat <- object$fit$frame
    } else {
      stop("Model data not found in `object$data` (and `object$fit$frame` was empty). Please keep training data in the fit object.")
    }
  }
  
  # Accept `effects = NULL` (brms-like defaults)
  if (is.null(effects)) {
    f_fixed <- .qbrms__remove_random_effects(object$original_formula)
    tt <- stats::terms(f_fixed, data = dat)
    tl <- attr(tt, "term.labels")
    
    # main numeric effects
    main <- tl[!grepl(":", tl, fixed = TRUE)]
    main <- main[main %in% names(dat)]
    main <- main[vapply(dat[main], is.numeric, TRUE)]
    
    # numeric:factor interactions
    inter <- tl[grepl(":", tl, fixed = TRUE)]
    inter <- inter[vapply(inter, function(lbl) {
      parts <- strsplit(lbl, ":", fixed = TRUE)[[1]]
      if (!all(parts %in% names(dat))) return(FALSE)
      is_num <- vapply(dat[parts], is.numeric, TRUE)
      sum(is_num) == 1L  # exactly one numeric, one non-numeric (factor)
    }, logical(1))]
    
    effects <- c(main, inter)
    
    if (!length(effects)) {
      stop("No suitable numeric predictors found to compute conditional effects; pass `effects=` explicitly.")
    }
  } else if (!is.character(effects) || !length(effects)) {
    stop("`effects` must be a non-empty character vector (or NULL for defaults).")
  } else {
    effects <- effects[nzchar(effects)]
  }
  
  if (!is.null(seed)) set.seed(as.integer(seed))
  ce_list <- list()
  
  for (current_effect in effects) {
    
    
    # Identify focal numeric x and optional by-factor from "A:B"
    byvar  <- NULL
    x_name <- current_effect
    if (grepl(":", current_effect, fixed = TRUE)) {
      parts <- strsplit(current_effect, ":", fixed = TRUE)[[1]]
      if (!all(parts %in% names(dat))) stop("Interaction parts not found in the model data.")
      is_num <- vapply(dat[parts], is.numeric, TRUE)
      if (sum(is_num) != 1L) stop("Interaction must be numeric:factor or factor:numeric.")
      x_name <- parts[is_num]
      byvar  <- parts[!is_num]
    }
    
    # Reference row and grid
    ref   <- .qbrms__apply_at(.qbrms__ref_row(dat), at)
    xgrid <- .qbrms__build_grid(dat, x_name, n_points)
    
    newd  <- ref[rep(1, length(xgrid)), , drop = FALSE]
    newd[[x_name]] <- xgrid
    
    if (!is.null(byvar)) {
      levs <- levels(as.factor(dat[[byvar]]))
      newd <- do.call(rbind, lapply(levs, function(L) {
        tmp <- newd
        tmp[[byvar]] <- factor(L, levels = levs)
        tmp
      }))
      rownames(newd) <- NULL
    }
    
    # ----- Build fixed-effects design FIRST (order matters) -----------------
    f_raw   <- object$original_formula
    f_fixed <- .qbrms__remove_random_effects(f_raw)
    
    # Extra guard against any lingering "(...|...)"
    f_txt <- paste(deparse(f_fixed), collapse = "")
    if (grepl("\\|", f_txt, perl = TRUE)) {
      f_fixed <- stats::as.formula(
        gsub("\\([^()]*\\|[^()]*\\)", "", paste(deparse(f_raw), collapse = ""))
      )
    }
    environment(f_fixed) <- environment(f_raw)
    
    tt       <- stats::delete.response(stats::terms(f_fixed, data = dat))
    Xtrain   <- stats::model.matrix(tt, dat)
    ctr_list <- attr(Xtrain, "contrasts")
    
    # Ensure factor levels AND contrasts in newd match the fit exactly
    for (nm in names(dat)) {
      if (is.factor(dat[[nm]])) {
        newd[[nm]] <- factor(newd[[nm]], levels = levels(dat[[nm]]))
        if (!is.null(ctr_list) && !is.null(ctr_list[[nm]])) {
          contrasts(newd[[nm]]) <- ctr_list[[nm]]
        }
      }
    }
    
    # Grid design using identical terms and contrasts as the fit
    Xgrid <- stats::model.matrix(tt, newd, contrasts.arg = ctr_list)
    
    # Sanitised copy for any helper that inspects $fit$formula/$fit$call
    obj2 <- object
    obj2$original_formula <- f_fixed
    if (!is.null(obj2$fit) && (is.list(obj2$fit) || is.environment(obj2$fit))) {
      if (!is.null(obj2$fit$formula) && inherits(obj2$fit$formula, "formula")) {
        obj2$fit$formula <- f_fixed
      }
      if (!is.null(obj2$fit$call) && (is.call(obj2$fit$call) || is.language(obj2$fit$call))) {
        cl <- as.list(obj2$fit$call)
        if ("formula" %in% names(cl)) obj2$fit$call[["formula"]] <- f_fixed
      }
    }
    obj2 <- .qbrms__sanitize_formulas_in_object(obj2, f_fixed)
    
    # Align to fitted fixed effects
    summ   <- object$fit$summary.fixed
    if (is.null(summ)) stop("No summary.fixed available in fit.")
    bnames <- rownames(summ)
    bmean  <- as.numeric(summ[, "mean"]); names(bmean) <- bnames
    
    if (!all(bnames %in% colnames(Xgrid))) {
      missing_cols <- setdiff(bnames, colnames(Xgrid))
      stop("Design matrix is missing columns: ", paste(missing_cols, collapse = ", "))
    }
    Xgrid <- Xgrid[, bnames, drop = FALSE]
    
    # Covariance for draws: prefer stored vcov on the sanitised object; else diagonal
    V <- try(.qbrms_get_vcov_fixed(obj2), silent = TRUE)
    if (inherits(V, "try-error") || is.null(V)) {
      sds <- as.numeric(summ[, "sd"])
      V   <- diag(sds^2, nrow = length(bnames), ncol = length(bnames))
      rownames(V) <- colnames(V) <- bnames
    } else {
      V <- V[bnames, bnames, drop = FALSE]
    }
    if (is.null(attr(V, "cov_source"))) attr(V, "cov_source") <- "fixed_vcov_from_fit"
    if (is.null(attr(V, "source")))     attr(V, "source")     <- "fixed_vcov_from_fit"
    
    # Joint draws of coefficients
    if (!requireNamespace("MASS", quietly = TRUE))
      stop("Package 'MASS' is required for conditional_effects.")
    beta_draws <- MASS::mvrnorm(n = ndraws, mu = bmean, Sigma = V)
    
    # Linear predictor draws
    eta_draws <- Xgrid %*% t(beta_draws)
    
    if (isTRUE(spaghetti)) {
      # Transform per-draw to response scale
      mu_draws <- .qbrms__invlink(object, eta_draws)
      
      nd <- ncol(mu_draws); nX <- nrow(mu_draws)
      df <- data.frame(
        x      = rep(newd[[x_name]], times = nd),
        .value = as.numeric(mu_draws),
        .draw  = rep(seq_len(nd), each = nX),
        check.names = FALSE
      )
      names(df)[names(df) == "x"] <- x_name
      if (!is.null(byvar)) {
        df[[byvar]] <- rep(newd[[byvar]], times = nd)
        df$..grp..  <- interaction(df$.draw, df[[byvar]], drop = TRUE)
      }
      
      attr(df, "kind") <- "spaghetti"
      df <- .qbrms__stamp_cov_attrs(df, object, V)
      ce_list[[current_effect]] <- df
      
    } else {
      # Summarise on the RESPONSE scale (matches brms behaviour)
      mu_draws <- .qbrms__invlink(object, eta_draws)
      alpha    <- (1 - prob) / 2
      
      df <- data.frame(
        x          = newd[[x_name]],
        estimate__ = rowMeans(mu_draws),
        lower__    = apply(mu_draws, 1, stats::quantile, probs = alpha),
        upper__    = apply(mu_draws, 1, stats::quantile, probs = 1 - alpha),
        check.names = FALSE
      )
      
      names(df)[names(df) == "x"] <- x_name
      if (!is.null(byvar)) df[[byvar]] <- newd[[byvar]]
      
      attr(df, "kind") <- "summary"
      df <- .qbrms__stamp_cov_attrs(df, object, V)
      ce_list[[current_effect]] <- df
    }
  }
  
  class(ce_list) <- c("qbrms_conditional_effects", "list")
  ce_list
}


# ---- Helper functions------------------------

#' Generate spaghetti draws with proper correlation structure
#' @keywords internal
.qbrms__generate_spaghetti_draws_corrected <- function(object, Xgrid, V_used, ndraws) {
  
  summ <- object$fit$summary.fixed
  if (is.null(summ)) stop("No summary.fixed available in fit.")
  
  bmean  <- as.numeric(summ[, "mean"])
  bnames <- rownames(summ)
  names(bmean) <- bnames
  
  V_enhanced <- .qbrms__enhance_covariance_matrix(object, V_used, verbose = FALSE)
  
  Xgrid_aligned <- Xgrid[, bnames, drop = FALSE]
  
  beta_draws <- .qbrms_rmvnorm(ndraws, bmean, V_enhanced)
  colnames(beta_draws) <- bnames
  
  mu_draws <- matrix(NA_real_, nrow = nrow(Xgrid_aligned), ncol = ndraws)
  for (i in 1:ndraws) {
    mu_draws[, i] <- as.numeric(Xgrid_aligned %*% beta_draws[i, ])
  }
  
  mu_draws
}

.qbrms__get_attr_or <- function(x, nm, default) {
  val <- attr(x, nm)
  if (is.null(val)) default else val
}

#' Enhance covariance matrix to ensure proper correlations
#' @keywords internal
.qbrms__enhance_covariance_matrix <- function(object, V_current, verbose = TRUE) {
  
  if (nrow(V_current) <= 1) return(V_current)
  
  cors <- stats::cov2cor(V_current)
  offdiag_cors <- cors[upper.tri(cors, diag = FALSE)]
  max_abs_cor  <- max(abs(offdiag_cors), na.rm = TRUE)
  
  if (isTRUE(verbose)) {
    cat("Current covariance source:", .qbrms__get_attr_or(V_current, "source", "unknown"), "\n")
    cat("Maximum absolute correlation:", round(max_abs_cor, 3), "\n")
  }
  
  source_type      <- .qbrms__get_attr_or(V_current, "source", "unknown")
  is_diagonal_only <- identical(source_type, "diag_sd_only")
  
  if (max_abs_cor < 0.1 || is_diagonal_only) {
    if (isTRUE(verbose)) cat("Weak/diagonal correlations detected, computing geometric correlations...\n")
    
    V_geometric <- .qbrms__compute_geometric_covariance(object, V_current)
    if (!is.null(V_geometric)) {
      cors_geom <- stats::cov2cor(V_geometric)
      offdiag_cors_geom <- cors_geom[upper.tri(cors_geom, diag = FALSE)]
      max_abs_cor_geom  <- max(abs(offdiag_cors_geom), na.rm = TRUE)
      
      if (max_abs_cor_geom > max_abs_cor) {
        if (isTRUE(verbose)) {
          cat("Using geometric covariance (correlations improved to", round(max_abs_cor_geom, 3), ")\n")
          if (nrow(cors_geom) == 2) cat("Intercept-slope correlation:", round(cors_geom[1, 2], 3), "\n")
        }
        attr(V_geometric, "source") <- paste0(source_type, "_enhanced_geometric")
        return(V_geometric)
      }
    }
    
    V_ols <- .qbrms__compute_ols_covariance(object)
    if (!is.null(V_ols)) {
      cors_ols <- stats::cov2cor(V_ols)
      offdiag_cors_ols <- cors_ols[upper.tri(cors_ols, diag = FALSE)]
      max_abs_cor_ols  <- max(abs(offdiag_cors_ols), na.rm = TRUE)
      
      if (max_abs_cor_ols > max_abs_cor) {
        if (isTRUE(verbose)) cat("Using OLS covariance (better correlations)\n")
        attr(V_ols, "source") <- paste0(source_type, "_enhanced_with_ols")
        return(V_ols)
      }
    }
  }
  
  V_current
}

#' Compute geometric covariance from design matrix properties
#' @keywords internal
#' Compute geometric covariance from design matrix properties
#' @keywords internal
.qbrms__compute_geometric_covariance <- function(object, V_current) {
  tryCatch({
    # Guardrails
    if (is.null(object) || is.null(object$data)) return(NULL)
    if (is.null(object$fit) || is.null(object$fit$summary.fixed)) return(NULL)
    
    # 1) Build the training design with the SAME terms and contrasts as the fit
    f       <- object$original_formula
    df      <- object$data
    f_fixed <- .qbrms__remove_random_effects(f)
    tt      <- stats::delete.response(stats::terms(f_fixed, data = df))
    Xtrain  <- stats::model.matrix(tt, df)              # design used at fit time
    contrs  <- attr(Xtrain, "contrasts")
    
    # Use the identical design (we could also rebuild with contrasts.arg = contrs)
    X <- Xtrain
    
    # 2) Align to the fitted fixed effects
    summ   <- object$fit$summary.fixed
    bnames <- rownames(summ)
    if (is.null(bnames)) return(NULL)
    
    common <- intersect(bnames, colnames(X))
    if (length(common) < 2L) return(NULL)              # not enough parameters to form correlations
    X_aligned <- X[, common, drop = FALSE]
    
    # 3) Compute (X'X)^{-1} robustly
    XtX <- crossprod(X_aligned)
    V_geom <- tryCatch(
      solve(XtX),
      error = function(e) {
        if (requireNamespace("MASS", quietly = TRUE)) {
          return(MASS::ginv(XtX))
        } else {
          return(NULL)
        }
      }
    )
    if (is.null(V_geom)) return(NULL)
    
    # 4) Rescale to match the current marginal variances on the common set
    #    (assumes V_current is keyed by coefficient names)
    cur_diag <- diag(V_current)
    names(cur_diag) <- rownames(V_current)
    current_vars <- cur_diag[common]
    geom_vars    <- diag(V_geom)
    # numerical guard against tiny or negative due to rounding
    current_vars[!is.finite(current_vars)] <- 0
    geom_vars[geom_vars <= .Machine$double.eps] <- .Machine$double.eps
    
    scale_factors <- sqrt(current_vars / geom_vars)
    S <- diag(scale_factors, nrow = length(scale_factors), ncol = length(scale_factors))
    
    V_scaled <- S %*% V_geom %*% S
    rownames(V_scaled) <- colnames(V_scaled) <- common
    
    # 5) If the fitted model has more coefficients than 'common', expand back
    if (length(bnames) > length(common)) {
      V_full <- matrix(0, nrow = length(bnames), ncol = length(bnames))
      rownames(V_full) <- colnames(V_full) <- bnames
      # keep existing variances on the diagonal from V_current if available
      if (!is.null(rownames(V_current)) && all(bnames %in% rownames(V_current))) {
        diag(V_full) <- diag(V_current)[bnames]
      } else {
        diag(V_full)[match(common, bnames)] <- diag(V_scaled)
      }
      V_full[common, common] <- V_scaled
      return(V_full)
    }
    
    V_scaled
  }, error = function(e) NULL)
}

#' Compute OLS covariance matrix as fallback
#' @keywords internal
.qbrms__compute_ols_covariance <- function(object) {
  tryCatch({
    f  <- object$original_formula
    df <- object$data
    
    f_fixed <- .qbrms__remove_random_effects(f)
    
    resp <- all.vars(f_fixed)[1]
    df[[resp]] <- suppressWarnings(as.numeric(df[[resp]]))
    
    lm_fit <- stats::lm(f_fixed, data = df)
    V_ols  <- stats::vcov(lm_fit)
    
    summ <- object$fit$summary.fixed
    if (!is.null(summ)) {
      bnames <- rownames(summ)
      common <- intersect(bnames, rownames(V_ols))
      if (length(common) > 0) {
        V_ols <- V_ols[common, common, drop = FALSE]
        V_full <- matrix(0, nrow = length(bnames), ncol = length(bnames))
        rownames(V_full) <- colnames(V_full) <- bnames
        V_full[common, common] <- V_ols
        return(V_full)
      }
    }
    
    V_ols
  }, error = function(e) NULL)
}

#' Remove random effects from formula
#' @keywords internal
.qbrms__remove_random_effects <- function(formula) {
  # Prefer lme4::nobars if available (best-in-class parser)
  if (requireNamespace("lme4", quietly = TRUE)) {
    return(lme4::nobars(formula))
  }
  # Fallback: remove ONLY terms of the form ( ... | ... ) 
  f_str <- paste(deparse(formula), collapse = "")
  # remove all "(...|...)" chunks
  f_str <- gsub("\\([^()]*\\|[^()]*\\)", "", f_str)
  # collapse multiple '+' and clean up dangling '+' around '~'
  f_str <- gsub("\\+\\s*\\+", "+", f_str)
  f_str <- gsub("~\\s*\\+", "~", f_str)
  f_str <- gsub("\\+\\s*$", "", f_str)
  # if RHS becomes empty, use ~ 1
  if (grepl("~\\s*$", f_str)) f_str <- sub("~\\s*$", "~ 1", f_str)
  stats::as.formula(f_str)
}


#' Compute leverage-aware confidence intervals
#' @keywords internal
.qbrms__compute_leverage_uncertainty <- function(Xgrid, V, object, prob = 0.95) {
  
  summ <- object$fit$summary.fixed
  if (is.null(summ)) stop("No summary.fixed available in model object")
  
  bmean <- as.numeric(summ[, "mean"])
  names(bmean) <- rownames(summ)
  
  common <- intersect(names(bmean), colnames(Xgrid))
  if (length(common) == 0) stop("No matching coefficient names between model and design matrix")
  
  X_aligned  <- Xgrid[, common, drop = FALSE]
  V_aligned  <- V[common, common, drop = FALSE]
  b_aligned  <- bmean[common]
  
  mu_mean <- as.numeric(X_aligned %*% b_aligned)
  
  se_pred <- numeric(nrow(X_aligned))
  for (i in 1:nrow(X_aligned)) {
    xi <- X_aligned[i, , drop = FALSE]
    se_pred[i] <- sqrt(as.numeric(xi %*% V_aligned %*% t(xi)))
  }
  
  alpha  <- (1 - prob) / 2
  z_crit <- stats::qnorm(1 - alpha)
  
  lower <- mu_mean - z_crit * se_pred
  upper <- mu_mean + z_crit * se_pred
  
  list(mean = mu_mean, lower = lower, upper = upper, se = se_pred)
}

# ---- Plot method ------------------------------------------------------------
#' @export
plot.qbrms_conditional_effects <- function(x, ...) {
  if (!length(x)) stop("Empty conditional effects object.")
  
  n_effects    <- length(x)
  effect_names <- names(x)
  
  # pick the x variable robustly
  pick_xvar <- function(df, eff) {
    if (eff %in% names(df)) return(eff)
    cand <- setdiff(names(df), c(".value", ".draw", "estimate__", "lower__", "upper__", "..grp.."))
    num  <- cand[vapply(df[cand], is.numeric, TRUE)]
    if (length(num)) return(num[1])
    if (length(cand)) return(cand[1])
    stop("Could not determine x variable for plotting.")
  }
  
  if (n_effects == 1) {
    eff <- effect_names[1]
    df  <- x[[eff]]
    
    # detect kind
    kind <- attr(df, "kind", exact = TRUE)
    if (is.null(kind) && all(c(".draw", ".value") %in% names(df))) kind <- "spaghetti"
    
    blue <- .qbrms__package_blue()
    xvar <- pick_xvar(df, eff)
    
    # optional grouping factor (for interactions)
    known_cols    <- c(".draw", ".value", "estimate__", "lower__", "upper__", "..grp..", xvar)
    by_candidates <- setdiff(names(df), known_cols)
    by_candidates <- by_candidates[vapply(df[by_candidates], function(z) is.factor(z) || is.character(z), TRUE)]
    byvar <- if (length(by_candidates)) by_candidates[1] else NULL
    
    # ---- Spaghetti branch
    if (identical(kind, "spaghetti")) {
      if (!is.null(byvar)) {
        if (!("..grp.." %in% names(df))) {
          df$..grp.. <- interaction(df$.draw, df[[byvar]], drop = TRUE)
        }
        p <- ggplot2::ggplot(
          df,
          ggplot2::aes(x = .data[[xvar]], y = .data[[".value"]],
                       group = .data[["..grp.."]], colour = .data[[byvar]])
        ) +
          ggplot2::geom_line(alpha = 0.25, linewidth = 0.25)
      } else {
        p <- ggplot2::ggplot(
          df,
          ggplot2::aes(x = .data[[xvar]], y = .data[[".value"]], group = .data[[".draw"]])
        ) +
          ggplot2::geom_line(alpha = 0.25, linewidth = 0.25, colour = blue)
      }
      p <- p + ggplot2::labs(x = xvar, y = "Fitted value",
                             title = sprintf("Conditional effects: %s (spaghetti)", eff))
      return(p)
    }
    
    # ---- Summary branch (ribbons or slices)
    must_have <- c("estimate__", "lower__", "upper__")
    if (!all(must_have %in% names(df)))
      stop("Unrecognised conditional effects data format for plotting.")
    
    is_discrete <- ("slice" %in% names(df)) && is.factor(df[[xvar]])
    slice_title <- attr(df, "slice_name", exact = TRUE)
    
    if (is.null(byvar)) {
      # continuous x, single group: ribbon only (no outer line)
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data[[xvar]], y = .data[["estimate__"]])
      ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data[["lower__"]], ymax = .data[["upper__"]]),
          alpha  = 0.30,
          fill   = .qbrms__package_blue(),
          colour = NA   
        ) +
        ggplot2::geom_line(linewidth = 0.9, colour = .qbrms__package_blue())
      
    } else if (is_discrete) {
      # brms-style slices: points + intervals, dodged
      dodge <- ggplot2::position_dodge(width = 0.5)
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data[[xvar]], y = .data[["estimate__"]],
                     colour = .data[["slice"]])
      ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = .data[["lower__"]], ymax = .data[["upper__"]]),
          width = 0.15, position = dodge, linewidth = 0.5
        ) +
        ggplot2::geom_point(position = dodge, size = 2.6) +
        ggplot2::labs(colour = if (is.null(slice_title)) "slice" else slice_title)
    } else {
      # continuous x with factor groups: overlapping ribbons, no lines
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[xvar]])) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            ymin  = .data[["lower__"]],
            ymax  = .data[["upper__"]],
            fill  = .data[[byvar]],
            group = .data[[byvar]]
          ),
          alpha  = 0.30,
          colour = NA                     # <- no ribbon outline
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = .data[["estimate__"]], colour = .data[[byvar]]),
          linewidth = 0.9
        ) +
        ggplot2::labs(fill = byvar, colour = byvar)
      
    }
    
    p <- p + ggplot2::labs(x = xvar, y = "Fitted value",
                           title = sprintf("Conditional effects: %s", eff))
    return(p)
  }
  
  # ---- Multiple effects
  if (requireNamespace("patchwork", quietly = TRUE)) {
    plot_list <- vector("list", n_effects)
    for (i in seq_along(effect_names)) {
      plot_list[[i]] <- plot.qbrms_conditional_effects(x[i])
    }
    return(patchwork::wrap_plots(plot_list, ncol = min(2, n_effects)))
  } else {
    warning("Multiple effects detected. Install 'patchwork' for combined plots, or plot individual effects manually.")
    plot_list <- vector("list", n_effects)
    names(plot_list) <- effect_names
    for (i in seq_along(effect_names)) {
      plot_list[[i]] <- plot.qbrms_conditional_effects(x[i])
    }
    return(plot_list)
  }
}
