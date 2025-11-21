# =============================================================================
# R/kfold_cv.R  -  Safe K-fold cross-validation for qbrms/qbrmO
# =============================================================================

#' K-fold cross-validation for qbrms models (ordinal and standard families)
#'
#' @description
#' Performs K-fold cross-validation either from a fitted model or from
#' `formula + data`. For **ordinal (cumulative/ordinal)** families, you can
#' choose the re-fit engine used inside CV: TMB (`qbrmO`) or a robust
#' fallback using `MASS::polr` that avoids TMB compilation in each fold.
#' Your original fitted model is unchanged.
#'
#' @param object Either a fitted qbrms/qbrmO object **or** a formula.
#' @param data Required only if `object` is a formula. Ignored if `object` is a fit.
#' @param family Optional family override (used if `object` is a formula; fits use their own).
#' @param K Number of folds (default 10).
#' @param folds Optional integer vector of length `nrow(data)` giving fold IDs.
#' @param seed Optional seed for stratified folds.
#' @param stratify Logical; stratify on response if factor/ordered (default TRUE).
#' @param parallel Logical; use `future.apply` if available (default FALSE).
#' @param workers Optional workers when parallel and no plan is set.
#' @param keep_fits Logical; keep per-fold fits (default FALSE).
#' @param engine Ordinal CV engine: `"auto"` (default), `"tmb"`, or `"polr"`.
#'   Only used for **ordinal** families during CV refits. `"auto"` uses
#'   `getOption("qbrms.kfold.ordinal_engine", "polr")`.
#' @param verbose Logical; brief progress (default TRUE).
#' @param ... Passed to `qbrms()` when refitting folds (non-ordinal or `engine="tmb"`).
#'
#' @return An object of class `qbrms_kfold` with ELPD, pointwise elpd, SE, etc.
#' @export
kfold_cv <- function(object,
                     data       = NULL,
                     family     = gaussian(),
                     K          = 10,
                     folds      = NULL,
                     seed       = NULL,
                     stratify   = TRUE,
                     parallel   = FALSE,
                     workers    = NULL,
                     keep_fits  = FALSE,
                     engine     = c("auto","tmb","polr"),
                     verbose    = TRUE,
                     ...) {
  
  call_ <- match.call()
  engine <- match.arg(engine)
  
  # --------------------------------------------------------------------------
  # Accept fit or formula+data
  # --------------------------------------------------------------------------
  if (inherits(object, c("qbrms_fit","qbrms_multinomial_fit","qbrms_prior_fit","qbrmO_fit","list"))) {
    fit        <- object
    ext        <- .extract_from_fit(fit)
    formula    <- ext$formula
    data_full  <- ext$data
    fam_name   <- ext$fam_name
    fam_obj    <- ext$fam_obj
  } else {
    if (!inherits(object, "formula")) {
      stop("First argument must be either a fitted qbrms/qbrmO model or a formula.", call. = FALSE)
    }
    if (is.null(data)) {
      stop("When providing a formula, you must also supply `data=`.", call. = FALSE)
    }
    formula    <- object
    data_full  <- data
    fam_name   <- tryCatch(tolower(extract_family_name(convert_family_to_inla(family))),
                           error = function(e) tolower(as.character(family$family %||% "gaussian")))
    fam_obj    <- family
  }
  
  if (!is.data.frame(data_full) || nrow(data_full) < 2L) {
    stop("`data` must be a data frame with at least two rows.", call. = FALSE)
  }
  
  resp_name <- all.vars(formula)[1]
  if (is.null(resp_name) || !nzchar(resp_name) || !(resp_name %in% names(data_full))) {
    stop("Could not identify the response variable in `data`.", call. = FALSE)
  }
  
  n <- nrow(data_full)
  
  # Folds
  if (!is.null(folds)) {
    if (length(folds) != n || anyNA(folds)) {
      stop("`folds` must be an integer vector of length nrow(data) with no NAs.", call. = FALSE)
    }
    K <- length(unique(folds))
  } else {
    folds <- .make_folds(data_full[[resp_name]], K = K, seed = seed, stratify = stratify)
  }
  
  # Decide ordinal engine for CV refits
  ord_engine <- if (engine == "auto") {
    getOption("qbrms.kfold.ordinal_engine", "polr")
  } else engine
  ord_engine <- match.arg(ord_engine, c("tmb","polr"))
  
  use_parallel <- isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE)
  if (parallel && !use_parallel && verbose) {
    message("Parallel requested but future.apply not available; running sequentially.")
  }
  
  # --------------------------------------------------------------------------
  # Per-fold runner
  # --------------------------------------------------------------------------
  runner <- function(k) {
    test_idx  <- which(folds == k)
    train_idx <- setdiff(seq_len(n), test_idx)
    train_dat <- data_full[train_idx, , drop = FALSE]
    test_dat  <- data_full[test_idx,  , drop = FALSE]
    
    # ORDINAL family: optional safe engine via MASS::polr
    if (tolower(fam_name) %in% c("cumulative","ordinal") && ord_engine == "polr") {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        return(list(ok = FALSE,
                    msg = "MASS package (polr) is required for ordinal CV with engine='polr'.",
                    test_index = test_idx))
      }
      # Ensure ordered response with consistent levels
      levs_all <- levels(as.factor(data_full[[resp_name]]))
      train_dat[[resp_name]] <- ordered(train_dat[[resp_name]], levels = levs_all)
      test_dat[[resp_name]]  <- ordered(test_dat[[resp_name]],  levels = levs_all)
      
      # Fit polr (logistic link)
      fit_k <- tryCatch(
        MASS::polr(formula, data = train_dat, method = "logistic", Hess = FALSE, model = FALSE),
        error = function(e) e
      )
      if (inherits(fit_k, "error")) {
        return(list(ok = FALSE, msg = fit_k$message, test_index = test_idx))
      }
      
      # Predicted category probabilities
      pr <- tryCatch(stats::predict(fit_k, newdata = test_dat, type = "probs"),
                     error = function(e) e)
      
      if (inherits(pr, "error")) {
        return(list(ok = FALSE, msg = pr$message, test_index = test_idx))
      }
      if (!is.matrix(pr)) pr <- as.matrix(pr)
      colnames(pr) <- colnames(pr) %||% levels(test_dat[[resp_name]])
      
      ytest <- test_dat[[resp_name]]
      levs  <- levels(ytest)
      if (!all(levs %in% colnames(pr))) {
        return(list(ok = FALSE, msg = "Could not align predicted probabilities with response levels.",
                    test_index = test_idx))
      }
      pick <- cbind(seq_len(nrow(test_dat)), as.integer(ytest))
      p    <- pmax(pr[pick], 1e-12)
      return(list(ok = TRUE, test_index = test_idx, lppd = log(p),
                  fit = if (keep_fits) fit_k else NULL))
    }
    
    # Non-ordinal or ordinal with TMB engine -> refit via qbrms
    # For safety on macOS, keep TMB single-threaded inside CV
    old_env <- Sys.getenv(c("TMB_NUM_THREADS","OMP_THREAD_LIMIT","TMB_OPENMP_ENABLE"), unset = NA)
    on.exit({
      # restore env
      for (nm in names(old_env)) {
        if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else Sys.setenv(nm = old_env[[nm]])
      }
    }, add = TRUE)
    Sys.setenv(TMB_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
    
    fit_k <- tryCatch(
      qbrms(formula, data = train_dat, family = fam_obj %||% gaussian(), verbose = FALSE, ...),
      error = function(e) e
    )
    if (inherits(fit_k, "error")) {
      return(list(ok = FALSE, msg = fit_k$message, test_index = test_idx))
    }
    
    # Held-out log predictive density
    ll_test <- tryCatch(
      .heldout_log_pred_density(fit_k, formula, train_dat, test_dat, fam_name, resp_name),
      error = function(e) e
    )
    if (inherits(ll_test, "error")) {
      return(list(ok = FALSE, msg = ll_test$message, test_index = test_idx))
    }
    
    out <- list(ok = TRUE, test_index = test_idx, lppd = ll_test)
    if (keep_fits) out$fit <- fit_k
    out
  }
  
  if (use_parallel) {
    if (!is.null(workers)) {
      try({
        if (is.null(future::plan("list")$strategy)) {
          future::plan(future::multisession, workers = workers)
        }
      }, silent = TRUE)
    }
    res_list <- future.apply::future_lapply(seq_len(K), runner)
  } else {
    res_list <- lapply(seq_len(K), runner)
  }
  
  bad <- which(vapply(res_list, function(z) identical(z$ok, FALSE), logical(1)))
  if (length(bad)) {
    msgs <- vapply(res_list[bad], function(z) z$msg, character(1))
    stop("K-fold failed on fold(s) ", paste(bad, collapse = ", "), ":\n  ",
         paste(msgs, collapse = "\n  "), call. = FALSE)
  }
  
  lppd <- numeric(n)
  for (z in res_list) lppd[z$test_index] <- z$lppd
  
  elpd_sum <- sum(lppd)
  se_elpd  <- sqrt(length(lppd) * stats::var(lppd))
  
  fits <- NULL
  if (keep_fits) {
    fits <- lapply(res_list, function(z) z$fit)
    names(fits) <- paste0("fold_", seq_len(K))
  }
  
  out <- structure(
    list(
      elpd_sum       = elpd_sum,
      elpd_pointwise = lppd,
      se_elpd        = se_elpd,
      folds          = folds,
      K              = K,
      family_used    = fam_name,
      ordinal_engine = if (tolower(fam_name) %in% c("cumulative","ordinal")) ord_engine else NA_character_,
      call           = call_,
      fits           = fits
    ),
    class = "qbrms_kfold"
  )
  
  if (verbose) {
    cat("K-fold CV completed. ELPD:", sprintf("%.3f", elpd_sum),
        "   SE(ELPD):", sprintf("%.3f", se_elpd))
    if (tolower(fam_name) %in% c("cumulative","ordinal")) {
      cat("   [ordinal engine:", ord_engine, "]")
    }
    cat("\n")
  }
  
  out
}

#' @export
print.qbrms_kfold <- function(x, ...) {
  cat("qbrms K-fold cross-validation\n")
  cat("Family:", x$family_used, "\n")
  if (!is.na(x$ordinal_engine)) cat("Ordinal engine:", x$ordinal_engine, "\n")
  cat("Folds :", x$K, "\n")
  cat("ELPD  :", sprintf("%.3f", x$elpd_sum), "  SE:", sprintf("%.3f", x$se_elpd), "\n")
  invisible(x)
}

# -----------------------------------------------------------------------------
# Internal helpers (unchanged except where noted)
# -----------------------------------------------------------------------------
.extract_from_fit <- function(fit) {
  form <- try(fit$original_formula, silent = TRUE)
  if (inherits(form, "try-error") || is.null(form)) form <- try(fit$formula, silent = TRUE)
  if (inherits(form, "try-error") || is.null(form) || !inherits(form, "formula")) {
    stop("Could not extract formula from the fitted object.", call. = FALSE)
  }
  
  dat <- try(fit$data, silent = TRUE)
  if (inherits(dat, "try-error") || is.null(dat) || !is.data.frame(dat)) {
    stop("Could not extract original data from the fitted object.", call. = FALSE)
  }
  
  fam_name <- "gaussian"; fam_obj <- gaussian()
  cand <- try(fit$family, silent = TRUE)
  if (!inherits(cand, "try-error") && !is.null(cand)) {
    nm <- try(tolower(extract_family_name(cand)), silent = TRUE)
    if (!inherits(nm, "try-error") && is.character(nm)) {
      fam_name <- nm
    } else if (is.list(cand) && !is.null(cand$family)) {
      fam_name <- tolower(as.character(cand$family))
    } else if (is.character(cand)) {
      fam_name <- tolower(cand[1])
    }
  }
  
  fam_obj <- switch(fam_name,
                    "gaussian"   = gaussian(),
                    "binomial"   = binomial(),
                    "poisson"    = poisson(),
                    "t"          = student(),
                    "student"    = student(),
                    "studentt"   = student(),
                    "lognormal"  = lognormal(),
                    "cumulative" = cumulative(),
                    fam_obj
  )
  
  list(formula = form, data = dat, fam_name = fam_name, fam_obj = fam_obj)
}

.make_folds <- function(y, K, seed = NULL, stratify = TRUE) {
  n <- length(y)
  if (!is.null(seed)) set.seed(seed)
  if (stratify && (is.factor(y) || is.ordered(y))) {
    levs <- levels(as.factor(y))
    folds <- integer(n)
    for (lv in levs) {
      idx <- which(as.factor(y) == lv)
      folds[idx] <- sample(rep_len(seq_len(K), length(idx)))
    }
    return(folds)
  } else {
    sample(rep_len(seq_len(K), n))
  }
}

.heldout_log_pred_density <- function(fit, formula, train_dat, test_dat, fam_name, resp_name) {
  fam_name <- tolower(fam_name)
  
  if (fam_name %in% c("cumulative","ordinal")) {
    probs <- .predict_ordinal_probs(fit, formula, train_dat, test_dat)
    ytest <- test_dat[[resp_name]]
    
    if (!is.ordered(ytest) && !is.factor(ytest)) {
      stop("Ordinal K-fold requires the response to be an ordered factor or factor.", call. = FALSE)
    }
    levs <- levels(ytest)
    if (!all(levs %in% colnames(probs))) {
      train_resp <- train_dat[[resp_name]]
      train_levels <- if (is.factor(train_resp) || is.ordered(train_resp)) levels(train_resp) else levs
      levs <- train_levels
      if (!all(levs %in% colnames(probs))) {
        stop("Could not align predicted probabilities with response levels.", call. = FALSE)
      }
    }
    probs <- probs[, levs, drop = FALSE]
    pick <- cbind(seq_len(nrow(test_dat)), as.integer(ytest))
    p <- pmax(probs[pick], 1e-12)
    return(log(p))
  }
  
  coefs <- .extract_fixed_means(fit)
  Xtest <- stats::model.matrix(formula, test_dat)
  eta   <- drop(Xtest %*% coefs)
  y     <- test_dat[[resp_name]]
  
  if (fam_name == "gaussian") {
    sd_y <- .extract_gaussian_sd(fit, formula, train_dat, coefs)
    mu   <- eta
    return(stats::dnorm(y, mean = mu, sd = sd_y, log = TRUE))
  }
  
  if (fam_name == "binomial") {
    if (!all(y %in% c(0, 1))) {
      stop("Binomial plug-in expects a 0/1 response for held-out scoring.", call. = FALSE)
    }
    p <- stats::plogis(eta)
    p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
    return(y * log(p) + (1 - y) * log(1 - p))
  }
  
  if (fam_name == "poisson") {
    lambda <- exp(eta)
    return(stats::dpois(y, lambda = lambda, log = TRUE))
  }
  
  if (fam_name == "t") {
    pars <- .extract_t_params(fit, formula, train_dat, coefs)
    z  <- (y - eta) / pars$sd
    return(stats::dt(z, df = pars$df, log = TRUE) - log(pars$sd))
  }
  
  if (fam_name == "lognormal") {
    pars <- .extract_lognormal_sd(fit, formula, train_dat, coefs)
    if (any(y <= 0, na.rm = TRUE)) {
      stop("Lognormal plug-in requires positive held-out responses.", call. = FALSE)
    }
    return(stats::dlnorm(y, meanlog = eta, sdlog = pars$sdlog, log = TRUE))
  }
  
  stop("Plug-in held-out scoring is implemented for cumulative/ordinal, gaussian, binomial, poisson, t, and lognormal.", call. = FALSE)
}

.predict_ordinal_probs <- function(fit, formula, train_dat, newdata) {
  pr <- try(stats::predict(fit, newdata = newdata, type = "probs"), silent = TRUE)
  if (!inherits(pr, "try-error") && is.matrix(pr)) return(pr)
  pr2 <- try(stats::predict(fit, newdata = newdata, type = "response"), silent = TRUE)
  if (!inherits(pr2, "try-error") && is.matrix(pr2)) return(pr2)
  
  beta <- .extract_fixed_means(fit)
  thr  <- .extract_thresholds(fit)
  if (is.null(thr) || length(thr) < 1L) {
    stop("Could not extract ordinal thresholds from fitted object.", call. = FALSE)
  }
  Xnew <- stats::model.matrix(formula, newdata)
  eta  <- as.numeric(drop(Xnew %*% beta))
  Kcat <- length(thr) + 1L
  Fk   <- matrix(NA_real_, nrow = length(eta), ncol = Kcat)
  for (k in seq_len(Kcat)) {
    if (k == Kcat) Fk[, k] <- 1.0 else Fk[, k] <- stats::plogis(thr[k] - eta)
  }
  Fkm1 <- cbind(0, Fk[, 1:(Kcat - 1), drop = FALSE])
  pk   <- pmax(Fk - Fkm1, 1e-12)
  
  resp_train <- model.frame(formula, train_dat)[[1]]
  levs <- if (is.factor(resp_train) || is.ordered(resp_train)) levels(resp_train) else paste0("cat", seq_len(Kcat))
  colnames(pk) <- levs
  pk
}

.extract_fixed_means <- function(fit) {
  cand <- try(fit$fit$summary.fixed, silent = TRUE)
  if (!inherits(cand, "try-error") && is.data.frame(cand) && "mean" %in% names(cand)) {
    cf <- cand$mean; names(cf) <- rownames(cand); return(cf)
  }
  cand2 <- try(fit$coef, silent = TRUE)
  if (!inherits(cand2, "try-error") && is.numeric(cand2)) return(cand2)
  stop("Could not extract fixed-effect means from fitted object.", call. = FALSE)
}

.extract_thresholds <- function(fit) {
  df1 <- try(fit$fit$summary.fixed, silent = TRUE)
  if (!inherits(df1, "try-error") && is.data.frame(df1) && "mean" %in% names(df1)) {
    rn <- tolower(rownames(df1)); is_thr <- grepl("cut|thresh|threshold", rn)
    if (any(is_thr)) return(as.numeric(df1$mean[is_thr]))
  }
  df2 <- try(fit$fit$summary.hyperpar, silent = TRUE)
  if (!inherits(df2, "try-error") && is.data.frame(df2) && "mean" %in% names(df2)) {
    rn <- tolower(rownames(df2)); is_thr <- grepl("cut|thresh|threshold", rn)
    if (any(is_thr)) return(as.numeric(df2$mean[is_thr]))
  }
  thr <- try(fit$thresholds, silent = TRUE)
  if (!inherits(thr, "try-error") && is.numeric(thr)) return(as.numeric(thr))
  NULL
}

.extract_gaussian_sd <- function(fit, formula, train_dat, coefs) {
  hp <- try(fit$fit$summary.hyperpar, silent = TRUE)
  if (!inherits(hp, "try-error") && is.data.frame(hp)) {
    rn <- rownames(hp)
    is_prec <- grepl("Precision.*Gaussian", rn, ignore.case = TRUE)
    if (any(is_prec)) {
      prec <- hp$mean[which(is_prec)[1]]
      if (is.finite(prec) && prec > 0) return(1 / sqrt(prec))
    }
  }
  Xtr <- stats::model.matrix(formula, train_dat)
  mu  <- drop(Xtr %*% coefs)
  ytr <- model.frame(formula, train_dat)[[1]]
  stats::sd(ytr - mu, na.rm = TRUE)
}

.extract_t_params <- function(fit, formula, train_dat, coefs) {
  hp <- try(fit$fit$summary.hyperpar, silent = TRUE)
  if (!inherits(hp, "try-error") && is.data.frame(hp)) {
    rn <- tolower(rownames(hp))
    i_df <- which(grepl("degrees|df", rn) & grepl("t|student", rn))
    i_pr <- which(grepl("precision", rn)     & grepl("t|student", rn))
    if (length(i_df)) df <- hp$mean[i_df[1]]
    if (length(i_pr)) sd <- 1 / sqrt(hp$mean[i_pr[1]])
    if (exists("df") && exists("sd")) return(list(df = as.numeric(df), sd = as.numeric(sd)))
  }
  Xtr <- stats::model.matrix(formula, train_dat)
  mu  <- drop(Xtr %*% coefs)
  ytr <- model.frame(formula, train_dat)[[1]]
  sd_f <- stats::sd(ytr - mu, na.rm = TRUE) * 1.25
  list(df = 7, sd = sd_f)
}

.extract_lognormal_sd <- function(fit, formula, train_dat, coefs) {
  hp <- try(fit$fit$summary.hyperpar, silent = TRUE)
  if (!inherits(hp, "try-error") && is.data.frame(hp)) {
    rn <- tolower(rownames(hp))
    i_pr <- which(grepl("precision", rn) & grepl("lognormal", rn))
    if (length(i_pr)) {
      sdlog <- 1 / sqrt(hp$mean[i_pr[1]])
      return(list(sdlog = as.numeric(sdlog)))
    }
  }
  Xtr  <- stats::model.matrix(formula, train_dat)
  mu   <- drop(Xtr %*% coefs)
  ytr  <- model.frame(formula, train_dat)[[1]]
  if (any(ytr <= 0, na.rm = TRUE)) return(list(sdlog = 0.5))
  sdlog <- stats::sd(log(ytr) - mu, na.rm = TRUE)
  list(sdlog = if (is.finite(sdlog) && sdlog > 0) sdlog else 0.5)
}

`%||%` <- function(x, y) if (is.null(x)) y else x