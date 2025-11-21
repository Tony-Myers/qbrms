# =============================================================================
# R/compare_models.R
# =============================================================================

#' Compare qbrms models
#'
#' @description
#' Compares multiple fitted models using information criteria and simple
#' predictive metrics. Preference order for criterion = "auto" is:
#' LOO (from CPO) > WAIC > DIC. When information criteria are unavailable
#' it falls back to predictive metrics (RMSE/MAE).
#'
#' @param ... Two or more fitted model objects (qbrms_fit or qbrmO_fit)
#' @param criterion One of "auto","loo","waic","dic","all"
#' @param compare_predictions Logical; if TRUE, include RMSE/MAE comparison
#' @param weights Logical; if TRUE, compute weights when a single criterion is used
#' @return An object of class \code{"qbrms_comparison"}.
#' @export
compare_models <- function(...,
                           criterion = c("auto","loo","waic","dic","all"),
                           compare_predictions = TRUE,
                           weights = TRUE) {
  
  models <- list(...)
  if (length(models) < 2) stop("At least two models are required for comparison", call. = FALSE)
  
  model_names <- names(models)
  if (is.null(model_names) || !length(model_names)) {
    model_names <- paste0("Model", seq_along(models))
  } else {
    empty <- model_names == "" | is.na(model_names)
    if (any(empty)) model_names[empty] <- paste0("Model", which(empty))
  }
  names(models) <- model_names
  
  is_qbrms <- vapply(models, function(x) inherits(x, "qbrms_fit") || inherits(x, "qbrmO_fit"), logical(1))
  if (!all(is_qbrms)) {
    bad <- paste(model_names[!is_qbrms], collapse = ", ")
    stop("All objects must be qbrms model fits. Invalid: ", bad, call. = FALSE)
  }
  
  criterion <- match.arg(tolower(criterion[1]), c("auto","loo","waic","dic", "all"))
  
  # 1) Gather all available ICs
  ic <- .extract_information_criteria(models)
  
  # 2) Decide which criterion to use if "auto"
  criterion_used <- criterion
  if (criterion == "auto") {
    have <- function(col) isTRUE(col %in% colnames(ic$values)) && sum(is.finite(ic$values[[col]])) >= 2
    if      (have("loo"))  criterion_used <- "loo"
    else if (have("waic")) criterion_used <- "waic"
    else if (have("dic"))  criterion_used <- "dic"
    else                   criterion_used <- "dic"
  }
  
  # 3) Model weights if appropriate
  model_weights <- NULL
  if (weights && criterion_used != "all" && criterion_used %in% colnames(ic$values)) {
    model_weights <- .calculate_model_weights(ic$values, criterion_used)
  }
  
  # 4) Predictive comparison (fallback uses fitted values)
  pred_comparison <- NULL
  if (compare_predictions) pred_comparison <- .compare_predictions(models)
  
  # 5) Build comparison table
  comparison_table <- .create_comparison_table(ic, model_weights, pred_comparison, criterion_used)
  
  # 6) Best model by chosen criterion
  best_model <- NA_character_
  if (criterion_used %in% colnames(ic$values)) {
    vals <- ic$values[[criterion_used]]
    if (any(is.finite(vals))) best_model <- names(vals)[which.min(vals)]
  }
  
  structure(
    list(
      comparison_table = comparison_table,
      best_model       = best_model,
      weights          = model_weights,
      criteria         = ic,
      predictions      = pred_comparison,
      criterion_used   = criterion_used,
      model_names      = model_names,
      models           = models
    ),
    class = "qbrms_comparison"
  )
}

# ------------------------------- internals ------------------------------------

.extract_information_criteria <- function(models) {
  n  <- length(models)
  nm <- names(models)
  
  # Prefer namespace dispatch for our accessors
  call_ns <- function(fun, x) {
    if (exists(fun, where = asNamespace("qbrms"), inherits = FALSE)) {
      f <- get(fun, asNamespace("qbrms"))
      return(try(f(x), silent = TRUE))
    }
    if (exists(fun, mode = "function")) return(try(get(fun)(x), silent = TRUE))
    structure("not_found", class = "try-error")
  }
  
  as_num <- function(x, key = NULL) {
    if (inherits(x, "try-error")) return(NA_real_)
    if (is.numeric(x) && length(x) == 1) return(as.numeric(x))
    if (is.list(x)) {
      if (!is.null(key) && !is.null(x[[key]]) && is.numeric(x[[key]]))
        return(as.numeric(x[[key]][1]))
      for (k in c("looic","waic","dic","AIC","BIC","estimate","elpd_loo")) {
        if (!is.null(x[[k]]) && is.numeric(x[[k]])) return(as.numeric(x[[k]][1]))
      }
    }
    NA_real_
  }
  
  loo  <- rep(NA_real_, n)
  waic <- rep(NA_real_, n)
  dic  <- rep(NA_real_, n)
  aic  <- rep(NA_real_, n)
  bic  <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    m <- models[[i]]
    l <- call_ns("loo",  m);  loo[i]  <- as_num(l, "looic")
    w <- call_ns("waic", m);  waic[i] <- as_num(w, "waic")
    d <- call_ns("dic",  m);  dic[i]  <- as_num(d, "dic")
    a <- try(stats::AIC(m), silent = TRUE); aic[i] <- if (inherits(a, "try-error")) NA_real_ else as.numeric(a)
    b <- try(stats::BIC(m), silent = TRUE); bic[i] <- if (inherits(b, "try-error")) NA_real_ else as.numeric(b)
  }
  
  cols <- list(loo = loo, waic = waic, dic = dic, aic = aic, bic = bic)
  keep <- vapply(cols, function(v) any(is.finite(v)), logical(1))
  values <- as.data.frame(cols[keep], row.names = nm, check.names = FALSE)
  
  list(values = values)
}

.calculate_model_weights <- function(values_df, criterion) {
  if (!criterion %in% colnames(values_df)) return(NULL)
  v <- values_df[[criterion]]
  if (!any(is.finite(v))) return(NULL)
  d <- v - min(v, na.rm = TRUE)
  w <- exp(-0.5 * d); w <- w / sum(w, na.rm = TRUE); w[is.na(w)] <- 0
  stats::setNames(w, rownames(values_df))
}

.compare_predictions <- function(models) {
  nm <- names(models)
  RMSE <- rep(NA_real_, length(models))
  MAE  <- rep(NA_real_, length(models))
  
  get_formula <- function(m) {
    f <- try(stats::formula(m), silent = TRUE)
    if (!inherits(f, "try-error")) return(f)
    if (!is.null(m$formula)) return(m$formula)
    NULL
  }
  get_data <- function(m) {
    if (!is.null(m$data)) return(m$data)
    f <- get_formula(m)
    if (!is.null(f) && !is.null(environment(f)$data)) return(environment(f)$data)
    mf <- try(stats::model.frame(m), silent = TRUE)
    if (!inherits(mf, "try-error")) return(mf)
    NULL
  }
  
  for (i in seq_along(models)) {
    m  <- models[[i]]
    f  <- get_formula(m)
    df <- get_data(m)
    if (is.null(f) || is.null(df)) next
    
    resp <- all.vars(f)[1]
    if (!resp %in% names(df)) next
    y <- try(df[[resp]], silent = TRUE)
    if (inherits(y, "try-error")) next
    
    p <- try(stats::predict(m, type = "response"), silent = TRUE)
    if (inherits(p, "try-error") || length(p) != length(y)) next
    
    y <- suppressWarnings(as.numeric(y))
    p <- suppressWarnings(as.numeric(p))
    ok <- is.finite(y) & is.finite(p)
    if (!any(ok)) next
    
    err <- y[ok] - p[ok]
    RMSE[i] <- sqrt(mean(err^2))
    MAE[i]  <- mean(abs(err))
  }
  
  data.frame(Model = nm, RMSE = RMSE, MAE = MAE, row.names = NULL, check.names = FALSE)
}

.create_comparison_table <- function(ic_results, model_weights, pred_comparison, criterion_used) {
  tbl <- ic_results$values
  if (!is.null(model_weights) && length(model_weights)) {
    tbl$weight <- model_weights[rownames(tbl)]
  }
  if (!is.null(pred_comparison) && nrow(pred_comparison)) {
    pc <- pred_comparison
    rownames(pc) <- pc$Model
    pc$Model <- NULL
    for (nm in intersect(colnames(pc), c("RMSE","MAE"))) {
      tbl[[nm]] <- pc[rownames(tbl), nm]
    }
  }
  attr(tbl, "criterion_used") <- criterion_used
  tbl
}