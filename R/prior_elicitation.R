# =============================================================================
# R/prior_elicitation.R  (link-aware)
# =============================================================================

# Small helpers ---------------------------------------------------------------

.logit  <- function(p) log(p/(1-p))
.ilogit <- function(eta) 1/(1+exp(-eta))
.z_from_level <- function(level) stats::qnorm((1+level)/2)

# Map your family object/string to a canonical name and default link
.qbrms_family_info <- function(family) {
  nm <- tryCatch({
    if (is.character(family)) tolower(family)
    else if (inherits(family, "family")) tolower(family$family)
    else if (is.list(family) && !is.null(family$family)) tolower(family$family)
    else "gaussian"
  }, error = function(e) "gaussian")
  
  # default links used here
  link <- switch(nm,
                 "binomial"   = "logit",
                 "beta"       = "logit",
                 "poisson"    = "log",
                 "nbinomial"  = "log",
                 "negbinomial"= "log",
                 "gamma"      = "log",
                 "lognormal"  = "log",
                 "gaussian"   = "identity",
                 "normal"     = "identity",
                 "skew_normal"= "identity",
                 "student_t"  = "identity",
                 "t"          = "identity",
                 "sn"         = "identity",
                 "asymmetric_laplace" = "identity",
                 "multinomial" = "multinomial",
                 "identity"
  )
  list(name = nm, link = link)
}

# Public API ------------------------------------------------------------------

#' Build priors from elicited beliefs (GLM-aware)
#'
#' @title Prior Build from Beliefs
#' @param formula Model formula
#' @param data Data frame
#' @param family Model family
#' @param beliefs List of beliefs about parameters
#' @param outcome_location Expected outcome location
#' @param outcome_interval Expected outcome interval
#' @param outcome_level Confidence level for outcome interval
#' @param outcome_sd Outcome standard deviation
#' @param standardise Whether to standardise predictors
#' @param plausible_range Plausible range for outcomes
#' @param target_coverage Target coverage probability
#' @param tune Whether to tune priors
#' @param seed Random seed
#' @export
prior_build_from_beliefs <- function(formula,
                                     data,
                                     family = gaussian(),
                                     beliefs = list(),
                                     outcome_location = NULL,
                                     outcome_interval = NULL,
                                     outcome_level = 0.95,
                                     outcome_sd = NULL,
                                     standardise = TRUE,
                                     plausible_range = NULL,
                                     target_coverage = 0.8,
                                     tune = FALSE,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.character(formula)) formula <- stats::as.formula(formula)
  
  fam <- .qbrms_family_info(family)
  if (identical(fam$name, "multinomial")) {
    stop("prior_build_from_beliefs does not yet support multinomial models.")
  }
  
  # Identify predictors in RHS -------------------------------------------------
  f_terms <- all.vars(formula)
  if (length(f_terms) < 2) stop("Formula must contain at least one predictor")
  resp <- f_terms[1L]
  preds_all <- unique(f_terms[-1L])
  
  # Standardise if requested, update formula ----------------------------------
  used_preds <- preds_all
  data_out <- data
  frm_out  <- formula
  ref_x    <- list()
  if (standardise) {
    for (v in preds_all) {
      if (is.numeric(data[[v]])) {
        m <- mean(data[[v]], na.rm = TRUE)
        s <- stats::sd(data[[v]], na.rm = TRUE)
        data_out[[paste0(v, "_s")]] <- if (!is.finite(s) || s <= 0) data[[v]] - m else as.numeric((data[[v]] - m) / s)
        ref_x[[v]] <- list(mean = 0, sd = if (!is.finite(s) || s <= 0) NA_real_ else 1, used = paste0(v, "_s"))
      } else {
        ref_x[[v]] <- list(mean = NA_real_, sd = NA_real_, used = v)
      }
    }
    rhs <- deparse(formula[[3L]])
    for (v in preds_all) {
      if (paste0(v, "_s") %in% names(data_out)) {
        rhs <- gsub(paste0("\\b", v, "\\b"), paste0(v, "_s"), rhs)
        used_preds[used_preds == v] <- paste0(v, "_s")
      }
    }
    frm_out <- stats::as.formula(paste(deparse(formula[[2L]]), "~", rhs))
  } else {
    for (v in preds_all) {
      if (is.numeric(data[[v]])) {
        ref_x[[v]] <- list(mean = mean(data[[v]], na.rm = TRUE),
                           sd   = stats::sd(data[[v]], na.rm = TRUE),
                           used = v)
      } else {
        ref_x[[v]] <- list(mean = NA_real_, sd = NA_real_, used = v)
      }
    }
  }
  
  # Intercept centre on response scale ----------------------------------------
  if (is.null(outcome_location)) {
    if (!is.null(plausible_range) && length(plausible_range) == 2L) {
      outcome_location <- mean(range(plausible_range))
    } else if (resp %in% names(data) && is.numeric(data[[resp]])) {
      outcome_location <- mean(data[[resp]], na.rm = TRUE)
    } else {
      # Default neutral values by family
      outcome_location <- switch(fam$link,
                                 "logit" = 0.5,
                                 "log"   = 1.0,
                                 0.0
      )
    }
  }
  
  # Translate beliefs to slope priors -----------------------------------------
  priors <- list()
  prior_rows <- list()
  
  for (nm in names(beliefs)) {
    bel <- beliefs[[nm]]
    used_name <- if (standardise) paste0(nm, "_s") else nm
    if (!(used_name %in% used_preds)) {
      warning("Belief supplied for predictor '", nm, "' which is not in the current formula; skipping")
      next
    }
    dx     <- bel$dx %||% 1
    level  <- bel$level %||% 0.95
    z      <- .z_from_level(level)
    
    # Optional more natural specifications for non-identity links
    odds_ratio_interval <- bel$odds_ratio_interval
    ratio_interval      <- bel$ratio_interval
    baseline            <- bel$baseline
    
    if (fam$link == "identity") {
      # Effect on y is linear: delta_y ~ Normal(mu, sd) for +dx
      if (is.null(bel$interval)) stop("Provide 'interval' for identity-link beliefs.")
      L <- min(bel$interval); U <- max(bel$interval)
      mu_beta <- ((L + U) / 2) / dx
      sd_beta <- (U - L) / (2 * z * dx)
      dist    <- bel$dist %||% "normal"; df <- bel$df %||% 7
      if (identical(dist, "student_t")) {
        priors[[length(priors)+1L]] <- prior(student_t_prior(df, mu_beta, sd_beta), class = "b", coef = used_name)
        prior_rows[[length(prior_rows)+1L]] <- data.frame(term = used_name, class = "b", dist = "student_t",
                                                          mean = mu_beta, sd = sd_beta, df = df)
      } else {
        priors[[length(priors)+1L]] <- prior(normal(mu_beta, sd_beta), class = "b", coef = used_name)
        prior_rows[[length(prior_rows)+1L]] <- data.frame(term = used_name, class = "b", dist = "normal",
                                                          mean = mu_beta, sd = sd_beta, df = NA_real_)
      }
      
    } else if (fam$link == "logit") {
      # Prefer odds-ratio elicitation if supplied
      if (!is.null(odds_ratio_interval)) {
        L <- min(odds_ratio_interval); U <- max(odds_ratio_interval)
        mu_beta <- log(sqrt(L * U)) / dx
        sd_beta <- (log(U) - log(L)) / (2 * z * dx)
      } else {
        # Otherwise use small-change approximation around a baseline probability
        p0 <- if (!is.null(baseline)) baseline else outcome_location
        if (!is.finite(p0) || p0 <= 0 || p0 >= 1) p0 <- 0.5
        if (is.null(bel$interval)) stop("Provide 'interval' for probability change on logit link, or 'odds_ratio_interval'.")
        L <- min(bel$interval); U <- max(bel$interval)
        mu_beta <- ((L + U) / 2) / (dx * p0 * (1 - p0))
        sd_beta <- (U - L) / (2 * z * dx * p0 * (1 - p0))
      }
      priors[[length(priors)+1L]] <- prior(normal(mu_beta, sd_beta), class = "b", coef = used_name)
      prior_rows[[length(prior_rows)+1L]] <- data.frame(term = used_name, class = "b", dist = "normal",
                                                        mean = mu_beta, sd = sd_beta, df = NA_real_)
      
    } else if (fam$link == "log") {
      # Prefer multiplicative ratio elicitation on the mean
      if (is.null(ratio_interval) && is.null(bel$interval)) {
        stop("Provide 'ratio_interval' (preferred) or 'interval' for change on the mean with log link.")
      }
      if (!is.null(ratio_interval)) {
        L <- min(ratio_interval); U <- max(ratio_interval)
        if (L <= 0) stop("ratio_interval must be positive.")
        mu_beta <- log(sqrt(L * U)) / dx
        sd_beta <- (log(U) - log(L)) / (2 * z * dx)
      } else {
        # Small-change approximation: delta_mu ≈ mu0 * beta * dx
        mu0 <- if (!is.null(baseline)) baseline else outcome_location
        if (!is.finite(mu0) || mu0 <= 0) mu0 <- 1
        L <- min(bel$interval); U <- max(bel$interval)
        mu_beta <- ((L + U) / 2) / (dx * mu0)
        sd_beta <- (U - L) / (2 * z * dx * mu0)
      }
      priors[[length(priors)+1L]] <- prior(normal(mu_beta, sd_beta), class = "b", coef = used_name)
      prior_rows[[length(prior_rows)+1L]] <- data.frame(term = used_name, class = "b", dist = "normal",
                                                        mean = mu_beta, sd = sd_beta, df = NA_real_)
    } else {
      stop("Link '", fam$link, "' is not currently supported in the elicitation helper.")
    }
  }
  
  # Intercept prior on response scale, then expressed on the linear predictor --
  # We keep a Normal prior on the linear predictor that corresponds to a chosen
  # mean and scale on the response scale under the link.
  # Outcome centre:
  y_loc <- outcome_location
  # Outcome SD from interval or SD:
  if (!is.null(outcome_interval) && length(outcome_interval) == 2L) {
    int_sd_resp <- diff(range(outcome_interval)) / (2 * .z_from_level(outcome_level))
  } else {
    int_sd_resp <- outcome_sd %||% {
      if (!is.null(plausible_range) && length(plausible_range) == 2L) diff(range(plausible_range))/4 else 10
    }
  }
  # Map response-scale location and SD to linear predictor SD at the centre
  if (fam$link == "identity") {
    mu_eta <- y_loc
    sd_eta <- int_sd_resp
  } else if (fam$link == "logit") {
    p0     <- min(max(y_loc, 1e-6), 1-1e-6)
    mu_eta <- .logit(p0)
    # delta p ≈ p0(1-p0) delta eta  =>  sd_eta ≈ sd_p / [p0(1-p0)]
    sd_eta <- int_sd_resp / (p0 * (1 - p0))
  } else if (fam$link == "log") {
    mu0    <- if (is.finite(y_loc) && y_loc > 0) y_loc else 1
    mu_eta <- log(mu0)
    # delta mu ≈ mu0 * delta eta   =>  sd_eta ≈ sd_mu / mu0
    sd_eta <- int_sd_resp / mu0
  } else {
    mu_eta <- y_loc
    sd_eta <- int_sd_resp
  }
  
  priors[[length(priors)+1L]] <- prior(normal(mu_eta, sd_eta), class = "Intercept")
  prior_rows[[length(prior_rows)+1L]] <- data.frame(term="(Intercept)", class="Intercept", dist="normal",
                                                    mean=mu_eta, sd=sd_eta, df=NA_real_)
  
  if (length(priors) == 0L) warning("No priors constructed; returning defaults")
  prior_list <- do.call(c, priors)
  
  # Optional quick prior-predictive check + mild tuning ------------------------
  diag_obj <- NULL
  if (!is.null(plausible_range) && length(plausible_range) == 2L) {
    fit0 <- qbrm(formula = frm_out, data = data_out, family = family,
                 prior = prior_list, sample_prior = "only", verbose = FALSE)
    diag_obj <- prior_pp_diagnostics(fit0, include_observed = FALSE,
                                     plausible_lower  = min(plausible_range),
                                     plausible_upper  = max(plausible_range))
    if (isTRUE(tune)) {
      cov <- diag_obj$plausible_coverage
      tgt <- target_coverage
      iter <- 0L
      while (is.finite(cov) && iter < 6L && abs(cov - tgt) > 0.02) {
        scale_factor <- sqrt(cov / tgt)
        # rebuild scaled priors
        for (i in seq_along(priors)) {
          row <- prior_rows[[i]]
          if (is.na(row$sd)) next
          new_sd <- max(1e-6, as.numeric(row$sd) * scale_factor)
          prior_rows[[i]]$sd <- new_sd
          if (row$class == "Intercept") {
            priors[[i]] <- prior(normal(row$mean, new_sd), class = "Intercept")
          } else if (row$dist == "student_t") {
            priors[[i]] <- prior(student_t_prior(row$df, row$mean, new_sd), class = "b", coef = row$term)
          } else {
            priors[[i]] <- prior(normal(row$mean, new_sd), class = "b", coef = row$term)
          }
        }
        prior_list <- do.call(c, priors)
        fit0 <- qbrm(formula = frm_out, data = data_out, family = family,
                     prior = prior_list, sample_prior = "only", verbose = FALSE)
        diag_obj <- prior_pp_diagnostics(fit0, include_observed = FALSE,
                                         plausible_lower  = min(plausible_range),
                                         plausible_upper  = max(plausible_range))
        cov  <- diag_obj$plausible_coverage
        iter <- iter + 1L
      }
    }
  }
  
  prior_df <- do.call(rbind, prior_rows); rownames(prior_df) <- NULL
  
  out <- list(
    priors      = prior_list,
    formula     = frm_out,
    data        = data_out,
    reference   = list(standardised = standardise, outcome_location = outcome_location,
                       family = fam$name, link = fam$link),
    diagnostics = diag_obj,
    summary     = prior_df
  )
  class(out) <- "qbrms_prior_build"
  out
}

#' @export
print.qbrms_prior_build <- function(x, digits = 3, ...) {
  cat("Prior build from elicited beliefs (link-aware)\n")
  cat("  Family:", x$reference$family, " Link:", x$reference$link, "\n")
  cat("  Standardised predictors:", if (isTRUE(x$reference$standardised)) "yes" else "no", "\n")
  cat(sprintf("  Intercept centre (response scale): %s\n",
              format(round(x$reference$outcome_location, digits), nsmall = digits)))
  cat("\nSuggested priors:\n")
  df <- x$summary
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    if (row$class == "Intercept") {
      cat(sprintf("  Intercept ~ Normal(mean=%s, sd=%s)\n",
                  format(round(row$mean, digits), nsmall = digits),
                  format(round(row$sd,   digits), nsmall = digits)))
    } else if (identical(row$dist, "student_t")) {
      cat(sprintf("  %s: b ~ Student-t(df=%s, mean=%s, scale=%s)\n",
                  row$term,
                  format(row$df, nsmall = 0),
                  format(round(row$mean, digits), nsmall = digits),
                  format(round(row$sd,   digits), nsmall = digits)))
    } else {
      cat(sprintf("  %s: b ~ Normal(mean=%s, sd=%s)\n",
                  row$term,
                  format(round(row$mean, digits), nsmall = digits),
                  format(round(row$sd,   digits), nsmall = digits)))
    }
  }
  if (!is.null(x$diagnostics)) {
    cat("\nPrior predictive coverage of the user-declared plausible outcome range:\n")
    pb <- x$diagnostics$plausible_bounds
    cov_str <- if (is.finite(x$diagnostics$plausible_coverage))
      format(round(x$diagnostics$plausible_coverage, 3), nsmall = 3) else "NA"
    cat(sprintf("  Range [%s, %s]   coverage = %s   verdict = %s\n",
                as.character(pb[["lower"]]), as.character(pb[["upper"]]),
                cov_str, x$diagnostics$verdict$status))
  }
  invisible(x)
}

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# -----------------------------------------------------------------------------
# Turn a qbrms_prior_build into qbrms-brms-style prior() code
# -----------------------------------------------------------------------------

#' Format priors as qbrms prior() code
#' @param build A 'qbrms_prior_build' returned by prior_build_from_beliefs()
#' @param object_name Name of the object on the left-hand side (default "priors")
#' @param digits Number of decimal places to print
#' @param include_comments Logical; if TRUE, prepend a short comment header
#' @return A single character string containing formatted R code
#' @export
prior_code <- function(build,
                       object_name = "priors",
                       digits = 3,
                       include_comments = TRUE) {
  if (!inherits(build, "qbrms_prior_build")) {
    stop("build must be a 'qbrms_prior_build' object (from prior_build_from_beliefs())")
  }
  df <- build$summary
  if (is.null(df) || nrow(df) == 0) stop("No prior summary available in 'build'")
  
  fmt <- function(x) format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE)
  esc <- function(s) gsub('"', '\\"', as.character(s), fixed = TRUE)
  
  lines <- character(0)
  
  if (isTRUE(include_comments)) {
    hdr <- sprintf(
      "# Priors for qbrms (%s link); coefficients are on the linear predictor (link) scale",
      build$reference$link %||% "identity"
    )
    lines <- c(lines, hdr)
  }
  
  lines <- c(lines, sprintf("%s <- c(", object_name))
  
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    if (identical(row$class, "Intercept")) {
      rhs <- sprintf('prior(normal(%s, %s), class = "Intercept")',
                     fmt(row$mean), fmt(row$sd))
    } else if (identical(row$dist, "student_t")) {
      rhs <- sprintf('prior(student_t_prior(%s, %s, %s), class = "b", coef = "%s")',
                     fmt(row$df), fmt(row$mean), fmt(row$sd), esc(row$term))
    } else {
      rhs <- sprintf('prior(normal(%s, %s), class = "b", coef = "%s")',
                     fmt(row$mean), fmt(row$sd), esc(row$term))
    }
    comma <- if (i < nrow(df)) "," else ""
    lines <- c(lines, paste0("  ", rhs, comma))
  }
  
  lines <- c(lines, ")")
  paste(lines, collapse = "\n")
}

#' Print method for qbrms_prior_code objects
#'
#' @param x A qbrms_prior_code object
#' @param ... Additional arguments passed to cat
#' @export
print.qbrms_prior_code <- function(x, ...) {
  cat(as.character(x), sep = "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Pretty printer for qbrms_prior_build
# ---------------------------------------------------------------------------

#' @export
print.qbrms_prior_build <- function(x,
                                    digits = 3,
                                    show_data = FALSE,
                                    show_code = TRUE,
                                    code_object_name = "priors",
                                    max_terms = 12,
                                    ...) {
  # helper
  fmt <- function(v) format(round(as.numeric(v), digits), nsmall = digits, trim = TRUE)
  deparse1 <- function(z) paste(deparse(z), collapse = "")
  
  cat("Prior build from elicited beliefs\n")
  cat("  Family:", x$reference$family, "  Link:", x$reference$link, "\n", sep = "")
  cat("  Standardised predictors:", if (isTRUE(x$reference$standardised)) "yes" else "no", "\n", sep = "")
  cat("  Formula: ", deparse1(x$formula), "\n", sep = "")
  
  # Summary table of priors
  df <- x$summary
  if (is.null(df) || nrow(df) == 0L) {
    cat("\n<no priors found>\n")
    return(invisible(x))
  }
  
  cat("\nSuggested priors (on the linear predictor scale):\n")
  # Build compact lines
  make_line <- function(row) {
    if (identical(row$class, "Intercept")) {
      sprintf("  %-14s ~ Normal(mean=%s, sd=%s)",
              "(Intercept)", fmt(row$mean), fmt(row$sd))
    } else if (identical(row$dist, "student_t")) {
      sprintf("  %-14s ~ Student-t(df=%s, mean=%s, scale=%s)",
              row$term, fmt(row$df), fmt(row$mean), fmt(row$sd))
    } else {
      sprintf("  %-14s ~ Normal(mean=%s, sd=%s)",
              row$term, fmt(row$mean), fmt(row$sd))
    }
  }
  
  n_show <- min(nrow(df), max_terms)
  for (i in seq_len(n_show)) cat(make_line(df[i, ]), "\n")
  if (nrow(df) > n_show) {
    cat(sprintf("  ... and %d more\n", nrow(df) - n_show))
  }
  
  # Diagnostics summary (if available)
  if (!is.null(x$diagnostics)) {
    pb <- x$diagnostics$plausible_bounds
    cov <- x$diagnostics$plausible_coverage
    ver <- x$diagnostics$verdict$status
    cat("\nPrior-predictive coverage of plausible outcome range:\n")
    cat(sprintf("  Range [%s, %s]   coverage = %s   verdict = %s\n",
                as.character(pb[["lower"]]),
                as.character(pb[["upper"]]),
                if (is.finite(cov)) format(round(cov, 3), nsmall = 3) else "NA",
                ver))
  }
  
  # Ready-to-copy code block
  if (isTRUE(show_code)) {
    cat("\nCopy-pasteable qbrms prior() code:\n")
    code <- prior_code(x, object_name = code_object_name, digits = digits, include_comments = FALSE)
    cat(code, "\n", sep = "")
  }
  
  # Optional: small peek at the standardised columns only
  if (isTRUE(show_data)) {
    cat("\nData (first 6 rows, standardised predictors only):\n")
    std_cols <- grep("(_s)$", names(x$data), value = TRUE)
    if (length(std_cols)) print(utils::head(x$data[std_cols]))
    else cat("  <no standardised predictors found>\n")
  }
  
  invisible(x)
}
