# Joint-draw utilities (fixed-effects covariance + rmvnorm + beta draws)

.qbrms_get_vcov_fixed <- function(object) {
  okM <- function(m, p, rn) is.matrix(m) && nrow(m) == p && ncol(m) == p &&
    (is.null(rownames(m)) || all(rownames(m) %in% rn)) &&
    (is.null(colnames(m)) || all(colnames(m) %in% rn))
  
  # Check for fit object first
  if (!"fit" %in% names(object) || is.null(object$fit)) {
    return(NULL)
  }
  
  fit <- object$fit
  summ <- fit$summary.fixed
  
  # Handle case where summary.fixed might be NULL
  if (is.null(summ)) {
    # Try to extract from marginals.fixed if available
    if ("marginals.fixed" %in% names(fit) && !is.null(fit$marginals.fixed)) {
      marg <- fit$marginals.fixed
      param_names <- names(marg)
      
      # Extract means and sds from marginals
      means <- sapply(marg, function(m) {
        if (is.list(m) && all(c("x", "y") %in% names(m))) {
          dx <- diff(m$x); dx <- c(dx[1], dx)
          sum(m$x * m$y * dx) / sum(m$y * dx)
        } else NA
      })
      
      sds <- sapply(marg, function(m) {
        if (is.list(m) && all(c("x", "y") %in% names(m))) {
          dx <- diff(m$x); dx <- c(dx[1], dx)
          mean_val <- sum(m$x * m$y * dx) / sum(m$y * dx)
          var_val <- sum((m$x - mean_val)^2 * m$y * dx) / sum(m$y * dx)
          sqrt(var_val)
        } else NA
      })
      
      if (!any(is.na(means)) && !any(is.na(sds))) {
        summ <- data.frame(mean = means, sd = sds, row.names = param_names)
      }
    }
    
    if (is.null(summ)) return(NULL)
  }
  
  rn <- rownames(summ); p <- length(rn)
  
  # Try standard INLA covariance sources
  cand <- list(
    fit$cov.fixed,
    if("misc" %in% names(fit)) fit$misc$cov.fixed else NULL,
    tryCatch(fit$sdreport$cov.fixed, error = function(e) NULL)
  )
  
  for (V in cand) {
    if (!is.null(V) && okM(V, p, rn)) {
      if (!is.null(rownames(V))) V <- V[rn, rn, drop = FALSE]
      attr(V, "source") <- "model_cov"
      return(V)
    }
  }
  
  # Try correlation + standard deviations
  R <- tryCatch(fit$misc$cor.fixed, error = function(e) NULL)
  sds <- as.numeric(summ[, "sd"])
  if (is.matrix(R) && nrow(R) == p && ncol(R) == p) {
    V <- diag(sds) %*% R %*% diag(sds)
    rownames(V) <- colnames(V) <- rn
    attr(V, "source") <- "model_cor+sd"
    return(V)
  }
  
  # OLS fallback - CRITICAL: This provides proper intercept-slope covariance
  f <- object$original_formula
  df <- object$data
  resp <- all.vars(f)[1]
  df[[resp]] <- suppressWarnings(as.numeric(df[[resp]]))
  
  # Robust OLS fitting
  Vlm <- tryCatch({
    ols_fit <- stats::lm(f, data = df)
    vcov(ols_fit)
  }, error = function(e) NULL)
  
  # Build covariance matrix with INLA standard deviations
  Vpad <- diag(sds^2)
  rownames(Vpad) <- colnames(Vpad) <- rn
  
  if (!is.null(Vlm)) {
    keep <- intersect(rn, colnames(Vlm))
    if (length(keep) > 0) {
      # Use OLS covariance structure with INLA marginal variances
      ols_sub <- Vlm[keep, keep, drop = FALSE]
      
      # Scale OLS correlations with INLA standard deviations
      if (length(keep) > 1) {
        ols_cor <- cov2cor(ols_sub)
        inla_sds <- sds[match(keep, rn)]
        scaled_cov <- diag(inla_sds) %*% ols_cor %*% diag(inla_sds)
        Vpad[keep, keep] <- scaled_cov
      } else {
        Vpad[keep, keep] <- sds[match(keep, rn)]^2
      }
      
      attr(Vpad, "source") <- "ols_vcov_fallback"
      return(Vpad)
    }
  }
  
  # Final fallback - diagonal only (this produces flat lines!)
  attr(Vpad, "source") <- "diag_sd_only"
  return(Vpad)
}

.qbrms_rmvnorm <- function(n, mean, sigma) {
  p <- length(mean)
  
  # Ensure sigma is symmetric and positive definite
  sigma <- as.matrix(sigma)
  if (nrow(sigma) != p || ncol(sigma) != p) {
    stop("Covariance matrix dimension mismatch")
  }
  
  # Make symmetric (important for numerical stability)
  sigma <- (sigma + t(sigma)) / 2
  
  # Add small ridge to diagonal for numerical stability
  sigma <- sigma + diag(1e-10, p)
  
  # Try Cholesky decomposition first (more stable when it works)
  L <- tryCatch({
    chol(sigma)
  }, error = function(e) {
    # If Cholesky fails, use eigendecomposition
    eigen_result <- eigen(sigma, symmetric = TRUE)
    eigenvalues <- pmax(eigen_result$values, 1e-10)  # Ensure positive
    eigenvectors <- eigen_result$vectors
    
    # Construct matrix square root: sigma = V * diag(lambda) * V'
    # So sigma^(1/2) = V * diag(sqrt(lambda))
    eigenvectors %*% diag(sqrt(eigenvalues), p)
  })
  
  # Generate standard normal random variables
  Z <- matrix(rnorm(n * p), n, p)
  
  # Transform to target distribution
  if (is.matrix(L) && nrow(L) == p) {
    # If L is from Cholesky, it's upper triangular
    if (all(L[lower.tri(L)] == 0)) {
      # Cholesky case: X = Z * L + mean
      samples <- Z %*% L
    } else {
      # Eigendecomposition case: X = Z * L + mean  
      samples <- Z %*% t(L)
    }
  } else {
    # Fallback: direct matrix multiplication
    samples <- Z %*% t(L)
  }
  
  # Add mean
  samples <- sweep(samples, 2, mean, "+")
  
  return(samples)
}

.qbrms_beta_draws <- function(object, X, ndraws = 200, V_override = NULL) {
  summ <- object$fit$summary.fixed
  if (is.null(summ)) stop("No summary.fixed available in fit.")
  bmean  <- as.numeric(summ[, "mean"])
  bnames <- rownames(summ); names(bmean) <- bnames
  
  X_use <- X
  miss <- setdiff(bnames, colnames(X_use))
  if (length(miss)) for (nm in miss) X_use[[nm]] <- 0
  X_use <- as.matrix(X_use[, bnames, drop = FALSE])
  
  V <- if (!is.null(V_override)) V_override else .qbrms_get_vcov_fixed(object)
  if (is.null(V)) {
    V <- diag(as.numeric(summ[, "sd"])^2)
    rownames(V) <- colnames(V) <- bnames
    attr(V, "source") <- "diag_sd_only"
  } else {
    V <- V[bnames, bnames, drop = FALSE]
  }
  
  # Ensure V is properly ordered and symmetric
  V <- as.matrix(V)
  V <- (V + t(V)) / 2  # Force symmetry
  
  # Generate samples with enhanced function
  beta <- .qbrms_rmvnorm(ndraws, bmean, V)  # ndraws x p
  colnames(beta) <- bnames
  
  # Verification: check if correlation is preserved (optional, for debugging)
  if (ncol(beta) > 1 && length(bnames) > 1) {
    expected_cor <- cov2cor(V)
    actual_cor <- cor(beta)
    max_diff <- max(abs(expected_cor - actual_cor))
    if (max_diff > 0.1) {
      warning("Large correlation discrepancy detected (max diff: ", round(max_diff, 3), 
              "). Consider numerical issues.")
    }
  }
  
  list(X = X_use, beta = beta, V_used = V)
}