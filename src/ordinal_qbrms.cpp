// =============================================================================
// src/ordinal_qbrms.cpp - TMB ordinal regression template for qbrms
// TMB C++ template for ordinal regression with proper threshold handling
// =============================================================================
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs
  DATA_IVECTOR(y);             // Ordinal response (1, 2, ..., K)
  DATA_MATRIX(X);              // Fixed effects design matrix
  DATA_MATRIX(Z);              // Random effects design matrix
  DATA_VECTOR(thresh_mean);    // Prior means for thresholds
  DATA_VECTOR(thresh_sd);      // Prior SDs for thresholds
  DATA_VECTOR(beta_mean);      // Prior means for coefficients
  DATA_VECTOR(beta_sd);        // Prior SDs for coefficients
  DATA_SCALAR(re_sd_shape);    // Random effect sd Gamma shape
  DATA_SCALAR(re_sd_rate);     // Random effect sd Gamma rate
  DATA_INTEGER(has_random_effects); // 1 if random effects present, 0 otherwise

  // Parameter inputs
  PARAMETER_VECTOR(threshold_raw);  // Unconstrained thresholds
  PARAMETER_VECTOR(beta);           // Regression coefficients
  PARAMETER(log_re_sd);             // Log random effect standard deviation
  PARAMETER_VECTOR(u);              // Random effects

  int n = y.size();
  int n_beta = beta.size();
  int n_thresh = threshold_raw.size();

  // --------------------------------------------------------------------------
  // Constrain thresholds to be ordered using cumulative sum of exponentials
  // --------------------------------------------------------------------------
  vector<Type> threshold(n_thresh);
  if (n_thresh > 0) {
    threshold(0) = threshold_raw(0);
    for (int k = 1; k < n_thresh; k++) {
      threshold(k) = threshold(k - 1) + exp(threshold_raw(k));
    }
  }

  // --------------------------------------------------------------------------
  // Random effect standard deviation (sd) with Gamma prior on precision
  // --------------------------------------------------------------------------
  Type re_sd = exp(log_re_sd);

  // Negative log-likelihood
  Type nll = 0.0;

  // --------------------------------------------------------------------------
  // Priors for thresholds (Gaussian)
  // --------------------------------------------------------------------------
  for (int k = 0; k < n_thresh; k++) {
    Type z = (threshold(k) - thresh_mean(k)) / thresh_sd(k);
    nll -= dnorm(z, Type(0.0), Type(1.0), true);
  }

  // Priors for fixed effects (Gaussian)
  for (int j = 0; j < n_beta; j++) {
    Type z = (beta(j) - beta_mean(j)) / beta_sd(j);
    nll -= dnorm(z, Type(0.0), Type(1.0), true);
  }

  // Prior for random effect standard deviation (Gamma on 1/sd^2)
  if (has_random_effects > 0) {
    Type prec = Type(1.0) / (re_sd * re_sd);
    nll -= dgamma(prec, re_sd_shape, re_sd_rate, true);
  }

  // --------------------------------------------------------------------------
  // Likelihood contribution
  // --------------------------------------------------------------------------
  for (int i = 0; i < n; i++) {
    int yi = y(i);
    if (yi <= 0) continue;  // handle any missing or invalid values defensively

    // Linear predictor: fixed + random
    Type eta = Type(0.0);
    for (int j = 0; j < n_beta; j++) {
      eta += X(i, j) * beta(j);
    }

    if (has_random_effects > 0 && u.size() > 0) {
      // Random effects via Z * u
      Type re_contrib = Type(0.0);
      for (int j = 0; j < u.size(); j++) {
        re_contrib += Z(i, j) * u(j);
      }
      eta += re_contrib;
    }

    // Cumulative probabilities (logit link)
    Type logit_lower = -INFINITY;
    Type logit_upper = INFINITY;

    if (yi > 1) {
      logit_lower = threshold(yi - 2) - eta;
    }
    if (yi <= n_thresh) {
      logit_upper = threshold(yi - 1) - eta;
    }

    Type p_lower = (yi > 1)        ? invlogit(logit_lower) : Type(0.0);
    Type p_upper = (yi <= n_thresh) ? invlogit(logit_upper) : Type(1.0);

    Type p_cat = p_upper - p_lower;
    p_cat = CppAD::CondExpGt(p_cat, Type(1e-16), p_cat, Type(1e-16)); // avoid log(0)
    nll -= log(p_cat);
  }

  // Random effects likelihood contribution
  if (has_random_effects > 0 && u.size() > 0) {
    for (int j = 0; j < u.size(); j++) {
      nll -= dnorm(u(j), Type(0.0), re_sd, true);
    }
  }

  // --------------------------------------------------------------------------
  // REPORT / ADREPORT for posterior summaries
  // --------------------------------------------------------------------------
  REPORT(threshold);  // Report constrained thresholds
  if (has_random_effects > 0) {
    REPORT(re_sd);
  }

  // ADREPORT for standard errors
  ADREPORT(threshold);  // Report constrained thresholds
  ADREPORT(beta);
  if (has_random_effects > 0) {
    ADREPORT(re_sd);
  }

  return nll;
}
