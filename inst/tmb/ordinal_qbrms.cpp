// =============================================================================
// inst/tmb/ordinal_qbrms.cpp - CORRECTED VERSION
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
  DATA_SCALAR(re_sd_shape);    // Random effects SD prior shape
  DATA_SCALAR(re_sd_rate);     // Random effects SD prior rate
  DATA_INTEGER(has_random_effects); // Whether to include random effects
  
  // Parameters
  PARAMETER_VECTOR(threshold_raw); // Raw threshold parameters
  PARAMETER_VECTOR(beta);          // Fixed effect coefficients
  PARAMETER(log_re_sd);            // Log of random effects standard deviation
  PARAMETER_VECTOR(u);             // Random effects (if present)
  
  // Initialize negative log-likelihood
  Type nll = 0;
  
  // Get dimensions
  int n_obs = y.size();
  int n_thresh = threshold_raw.size();
  int n_cats = n_thresh + 1;
  int n_groups = u.size();
  int n_beta = beta.size();
  
  // CRITICAL FIX 1: Proper threshold parameterization with strict ordering
  vector<Type> threshold(n_thresh);
  threshold(0) = threshold_raw(0);  // First threshold unconstrained
  for(int k = 1; k < n_thresh; k++) {
    // Force strict increasing order: threshold[k] = threshold[k-1] + exp(threshold_raw[k])
    threshold(k) = threshold(k-1) + exp(threshold_raw(k));
  }
  
  // Prior contributions for RAW thresholds (not constrained ones)
  for(int k = 0; k < n_thresh; k++) {
    nll -= dnorm(threshold_raw(k), thresh_mean(k), thresh_sd(k), true);
  }
  
  // Prior contributions for fixed effects
  for(int j = 0; j < n_beta; j++) {
    nll -= dnorm(beta(j), beta_mean(j), beta_sd(j), true);
  }
  
  // Random effects priors
  Type re_sd = exp(log_re_sd);
  if(has_random_effects > 0 && n_groups > 0) {
    Type re_precision = Type(1.0) / (re_sd * re_sd);
    nll -= dgamma(re_precision, re_sd_shape, Type(1.0)/re_sd_rate, true);
    nll -= log_re_sd;
    nll -= log(Type(2.0)) + log(re_precision);
    
    for(int i = 0; i < n_groups; i++) {
      nll -= dnorm(u(i), Type(0), re_sd, true);
    }
  }
  
  // CRITICAL FIX 2: Correct ordinal likelihood computation
  vector<Type> fitted_probs(n_obs);
  
  for(int i = 0; i < n_obs; i++) {
    // Linear predictor for fixed effects
    Type eta = Type(0);
    for(int j = 0; j < X.cols(); j++) {
      eta += X(i, j) * beta(j);
    }
    
    // Add random effects if present
    if(has_random_effects > 0 && n_groups > 0) {
      for(int g = 0; g < Z.cols(); g++) {
        if(g < n_groups) {
          eta += Z(i, g) * u(g);
        }
      }
    }
    
    // CRITICAL FIX 3: Correct cumulative probability computation
    // P(Y <= j) = logit^(-1)(threshold[j] - eta)
    vector<Type> cumprob(n_thresh + 2);
    cumprob(0) = Type(0.0);  // P(Y <= 0) = 0
    
    for(int k = 0; k < n_thresh; k++) {
      cumprob(k + 1) = invlogit(threshold(k) - eta);
    }
    cumprob(n_thresh + 1) = Type(1.0);  // P(Y <= K) = 1
    
    // CRITICAL FIX 4: Proper category probabilities
    // P(Y = j) = P(Y <= j) - P(Y <= j-1)
    vector<Type> prob(n_cats);
    for(int k = 0; k < n_cats; k++) {
      prob(k) = cumprob(k + 1) - cumprob(k);
      
      // Ensure numerical stability - probabilities must be positive
      prob(k) = (prob(k) > Type(1e-12)) ? prob(k) : Type(1e-12);
    }
    
    // Normalize to ensure sum = 1 (numerical safety)
    Type prob_sum = Type(0);
    for(int k = 0; k < n_cats; k++) {
      prob_sum += prob(k);
    }
    for(int k = 0; k < n_cats; k++) {
      prob(k) = prob(k) / prob_sum;
    }
    
    // Get observed category (1-indexed to 0-indexed)
    int obs_cat = asDouble(y(i)) - 1;
    
    // Safety check
    if(obs_cat < 0 || obs_cat >= n_cats) {
      continue;  // Skip invalid observations
    }
    
    // Store fitted probability and add to likelihood
    fitted_probs(i) = prob(obs_cat);
    nll -= log(prob(obs_cat));
  }
  
  // Reports for debugging
  REPORT(fitted_probs);
  REPORT(threshold);  // Report constrained thresholds
  if(has_random_effects > 0) {
    REPORT(re_sd);
  }
  
  // ADREPORT for standard errors
  ADREPORT(threshold);  // Report constrained thresholds
  ADREPORT(beta);
  if(has_random_effects > 0) {
    ADREPORT(re_sd);
  }
  
  return nll;
}
