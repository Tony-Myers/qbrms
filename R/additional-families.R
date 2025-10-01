# =============================================================================
# R/additional-families.R - Extended Family Constructors for qbrms
# =============================================================================
# Additional INLA family constructors that extend qbrms beyond basic families
# NOTE: Family conversion logic is in families.R - this file contains ONLY constructors

#' @title Additional Statistical Families for qbrms
#' @description Extended collection of statistical family constructors 
#' optimized for INLA integration and CRAN compatibility
#' @name additional_families
#' @keywords internal
NULL

# =============================================================================
# SECTION 1: COMPOSITIONAL AND CONSTRAINED FAMILIES
# =============================================================================

#' Simplex Family for Compositional Data
#'
#' @description
#' The Simplex family for modeling data constrained to the unit interval (0,1)
#' with support for precision parameter modeling. Particularly useful for 
#' proportions, percentages, and beta-like regression when more flexibility 
#' is needed than standard beta regression.
#'
#' @param link Link function for the mean (default: "logit")
#' @param link.precision Link function for precision parameter (default: "log")
#'
#' @return A family object for use with qbrms()
#' 
#' @details
#' The Simplex distribution is parameterized with a mean μ ∈ (0,1) and precision 
#' parameter φ > 0. It provides more flexible tail behavior than beta regression 
#' and is particularly useful for highly skewed proportional data.
#' 
#' Requires response variable to be in interval (0,1). Values exactly 0 or 1 
#' should be adjusted slightly (e.g., using (y * (n-1) + 0.5) / n transformation).
#'
#' @examples
#' \dontrun{
#' # Basic simplex regression
#' fit <- qbrms(proportion ~ treatment, data = data, family = simplex())
#' 
#' # With precision modeling
#' fit <- qbrms(bf(prop ~ predictor, phi ~ group), data = data, family = simplex())
#' }
#'
#' @export
simplex <- function(link = "logit", link.precision = "log") {
  structure(list(
    family = "simplex",
    link = link,
    link.precision = link.precision
  ), class = "family")
}

#' Beta Binomial Family for Overdispersed Binary Data
#'
#' @description
#' Beta-binomial family for modeling overdispersed binomial data. Accounts for
#' extra-binomial variation through a beta mixing distribution.
#'
#' @param link Link function for probability parameter (default: "logit")
#'
#' @return A family object for use with qbrms()
#' 
#' @details
#' Useful when binomial models show overdispersion. The beta-binomial allows
#' the success probability to vary according to a beta distribution, naturally
#' modeling correlation within clusters or overdispersion.
#'
#' @examples
#' \dontrun{
#' # Overdispersed binary data
#' fit <- qbrms(cbind(successes, failures) ~ treatment, 
#'              data = clinical_data, family = beta_binomial())
#' }
#'
#' @export
beta_binomial <- function(link = "logit") {
  structure(list(
    family = "beta_binomial",
    link = link
  ), class = "family")
}

#' Alternative Beta Parameterizations
#'
#' @description
#' Alternative parameterizations of the beta distribution for specialized use cases.
#' 
#' @param link Link function (default: "logit")
#'
#' @details
#' - beta0: Beta distribution parameterized for boundary behavior
#' - beta1: Alternative beta parameterization for different tail behavior
#' - logitbeta: Beta with explicit logit parameterization
#'
#' @name beta_variants
NULL

#' @rdname beta_variants
#' @export
beta0 <- function(link = "logit") {
  structure(list(family = "beta0", link = link), class = "family")
}

#' @rdname beta_variants  
#' @export
beta1 <- function(link = "logit") {
  structure(list(family = "beta1", link = link), class = "family")
}

#' @rdname beta_variants
#' @export
logitbeta <- function(link = "logit") {
  structure(list(family = "logitbeta", link = link), class = "family")
}

# =============================================================================
# SECTION 2: ZERO-INFLATED AND HURDLE FAMILIES
# =============================================================================

#' Zero-Inflated Poisson Family
#'
#' @description
#' Zero-inflated Poisson for count data with excess zeros.
#'
#' @param link Link function for the mean (default: "log")
#' @param link.zi Link function for zero-inflation parameter (default: "logit") 
#'
#' @return A family object for use with qbrms()
#'
#' @examples
#' \dontrun{
#' # Count data with many zeros
#' fit <- qbrms(count ~ predictor, data = ecology_data, family = zero_inflated_poisson())
#' }
#'
#' @export
zero_inflated_poisson <- function(link = "log", link.zi = "logit") {
  structure(list(
    family = "zero_inflated_poisson", 
    link = link,
    link.zi = link.zi
  ), class = "family")
}
#' @rdname zero_inflated_poisson
#' @export  
zip <- function() {
  structure(list(family = "zero_inflated_poisson"), class = "family")
}

#' Zero-Inflated Negative Binomial Family
#'
#' @description
#' Zero-inflated negative binomial for count data with excess zeros and overdispersion.
#'
#' @param link Link function for the mean (default: "log")
#' @param link.zi Link function for zero-inflation parameter (default: "logit") 
#'
#' @return A family object for use with qbrms()
#'
#' @details
#' Combines a point mass at zero with a negative binomial distribution.
#' Useful for count data with both overdispersion and excess zeros beyond
#' what a negative binomial alone can handle.
#'
#' @examples
#' \dontrun{
#' # Count data with many zeros and overdispersion
#' fit <- qbrms(count ~ predictor, data = ecology_data, 
#'              family = zero_inflated_negbinomial())
#' 
#' # Model zero-inflation probability
#' fit <- qbrms(bf(count ~ habitat, zi ~ elevation), 
#'              data = ecology_data, family = zinb())
#' }
#'
#' @export
zero_inflated_negbinomial <- function(link = "log", link.zi = "logit") {
  structure(list(
    family = "zero_inflated_negbinomial", 
    link = link,
    link.zi = link.zi
  ), class = "family")
}

#' @rdname zero_inflated_negbinomial
#' @export
zinb <- function(link = "log", link.zi = "logit") {
  zero_inflated_negbinomial(link = link, link.zi = link.zi)
}

#' Hurdle Families for Two-Part Models
#'
#' @description
#' Hurdle models for count data, separating zero vs non-zero process.
#'
#' @param link Link function for count component (default: "log")
#' @param link.hu Link function for hurdle component (default: "logit")
#'
#' @return A family object for use with qbrms()
#'
#' @details
#' Two-part models: first part models zero vs non-zero probability,
#' second part models positive counts. Different from zero-inflation
#' as zeros only come from the hurdle process.
#'
#' @name hurdle_families
NULL

#' @rdname hurdle_families
#' @export
hurdle_poisson <- function(link = "log", link.hu = "logit") {
  structure(list(
    family = "hurdle_poisson",
    link = link,
    link.hu = link.hu
  ), class = "family")
}

#' @rdname hurdle_families
#' @export
hurdle_negbinomial <- function(link = "log", link.hu = "logit") {
  structure(list(
    family = "hurdle_negbinomial", 
    link = link,
    link.hu = link.hu
  ), class = "family")
}

# =============================================================================
# SECTION 3: SURVIVAL AND EXTREME VALUE FAMILIES
# =============================================================================

#' Weibull Survival Family
#'
#' @description
#' Weibull distribution for survival analysis and time-to-event data.
#'
#' @param link Link function for scale parameter (default: "log") 
#' @param link.shape Link function for shape parameter (default: "log")
#'
#' @return A family object for use with qbrms()
#'
#' @details
#' Standard Weibull parameterization for survival modeling. Can handle
#' both censored and uncensored data. Shape parameter determines whether
#' hazard increases (>1), decreases (<1), or remains constant (=1, exponential).
#'
#' @examples
#' \dontrun{
#' # Survival analysis
#' fit <- qbrms(time ~ treatment + age, data = survival_data, family = weibull())
#' 
#' # Model both scale and shape
#' fit <- qbrms(bf(time ~ treatment, shape ~ age), 
#'              data = survival_data, family = weibull())
#' }
#'
#' @export
weibull <- function(link = "log", link.shape = "log") {
  structure(list(
    family = "weibull",
    link = link,
    link.shape = link.shape
  ), class = "family")
}
#' Exponential Survival Family
#'
#' @description
#' Exponential distribution for survival analysis (special case of Weibull with shape=1).
#'
#' @param link Link function for rate parameter (default: "log")
#'
#' @return A family object for use with qbrms()
#'
#' @examples
#' \dontrun{
#' fit <- qbrms(time ~ treatment, data = survival_data, family = exponential())
#' fit2 <- qbrms(time ~ treatment, data = survival_data, family = exponential(link = "log"))
#' }
#'
#' @export
exponential <- function(link = "log") {
  structure(list(family = "exponential", link = link), class = "family")
}


#' Generalized Extreme Value Family
#'
#' @description  
#' Generalized Extreme Value (GEV) distribution for modeling extreme events.
#'
#' @param link Link function for location parameter (default: "identity")
#' @param link.sigma Link function for scale parameter (default: "log")
#' @param link.xi Link function for shape parameter (default: "identity") 
#'
#' @return A family object for use with qbrms()
#'
#' @details
#' Three-parameter family encompassing Gumbel (ξ=0), Fréchet (ξ>0), and 
#' Weibull (ξ<0) distributions. Useful for extreme value analysis such as
#' modeling maximum temperatures, flood levels, or financial extremes.
#'
#' @examples
#' \dontrun{
#' # Extreme precipitation modeling  
#' fit <- qbrms(max_precip ~ elevation + latitude, 
#'              data = climate_data, family = gev())
#' }
#'
#' @export
gev <- function(link = "identity", link.sigma = "log", link.xi = "identity") {
  structure(list(
    family = "gev",
    link = link,
    link.sigma = link.sigma,
    link.xi = link.xi
  ), class = "family")
}

#' @rdname gev
#' @export
gumbel <- function(link = "identity", link.sigma = "log") {
  structure(list(family = "gumbel", link = link, link.sigma = link.sigma), class = "family")
}

# =============================================================================
# SECTION 4: CIRCULAR AND DIRECTIONAL FAMILIES
# =============================================================================

#' Circular Normal Family for Directional Data
#'
#' @description
#' Circular normal (von Mises) distribution for modeling angular/directional data.
#'
#' @param link Link function for mean direction (default: "tan_half")
#' @param link.kappa Link function for concentration parameter (default: "log")
#'
#' @return A family object for use with qbrms()
#'
#' @details
#' Used for circular data such as wind directions, animal movement bearings,
#' or any periodic data on [0, 2π). The concentration parameter κ controls
#' how tightly distributed the data is around the mean direction.
#'
#' @examples
#' \dontrun{
#' # Wind direction analysis
#' fit <- qbrms(wind_direction ~ temperature + pressure, 
#'              data = weather_data, family = von_mises())
#' 
#' # Animal movement analysis  
#' fit <- qbrms(bearing ~ time_of_day, data = tracking_data, family = circular_normal())
#' }
#'
#' @export
circular_normal <- function(link = "tan_half", link.kappa = "log") {
  structure(list(
    family = "circular_normal",
    link = link,
    link.kappa = link.kappa
  ), class = "family")
}

#' @rdname circular_normal
#' @export
von_mises <- function(link = "tan_half", link.kappa = "log") {
  circular_normal(link = link, link.kappa = link.kappa)
}

# =============================================================================
# SECTION 5: ROBUST AND HEAVY-TAILED FAMILIES
# =============================================================================

#' Laplace (Double Exponential) Family
#'
#' @description
#' Laplace distribution for robust regression with heavier tails than normal.
#'
#' @param link Link function for location parameter (default: "identity")
#' @param link.sigma Link function for scale parameter (default: "log")
#'
#' @return A family object for use with qbrms()
#'
#' @details
#' Symmetric heavy-tailed distribution useful for robust regression when
#' outliers are expected. Equivalent to L1 regression (median regression)
#' in the maximum likelihood framework.
#'
#' @examples
#' \dontrun{
#' # Robust regression
#' fit <- qbrms(y ~ x, data = data_with_outliers, family = laplace())
#' }
#'
#' @export
laplace <- function(link = "identity", link.sigma = "log") {
  structure(list(family = "laplace", link = link, link.sigma = link.sigma), class = "family")
}

#' @rdname laplace
#' @export
double_exponential <- function(link = "identity", link.sigma = "log") {
  laplace(link = link, link.sigma = link.sigma)
}

#' Generalized t Family
#'
#' @description
#' Generalized t-distribution with flexible tail behavior.
#'
#' @param link Link function for location parameter (default: "identity")
#' @param link.sigma Link function for scale parameter (default: "log") 
#' @param link.nu Link function for degrees of freedom (default: "log")
#'
#' @return A family object for use with qbrms()
#'
#' @details
#' Extends Student's t with additional flexibility. As ν → ∞, approaches
#' normal distribution. Lower ν values give heavier tails for robustness.
#'
#' @export
gen_student_t <- function(link = "identity", link.sigma = "log", link.nu = "log") {
  structure(list(
    family = "gen_student_t",
    link = link,
    link.sigma = link.sigma,
    link.nu = link.nu
  ), class = "family")
}

# =============================================================================
# SECTION 6: SMOOTHING AND SPLINE-LIKE FAMILIES  
# =============================================================================

#' Random Walk Families for Smoothing
#'
#' @description
#' Random walk prior families for smooth function estimation.
#'
#' @param scale.model Scaling model for precision parameter
#' @param diagonal Diagonal precision matrix structure
#'
#' @return A family object for smooth terms
#'
#' @details
#' - rw1: First-order random walk (piecewise constant smoothing)
#' - rw2: Second-order random walk (piecewise linear smoothing)
#' 
#' Used internally for spline-like smooth function approximation in INLA.
#'
#' @name random_walk_families
NULL

#' @rdname random_walk_families
#' @export
rw1 <- function(scale.model = "log", diagonal = 1e-6) {
  structure(list(
    family = "rw1",
    scale.model = scale.model,
    diagonal = diagonal
  ), class = "family")
}

#' @rdname random_walk_families
#' @export
rw2 <- function(scale.model = "log", diagonal = 1e-6) {
  structure(list(
    family = "rw2", 
    scale.model = scale.model,
    diagonal = diagonal
  ), class = "family")
}

#' Independent Identically Distributed Random Effects
#'
#' @description
#' IID random effects with normal distribution.
#'
#' @param scale.model Scaling model for precision (default: "log")
#' @param diagonal Diagonal element for numerical stability (default: 1e-6)
#'
#' @return A family object for random effects
#'
#' @details
#' Standard IID normal random effects. Used for grouping factors and
#' exchangeable random effects structures.
#'
#' @examples
#' \dontrun{
#' # Used internally in mixed effects models
#' fit <- qbrms(y ~ x + (1|group), data = data, family = gaussian())
#' }
#'
#' @export
iid <- function(scale.model = "log", diagonal = 1e-6) {
  structure(list(
    family = "iid",
    scale.model = scale.model, 
    diagonal = diagonal
  ), class = "family")
}

#' Beta Family Constructor
#'
#' @description
#' Beta distribution for modeling proportional data in the interval (0,1).
#' Standard beta regression parameterization with mean and precision parameters.
#'
#' @param link Link function for the mean parameter (default: "logit")
#' @param link.phi Link function for the precision parameter (default: "log")
#'
#' @return A family object suitable for use with qbrms()
#'
#' @details
#' The beta distribution is parameterized with mean μ ∈ (0,1) and precision 
#' φ > 0. The variance is μ(1-μ)/(1+φ), so larger φ means smaller variance.
#' 
#' Response variable must be in interval (0,1). Values exactly 0 or 1 should 
#' be adjusted slightly before modeling.
#'
#' @examples
#' \dontrun{
#' # Basic beta regression
#' fit <- qbrms(proportion ~ treatment + age, data = data, family = Beta())
#' 
#' # Model both mean and precision
#' fit <- qbrms(bf(prop ~ predictor, phi ~ group), data = data, family = Beta())
#' 
#' # Alternative link functions
#' fit <- qbrms(prop ~ x, data = data, family = Beta("probit", "identity"))
#' }
#'
#' @export
Beta <- function(link = "logit", link.phi = "log") {
  
  # Validate link functions
  valid_links <- c("logit", "probit", "cloglog", "cauchit", "log", "identity")
  if (!link %in% valid_links) {
    stop("Invalid link function '", link, "'. Must be one of: ", 
         paste(valid_links, collapse = ", "))
  }
  
  valid_phi_links <- c("log", "identity", "inverse", "sqrt")
  if (!link.phi %in% valid_phi_links) {
    stop("Invalid link.phi function '", link.phi, "'. Must be one of: ", 
         paste(valid_phi_links, collapse = ", "))
  }
  
  # Create beta family structure
  structure(
    list(
      family = "beta",
      link = link,
      link.phi = link.phi,
      linkfun = make.link(link)$linkfun,
      linkinv = make.link(link)$linkinv,
      mu.eta = make.link(link)$mu.eta,
      validmu = function(mu) all(mu > 0 & mu < 1),
      valideta = function(eta) TRUE,
      initialize = expression({
        if (any(y <= 0 | y >= 1)) 
          stop("Beta regression requires response values in (0, 1)")
        n <- rep.int(1, nobs)
        mustart <- (y + 0.5/6) / (1 + 1/6)  # Adjustment for boundary values
      })
    ),
    class = "family"
  )
}


#' Log-Normal Family Constructor
#'
#' @description
#' Log-normal distribution for modeling positive continuous data.
#' Useful when the log of the response is approximately normally distributed.
#'
#' @param link Link function (default: "identity")
#'
#' @return A family object suitable for use with qbrms()
#'
#' @details
#' The log-normal distribution is appropriate for positive data that may be 
#' right-skewed. It's equivalent to assuming log(Y) ~ Normal(μ, σ²).
#' 
#' Response variable must be strictly positive (> 0).
#'
#' @examples
#' \dontrun{
#' # Basic log-normal regression
#' fit <- qbrms(income ~ education + experience, data = salary_data, 
#'              family = lognormal())
#' 
#' # Alternative link functions
#' fit <- qbrms(response ~ predictor, data = data, family = lognormal("log"))
#' }
#'
#' @export
lognormal <- function(link = "identity") {
  
  # Validate link functions
  valid_links <- c("identity", "log", "inverse")
  if (!link %in% valid_links) {
    stop("Invalid link function '", link, "'. Must be one of: ", 
         paste(valid_links, collapse = ", "))
  }
  
  # Create lognormal family structure
  structure(
    list(
      family = "lognormal",
      link = link,
      linkfun = make.link(link)$linkfun,
      linkinv = make.link(link)$linkinv,
      mu.eta = make.link(link)$mu.eta,
      validmu = function(mu) all(mu > 0),
      valideta = function(eta) TRUE,
      initialize = expression({
        if (any(y <= 0)) 
          stop("Log-normal regression requires positive response values")
        n <- rep.int(1, nobs)
        mustart <- y
      })
    ),
    class = "family"
  )
}

# =============================================================================
# SECTION 7: FAMILY INFORMATION AND UTILITIES
# =============================================================================

#' Get Family Documentation 
#'
#' @description
#' Returns documentation and usage information for a family.
#'
#' @param family_name Name of the family
#' @return Character string with family information
#'
#' @export
family_info <- function(family_name) {
  
  info_map <- list(
    "simplex" = "Simplex family for proportional data in (0,1) with flexible precision modeling",
    "beta_binomial" = "Beta-binomial for overdispersed binary data",
    "zero_inflated_poisson" = "Zero-inflated Poisson for count data with excess zeros", 
    "zero_inflated_negbinomial" = "Zero-inflated negative binomial for overdispersed counts with excess zeros",
    "weibull" = "Weibull distribution for survival analysis and reliability modeling",
    "exponential" = "Exponential distribution for survival analysis (constant hazard)",
    "gev" = "Generalized extreme value distribution for modeling extreme events",
    "circular_normal" = "Circular normal (von Mises) for directional/angular data",
    "von_mises" = "von Mises (circular normal) for directional/angular data",
    "laplace" = "Laplace (double exponential) for robust regression with heavy tails",
    "gen_student_t" = "Generalized t-distribution with flexible tail behavior",
    "hurdle_poisson" = "Hurdle Poisson for two-part count models",
    "hurdle_negbinomial" = "Hurdle negative binomial for two-part overdispersed count models"
  )
  
  if (family_name %in% names(info_map)) {
    return(info_map[[family_name]])
  } else {
    return(paste("Family", family_name, "- see documentation for details"))
  }
}

#' List Available Extended Families
#'
#' @description
#' Returns a data frame of all available extended families with categories.
#'
#' @return Data frame with family names, categories, and brief descriptions
#' @export
list_extended_families <- function() {
  
  families_df <- data.frame(
    Family = c(
      # Compositional
      "simplex", "beta_binomial", "beta0", "beta1", "logitbeta",
      # Zero-inflated
      "zero_inflated_poisson", "zip", "zero_inflated_negbinomial", "zinb",
      # Hurdle
      "hurdle_poisson", "hurdle_negbinomial",
      # Survival
      "weibull", "exponential", "gev", "gumbel",
      # Circular
      "circular_normal", "von_mises",
      # Robust
      "laplace", "double_exponential", "gen_student_t",
      # Smoothing
      "rw1", "rw2", "iid"
    ),
    Category = c(
      # Compositional
      rep("Compositional", 5),
      # Zero-inflated
      rep("Zero-inflated", 4),
      # Hurdle
      rep("Hurdle", 2),
      # Survival
      rep("Survival", 4),
      # Circular
      rep("Circular", 2),
      # Robust
      rep("Robust", 3),
      # Smoothing
      rep("Smoothing", 3)
    ),
    Description = c(
      # Compositional
      "Proportional data (0,1)", "Overdispersed binomial", "Beta variant", "Beta variant", "Logit-beta",
      # Zero-inflated
      "Poisson with excess zeros", "ZIP alias", "NegBin with excess zeros", "ZINB alias",
      # Hurdle
      "Two-part Poisson", "Two-part negative binomial",
      # Survival
      "Weibull survival", "Exponential survival", "Extreme values", "Gumbel (GEV)",
      # Circular
      "Directional data", "von Mises alias",
      # Robust
      "Robust regression", "Laplace alias", "Heavy-tailed",
      # Smoothing
      "Random walk 1", "Random walk 2", "IID random effects"
    ),
    stringsAsFactors = FALSE
  )
  
  return(families_df)
}
