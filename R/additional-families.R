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
#' @param link Link function for the mean (default: "logit")
#' @param link.precision Link function for precision parameter (default: "log")
#' @return A family object
#' @export
simplex <- function(link = "logit", link.precision = "log") {
  structure(list(
    family = "simplex",
    link = link,
    link.precision = link.precision
  ), class = "family")
}

#' Beta Binomial Family for Overdispersed Binary Data
#' @param link Link function for probability parameter (default: "logit")
#' @return A family object
#' @export
beta_binomial <- function(link = "logit") {
  structure(list(
    family = "beta_binomial",
    link = link
  ), class = "family")
}

#' Alternative Beta Parameterizations
#' @param link Link function (default: "logit")
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
#' @param link Link function for mean (default: "log")
#' @param link.zi Link function for zero-inflation (default: "logit")
#' @return A family object
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
#' @param link Link function for mean (default: "log")
#' @param link.zi Link function for zero-inflation (default: "logit")
#' @return A family object
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
#' @param link Link function for count component
#' @param link.hu Link function for hurdle component
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

#' Gamma family (GLM-style)
#' @description Gamma family constructor to avoid conflict with base::gamma.
#' @param link Link function (default: "log")
#' @return A family object
#' @name Gamma_family
#' @rdname Gamma_family
#' @aliases Gamma
#' @export
Gamma <- function(link = "log") {
  structure(list(family = "gamma", link = link), class = "family")
}

# =============================================================================
# SECTION 3: SURVIVAL AND EXTREME VALUE FAMILIES
# =============================================================================

#' Weibull Survival Family
#' @param link Link function for scale (default: "log")
#' @param link.shape Link function for shape (default: "log")
#' @return A family object
#' @export
weibull <- function(link = "log", link.shape = "log") {
  structure(list(
    family = "weibull",
    link = link,
    link.shape = link.shape
  ), class = "family")
}

#' Generalized Extreme Value Family
#' @param link Link function for location
#' @param link.sigma Link function for scale
#' @param link.xi Link function for shape
#' @return A family object
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
#' @param link Link function for mean direction
#' @param link.kappa Link function for concentration
#' @return A family object
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
#' @param link Link function for location
#' @param link.sigma Link function for scale
#' @return A family object
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
#' @param link Link function for location
#' @param link.sigma Link function for scale
#' @param link.nu Link function for degrees of freedom
#' @return A family object
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

#' Random Walk Families
#' @param scale.model Scaling model for precision
#' @param diagonal Diagonal precision matrix structure
#' @return A family object
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

#' IID Random Effects
#' @param scale.model Scaling model for precision
#' @param diagonal Diagonal precision matrix structure
#' @return A family object
#' @export
iid <- function(scale.model = "log", diagonal = 1e-6) {
  structure(list(
    family = "iid",
    scale.model = scale.model, 
    diagonal = diagonal
  ), class = "family")
}

#' Beta Family Constructor (Capital B)
#' @param link Link function (default: "logit")
#' @param link.phi Link function for precision (default: "log")
#' @return A family object
#' @export
Beta <- function(link = "logit", link.phi = "log") {
  valid_links <- c("logit", "probit", "cloglog", "cauchit", "log", "identity")
  if (!link %in% valid_links) stop("Invalid link function")
  
  structure(list(
    family = "beta",
    link = link,
    link.phi = link.phi
  ), class = "family")
}

# =============================================================================
# SECTION 7: FAMILY INFORMATION AND UTILITIES
# =============================================================================

#' Get Family Documentation 
#' @param family_name Name of the family
#' @return Character string with family information
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
#' @return Data frame with family names, categories, and brief descriptions
#' @export
list_extended_families <- function() {
  data.frame(
    Family = c(
      "simplex", "beta_binomial", "beta0", "beta1", "logitbeta",
      "zero_inflated_poisson", "zip", "zero_inflated_negbinomial", "zinb",
      "hurdle_poisson", "hurdle_negbinomial",
      "weibull", "exponential", "gev", "gumbel",
      "circular_normal", "von_mises",
      "laplace", "double_exponential", "gen_student_t",
      "rw1", "rw2", "iid"
    ),
    Category = c(
      rep("Compositional", 5),
      rep("Zero-inflated", 4),
      rep("Hurdle", 2),
      rep("Survival", 4),
      rep("Circular", 2),
      rep("Robust", 3),
      rep("Smoothing", 3)
    ),
    Description = c(
      "Proportional data (0,1)", "Overdispersed binomial", "Beta variant", "Beta variant", "Logit-beta",
      "Poisson with excess zeros", "ZIP alias", "NegBin with excess zeros", "ZINB alias",
      "Two-part Poisson", "Two-part negative binomial",
      "Weibull survival", "Exponential survival", "Extreme values", "Gumbel (GEV)",
      "Directional data", "von Mises alias",
      "Robust regression", "Laplace alias", "Heavy-tailed",
      "Random walk 1", "Random walk 2", "IID random effects"
    ),
    stringsAsFactors = FALSE
  )
}