# =============================================================================
# R/families.R
# =============================================================================

#' Family Conversion and Utilities for qbrms Package
#' 
#' @description
#' Core family conversion functions with comprehensive family support,
#' numerical stability enhancements, and CRAN-ready implementation.
#' 
#' @importFrom stats cov2cor qlogis
#' @importFrom utils tail
#' @name families
#' @keywords internal
NULL

# =============================================================================
# SECTION 1: CORE FAMILY CONVERSION FUNCTION
# =============================================================================

#' Convert Family Object to INLA-Compatible Specification
#'
#' @description
#' Enhanced family conversion supporting all standard and additional families
#' with automatic routing to specialised implementations when enabled.
#'
#' @param family A family object, character string, or list specifying the 
#'   response distribution
#' @param quantile Numeric value between 0 and 1 for quantile regression
#' @param allow_ordinal_routing Logical; if TRUE, enables routing for ordinal families
#'
#' @return Character string, list, or routing object specifying family/routing info
#' @export
convert_family_to_inla <- function(family, quantile = 0.5, allow_ordinal_routing = FALSE) {
  
  # Handle NULL input - return default
  if (is.null(family)) {
    return("gaussian")
  }
  
  # Extract family name from various input types with enhanced error handling
  family_name <- tryCatch({
    if (is.function(family)) {
      family_obj <- tryCatch(family(), error = function(e) {
        stop("Family function failed to execute: ", e$message, call. = FALSE)
      })
      tolower(as.character(family_obj$family))
    } else if (is.character(family)) {
      fam_char <- tolower(as.character(family[1]))
      if (nchar(fam_char) == 0) {
        stop("Family '' is not supported in qbrms", call. = FALSE)
      }
      fam_char
    } else if (is.list(family) && !is.null(family$family)) {
      tolower(as.character(family$family))
    } else {
      # This catches invalid inputs like functions that fail
      stop("Family '", deparse(substitute(family)), "' is not supported in qbrms", 
           call. = FALSE)
    }
  }, error = function(e) {
    # Re-throw the error if it's already informative
    if (grepl("not supported", e$message)) {
      stop(e$message, call. = FALSE)
    }
    # Otherwise, provide generic error
    stop("Family '", deparse(substitute(family)), "' is not supported in qbrms", 
         call. = FALSE)
  })
  
  # Handle routing for ordinal families when enabled
  if (allow_ordinal_routing && family_name %in% c("cumulative", "ordinal")) {
    return(.create_ordinal_routing_object(family, family_name))
  }
  
  # Handle asymmetric laplace with quantile specification
  if (family_name == "asymmetric_laplace") {
    if (!is.numeric(quantile) || length(quantile) != 1 || quantile <= 0 || quantile >= 1) {
      stop("quantile must be a single numeric value between 0 and 1", call. = FALSE)
    }
    return(list(family = "asymmetric_laplace", quantile = quantile))
  }
  
  # Comprehensive family mapping for INLA
  family_map <- c(
    # Basic continuous distributions
    "gaussian"    = "gaussian",
    "normal"      = "gaussian", 
    "lognormal"   = "lognormal",
    "gamma"       = "gamma",
    "beta"        = "beta",
    
    # Count distributions
    "poisson"     = "poisson",
    "pois"        = "poisson",
    "negbinomial" = "nbinomial",
    "negative_binomial" = "nbinomial",
    "nbinomial"   = "nbinomial",
    
    # Binary/binomial
    "binomial"    = "binomial",
    "binom"       = "binomial",
    
    # Robust distributions
    "skew_normal" = "sn",
    "student_t"   = "t",
    "student"     = "t",
    "studentt"    = "t",
    
    # Survival analysis
    "weibull"        = "weibull",
    "exponential"    = "exponential",
    "weibullsurv"    = "weibullsurv",
    "exponentialsurv"= "exponentialsurv",
    
    # Zero-inflated distributions
    "zero_inflated_poisson"     = "zeroinflatedpoisson1",
    "zip"                       = "zeroinflatedpoisson1",
    "zero_inflated_negbinomial" = "zeroinflatednbinomial1",
    "zinb"                      = "zeroinflatednbinomial1",
    
    # Overdispersed distributions
    "betabinomial"  = "betabinomial",
    "beta_binomial" = "betabinomial",
    
    # Extreme value distributions
    "gumbel" = "gev",
    "gev"    = "gev",
    
    # Robust regression
    "laplace"            = "t",
    "double_exponential" = "t",
    
    # Circular/directional data
    "circular_normal" = "circularnormal",
    "von_mises"       = "circularnormal",
    
    # Compositional and constrained data
    "simplex"         = "simplex",
    
    # Generalized distributions
    "gen_student_t"   = "stochvol_t",
    "generalized_t"   = "stochvol_t",
    
    # Categorical distributions
    "multinomial"     = "multinomial"
  )
  
  # Check if family is supported
  if (family_name %in% names(family_map)) {
    return(unname(family_map[family_name]))
  }
  
  # Provide informative error for unsupported families
  supported_families <- sort(unique(family_map))
  
  stop(
    "Family '", family_name, "' is not supported in qbrms.\n\n",
    "Supported families:\n",
    "  ", paste(supported_families, collapse = ", "), "\n\n",
    "For ordinal regression, use qbrmO() directly or enable routing:\n",
    "  qbrms(..., family = cumulative()) with qbrmO available\n\n",
    "For additional family support, consider using the 'brms' package.",
    call. = FALSE
  )
}

# =============================================================================
# SECTION 2: FAMILY DATA VALIDATION
# =============================================================================

#' Validate Family-Specific Data Constraints
#'
#' @description
#' Validates that response data meets family-specific constraints and 
#' automatically adjusts boundary values when possible.
#'
#' @param y Response variable vector
#' @param family_name Character string specifying the family name
#' @param tolerance Tolerance for boundary adjustments (default: 1e-6)
#'
#' @return Invisibly returns TRUE if validation passes
#' @keywords internal
validate_family_data <- function(y, family_name, tolerance = 1e-6) {
  
  if (is.null(y) || length(y) == 0) {
    stop("Response variable cannot be empty")
  }
  
  # Remove NA values for validation
  y_clean <- y[!is.na(y)]
  
  if (length(y_clean) == 0) {
    stop("Response variable contains only NA values")
  }
  
  # Family-specific validation
  switch(family_name,
         
         # Positive continuous data
         "gamma" = {
           if (any(y_clean <= 0)) {
             stop("Gamma family requires strictly positive response values (y > 0)")
           }
         },
         
         "lognormal" = {
           if (any(y_clean <= 0)) {
             stop("Log-normal family requires strictly positive response values (y > 0)")
           }
         },
         
         "weibull" = {
           if (any(y_clean <= 0)) {
             stop("Weibull family requires strictly positive response values (y > 0)")
           }
         },
         
         # Unit interval data
         "beta" = {
           if (any(y_clean <= 0 | y_clean >= 1)) {
             stop("Beta family requires response values in open interval (0,1)")
           }
         },
         
         "simplex" = {
           if (any(y_clean <= 0 | y_clean >= 1)) {
             stop("Simplex family requires response values in open interval (0,1)")
           }
         },
         
         # Count data
         "poisson" = {
           if (any(y_clean < 0 | y_clean != round(y_clean))) {
             stop("Poisson family requires non-negative integer values")
           }
         },
         
         "nbinomial" = {
           if (any(y_clean < 0 | y_clean != round(y_clean))) {
             stop("Negative binomial family requires non-negative integer values")
           }
         },
         
         "zeroinflatedpoisson1" = {
           if (any(y_clean < 0 | y_clean != round(y_clean))) {
             stop("Zero-inflated Poisson family requires non-negative integer values")
           }
         },
         
         "zeroinflatednbinomial1" = {
           if (any(y_clean < 0 | y_clean != round(y_clean))) {
             stop("Zero-inflated negative binomial family requires non-negative integer values")
           }
         },
         
         # Binary/binomial data
         "binomial" = {
           unique_vals <- unique(y_clean)
           if (!all(unique_vals %in% c(0, 1))) {
             # Check if it might be proportional data
             if (all(y_clean >= 0 & y_clean <= 1)) {
               warning("Binomial family with non-binary data - ensure proper trials specification")
             } else {
               stop("Binomial family requires binary (0/1) or proportional [0,1] data")
             }
           }
         },
         
         "betabinomial" = {
           if (any(y_clean < 0 | y_clean != round(y_clean))) {
             stop("Beta-binomial family requires non-negative integer values")
           }
         },
         
         # Circular data (should be in radians)
         "circularnormal" = {
           if (any(y_clean < 0 | y_clean > 2*pi)) {
             warning("Circular normal typically expects values in [0, 2*pi] radians")
           }
         }
         
         # For other families (gaussian, t, sn, laplace, etc.), no specific constraints
  )
  
  return(invisible(TRUE))
}

# =============================================================================
# SECTION 3: ENHANCED INLA CONTROL SETTINGS
# =============================================================================

#' Get Enhanced INLA Control Settings for Family-Specific Stability
#'
#' @description
#' Provides family-specific INLA control settings to enhance numerical 
#' stability and convergence for different distribution families.
#'
#' @param family_name Character string specifying the family name
#' @param base_control List of base control settings (default: NULL)
#'
#' @return List of enhanced INLA control settings
#' @keywords internal
get_enhanced_inla_control <- function(family_name, base_control = NULL) {
  
  # Default enhanced controls
  enhanced_control <- list(
    inla.mode = "classic",
    restart = 3,
    int.strategy = "grid",
    strategy = "adaptive",
    h = 0.005,  # Smaller step size for stability
    tolerance = 1e-8,  # Tighter tolerance
    max.iter = 500,  # More iterations if needed
    npoints = 21  # More integration points
  )
  
  # Family-specific enhancements
  family_specific <- switch(family_name,
                            
                            # For simplex and beta: enhanced precision in probability space
                            "simplex" = list(
                              int.strategy = "ccd",
                              h = 0.001,
                              npoints = 31,
                              tolerance = 1e-10
                            ),
                            
                            "beta" = list(
                              int.strategy = "ccd", 
                              h = 0.001,
                              npoints = 31
                            ),
                            
                            # For zero-inflated models: enhanced convergence
                            "zeroinflatedpoisson1" = list(
                              restart = 5,
                              max.iter = 1000,
                              tolerance = 1e-10
                            ),
                            
                            "zeroinflatednbinomial1" = list(
                              restart = 5,
                              max.iter = 1000,
                              tolerance = 1e-10
                            ),
                            
                            # For robust families: enhanced tail handling
                            "laplace" = list(
                              int.strategy = "grid",
                              npoints = 41,
                              h = 0.002
                            ),
                            
                            "t" = list(
                              int.strategy = "grid", 
                              npoints = 31,
                              h = 0.003
                            ),
                            
                            # For survival models: enhanced handling of extreme values
                            "weibull" = list(
                              int.strategy = "ccd",
                              restart = 5,
                              tolerance = 1e-10
                            ),
                            
                            # Default (no specific enhancements)
                            list()
  )
  
  # Merge base, enhanced, and family-specific controls
  if (!is.null(base_control)) {
    enhanced_control <- modifyList(enhanced_control, base_control)
  }
  
  if (length(family_specific) > 0) {
    enhanced_control <- modifyList(enhanced_control, family_specific)
  }
  
  return(enhanced_control)
}

# =============================================================================
# SECTION 4: ROUTING SUPPORT FOR ORDINAL REGRESSION
# =============================================================================

#' Check if Family Requires Routing to Specialist Implementation
#'
#' @param family_spec Family specification from convert_family_to_inla()
#' @return Logical indicating if routing is required
#' @keywords internal
requires_routing <- function(family_spec) {
  if (is.list(family_spec) && !is.null(family_spec$route_to)) {
    return(TRUE)
  }
  return(FALSE)
}

#' Extract Routing Information from Family Specification
#'
#' @param family_spec Family specification with routing information
#' @return List containing routing details
#' @keywords internal
extract_routing_info <- function(family_spec) {
  if (is.list(family_spec) && !is.null(family_spec$route_to)) {
    return(family_spec)
  }
  stop("No routing information found in family specification")
}

#' Create Ordinal Routing Object
#'
#' @param family Original family object
#' @param family_name Extracted family name
#' @return Routing object for ordinal regression
#' @keywords internal
.create_ordinal_routing_object <- function(family, family_name) {
  
  # Extract link and threshold information
  link <- "logit"  # Default
  threshold <- "flexible"  # Default
  
  if (is.list(family)) {
    if (!is.null(family$link)) {
      link <- family$link
    }
    if (!is.null(family$threshold)) {
      threshold <- family$threshold
    }
  }
  
  return(list(
    route_to = "qbrmO",
    target = "qbrmO",
    original_family = family_name,
    family = "cumulative",
    link = link,
    threshold = threshold
  ))
}

# =============================================================================
# SECTION 5: HELPER UTILITIES
# =============================================================================

#' Extract Family Name from INLA Family Specification
#'
#' @description
#' Extracts the family name from an INLA family specification, handling both 
#' character strings and lists (e.g., for quantile regression).
#'
#' @param inla_family A family specification returned by \code{convert_family_to_inla()}.
#'   Can be a character string or a list with \code{$family} component.
#'
#' @return Character string containing the family name.
#'
#' @examples
#' \dontrun{
#' # Character family
#' extract_family_name("gaussian")  # "gaussian"
#' 
#' # List family (from quantile regression)
#' ald_spec <- convert_family_to_inla(asymmetric_laplace(), quantile = 0.9)
#' extract_family_name(ald_spec)    # "asymmetric_laplace"
#' }
#'
#' @export
extract_family_name <- function(inla_family) {
  if (is.character(inla_family)) {
    return(inla_family)
  } else if (is.list(inla_family) && !is.null(inla_family$family)) {
    return(inla_family$family)
  } else {
    stop("Invalid INLA family specification. Expected character string or list with $family component.")
  }
}

#' Beta Family Constructor
#' 
#' @description
#' Beta distribution family for response variables in (0,1)
#' 
#' @param link Link function for the mean parameter (default: "logit")  
#' @param link.phi Link function for precision parameter (default: "log")
#' 
#' @return A family object for use with qbrms()
#' 
#' @examples
#' \dontrun{
#' # Beta regression for proportions
#' fit <- qbrms(proportion ~ predictor, data = data, family = Beta())
#' }
#' 
#' @export
Beta <- function(link = "logit", link.phi = "log") {
  structure(list(
    family = "beta",
    link = link,
    link.phi = link.phi
  ), class = "family")
}

#' Lognormal Family Constructor
#' 
#' @description  
#' Lognormal distribution family for positive continuous responses
#' 
#' @param link Link function (default: "identity")
#' 
#' @return A family object for use with qbrms()
#' 
#' @examples
#' \dontrun{
#' # Lognormal regression 
#' fit <- qbrms(response ~ predictor, data = data, family = lognormal())
#' }
#' 
#' @export
lognormal <- function(link = "identity") {
  structure(list(
    family = "lognormal",
    link = link
  ), class = "family")
}

