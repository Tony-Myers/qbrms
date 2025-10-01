# R/zzz.R - Package startup and cleanup

.onAttach <- function(libname, pkgname) {
  # Check for INLA availability
  if (!requireNamespace("INLA", quietly = TRUE)) {
    packageStartupMessage(
      "qbrms: INLA package not found. Some functionality will use fallback methods.\n",
      "To install INLA: install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)"
    )
  }
  
  # Check for other important suggested packages
  missing_packages <- character(0)
  suggested_packages <- c("lme4", "quantreg", "ordinal", "posterior")
  
  for (pkg in suggested_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0 && length(missing_packages) < length(suggested_packages)) {
    packageStartupMessage(
      "Note: Some optional packages not available: ", 
      paste(missing_packages, collapse = ", "), 
      "\nSome functionality will use fallback methods."
    )
  }
}