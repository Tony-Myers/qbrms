# =============================================================================
# R/verbosity.R
# =============================================================================

#' Set qbrms verbosity for the current session
#'
#' @param verbose Logical. If TRUE, fitting prints progress; if FALSE, fitting is silent.
#' @return Invisibly returns the previous value.
#' @export
qbrms_set_verbosity <- function(verbose = FALSE) {
  old <- getOption("qbrms.verbose", FALSE)
  options(qbrms.verbose = isTRUE(verbose))
  invisible(old)
}

#' @keywords internal
.qbrms_capture_all <- function(expr) {
  if (isTRUE(getOption("qbrms.verbose", FALSE))) {
    return(eval.parent(substitute(expr)))
  }
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  zz <- file(tmp, open = "wt")
  on.exit({ try(sink(type = "message")); try(sink()); close(zz) }, add = TRUE)
  sink(zz); sink(zz, type = "message")
  res <- eval.parent(substitute(expr))
  res
}

# Silence TMB build chatter when verbose = FALSE
#' @keywords internal
.qbrms_silence_tmb <- function(expr) {
  if (isTRUE(getOption("qbrms.verbose", FALSE))) {
    return(eval.parent(substitute(expr)))
  }
  old_make <- Sys.getenv("MAKEFLAGS", unset = NA)
  on.exit({
    if (is.na(old_make)) Sys.unsetenv("MAKEFLAGS") else Sys.setenv(MAKEFLAGS = old_make)
  }, add = TRUE)
  # -s tells make to be silent
  Sys.setenv(MAKEFLAGS = "-s")
  .qbrms_capture_all(expr)
}

#' Get captured fit log from a qbrms object (if available)
#'
#' @param x A qbrms_fit / qbrmO_fit object returned by qbrm()/qbrms()
#' @return A character vector of captured console lines, or NULL if none
#' @export
qbrms_fit_log <- function(x) {
  attr(x, "qbrms.captured_output", exact = TRUE)
}

# Internal: run an expression silently unless verbosity is enabled.
# - does not swallow errors;
# - captures stdout/messages/warnings so users can inspect via qbrms_fit_log();
# - returns the expression's value unchanged, with an attribute storing the log.
.qbrms_silently <- function(expr) {
  verbose <- isTRUE(getOption("qbrms.verbose", FALSE))
  if (verbose) return(force(expr))
  
  out <- character(0)
  msgs <- character(0)
  warns <- character(0)
  
  res <- withCallingHandlers(
    tryCatch(
      {
        out <- utils::capture.output(
          value <- withRestarts(
            suppressWarnings(
              suppressPackageStartupMessages(
                suppressMessages(
                  force(expr)
                )
              )
            ),
            muffleMessage = function() {},
            muffleWarning = function() {}
          ),
          type = "output"
        )
        value
      },
      error = function(e) {
        # rethrow with captured output attached to aid debugging if needed
        attr(e, "qbrms.captured_output") <- out
        stop(e)
      }
    ),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  
  # Attach a compact log for later inspection
  log <- c(out, if (length(msgs)) paste0("[MESSAGE] ", msgs), if (length(warns)) paste0("[WARNING] ", warns))
  if (!is.null(res) && (is.list(res) || is.environment(res) || isS4(res))) {
    # store on the returned object when possible
    try(attr(res, "qbrms.captured_output") <- log, silent = TRUE)
  }
  res
}