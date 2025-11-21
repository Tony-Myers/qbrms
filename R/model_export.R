# =============================================================================
# R/model_export.R
# =============================================================================

#' Export Model Specification
#'
#' @description
#' Export model specifications to various formats for sharing, documentation,
#' or reproduction.
#'
#' @param model A fitted qbrms model object or qbrms_model_spec object
#' @param file Character string specifying output file path
#' @param format Character string specifying export format: "R" (R script),
#'   "markdown" (Rmd document), "text" (plain text), or "json" (JSON format)
#' @param include_data Logical; if TRUE, includes data summary in export
#'   (default: TRUE)
#' @param include_diagnostics Logical; if TRUE and model is fitted, includes
#'   diagnostic information (default: FALSE)
#'
#' @return Invisibly returns the export content as a character string
#'
#' @details
#' This function facilitates model sharing and documentation by exporting:
#' \itemize{
#'   \item Model formula and family specification
#'   \item Prior specifications (if any)
#'   \item Data summary and structure
#'   \item Model fitting code
#'   \item Results summary (for fitted models)
#'   \item Diagnostic information (if requested)
#' }
#'
#' The exported content can be used to:
#' \itemize{
#'   \item Share analyses with collaborators
#'   \item Document modelling decisions
#'   \item Create reproducible research reports
#'   \item Archive model specifications
#' }
#'
#' @examples
#' \dontrun{
#' # Export model specification
#' spec <- model_builder(data = mtcars, response = "mpg")
#' export_model(spec, "my_model_spec.R", format = "R")
#'
#' # Export fitted model
#' fit <- qbrms(mpg ~ hp + wt, data = mtcars, family = gaussian())
#' export_model(fit, "my_model.Rmd", format = "markdown", 
#'              include_diagnostics = TRUE)
#'
#' # Export as JSON
#' export_model(spec, "my_model.json", format = "json")
#' }
#'
#' @export
export_model <- function(model, file, format = c("R", "markdown", "text", "json"),
                         include_data = TRUE, include_diagnostics = FALSE) {
  
  format <- match.arg(format)
  
  # Determine model type
  is_fitted <- inherits(model, "qbrms_fit") || inherits(model, "qbrmO_fit")
  is_spec <- inherits(model, "qbrms_model_spec")
  
  if (!is_fitted && !is_spec) {
    stop("model must be a fitted qbrms model or qbrms_model_spec object",
         call. = FALSE)
  }
  
  # Generate content based on format
  content <- switch(format,
                    "R" = .export_as_r_script(model, include_data, include_diagnostics, is_fitted),
                    "markdown" = .export_as_markdown(model, include_data, include_diagnostics, is_fitted),
                    "text" = .export_as_text(model, include_data, include_diagnostics, is_fitted),
                    "json" = .export_as_json(model, include_data, is_fitted),
                    stop("Unsupported format", call. = FALSE)
  )
  
  # Write to file
  writeLines(content, file)
  
  message("Model exported to: ", file)
  invisible(content)
}

#' Export as R Script
#' @keywords internal
.export_as_r_script <- function(model, include_data, include_diagnostics, is_fitted) {
  
  lines <- character(0)
  
  # Header
  lines <- c(lines,
             "# =============================================================================",
             "# qbrms Model Specification",
             sprintf("# Generated: %s", Sys.time()),
             "# =============================================================================",
             "",
             "# Load required packages",
             "library(qbrms)",
             ""
  )
  
  # Extract model components
  if (is_fitted) {
    formula <- model$formula
    family_name <- model$family
    data_name <- "model_data"  # Placeholder
  } else {
    formula <- model$formula
    family_name <- model$family$name
    data_name <- deparse(substitute(model$data))
  }
  
  # Data summary
  if (include_data) {
    lines <- c(lines,
               "# Data Summary",
               "# ------------",
               if (is_fitted) {
                 c(sprintf("# Data: %d observations, %d variables", 
                           nrow(model$data), ncol(model$data)),
                   sprintf("# Response: %s", as.character(formula[[2]])))
               } else {
                 c(sprintf("# Data: %d observations, %d variables",
                           nrow(model$data), ncol(model$data)),
                   sprintf("# Response: %s", model$response_info$name),
                   sprintf("# Response type: %s", model$response_info$type))
               },
               ""
    )
  }
  
  # Model specification
  lines <- c(lines,
             "# Model Specification",
             "# -------------------",
             sprintf("# Formula: %s", deparse(formula)),
             sprintf("# Family: %s", family_name),
             ""
  )
  
  # Model code
  if (is_fitted) {
    model_code <- sprintf(
      "fit <- qbrms(\n  formula = %s,\n  data = %s,\n  family = %s()\n)",
      deparse(formula),
      data_name,
      family_name
    )
  } else {
    model_code <- model$model_code
  }
  
  lines <- c(lines,
             "# Model Fitting Code",
             "# ------------------",
             model_code,
             ""
  )
  
  # Results (if fitted)
  if (is_fitted) {
    lines <- c(lines,
               "# View Results",
               "# ------------",
               "summary(fit)",
               "plot(fit)",
               ""
    )
  }
  
  # Diagnostics (if requested and fitted)
  if (include_diagnostics && is_fitted) {
    lines <- c(lines,
               "# Model Diagnostics",
               "# -----------------",
               "diagnostics <- diagnose_model(fit)",
               "print(diagnostics)",
               "plot(diagnostics)",
               ""
    )
  }
  
  return(lines)
}

#' Export as Markdown
#' @keywords internal
.export_as_markdown <- function(model, include_data, include_diagnostics, is_fitted) {
  
  lines <- character(0)
  
  # YAML header
  lines <- c(lines,
             "---",
             "title: 'qbrms Model Specification'",
             sprintf("date: '%s'", Sys.Date()),
             "output: html_document",
             "---",
             "",
             "```{r setup, include=FALSE}",
             "knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)",
             "library(qbrms)",
             "```",
             ""
  )
  
  # Introduction
  lines <- c(lines,
             "## Model Overview",
             "",
             "This document contains the specification and results for a Bayesian regression model",
             "fitted using the qbrms package.",
             ""
  )
  
  # Extract model components
  if (is_fitted) {
    formula <- model$formula
    family_name <- model$family
  } else {
    formula <- model$formula
    family_name <- model$family$name
  }
  
  # Data section
  if (include_data) {
    lines <- c(lines,
               "## Data Description",
               ""
    )
    
    if (is_fitted) {
      lines <- c(lines,
                 sprintf("- **Observations**: %d", nrow(model$data)),
                 sprintf("- **Variables**: %d", ncol(model$data)),
                 sprintf("- **Response variable**: %s", as.character(formula[[2]])),
                 ""
      )
    } else {
      lines <- c(lines,
                 sprintf("- **Observations**: %d", nrow(model$data)),
                 sprintf("- **Variables**: %d", ncol(model$data)),
                 sprintf("- **Response variable**: %s", model$response_info$name),
                 sprintf("- **Response type**: %s", model$response_info$type),
                 ""
      )
    }
  }
  
  # Model specification
  lines <- c(lines,
             "## Model Specification",
             "",
             sprintf("**Formula**: `%s`", deparse(formula)),
             "",
             sprintf("**Family**: %s", family_name),
             ""
  )
  
  # Model code
  lines <- c(lines,
             "## Model Code",
             "",
             "```{r model_code, eval=FALSE}",
             if (is_fitted) {
               sprintf("fit <- qbrms(\n  formula = %s,\n  data = model_data,\n  family = %s()\n)",
                       deparse(formula), family_name)
             } else {
               model$model_code
             },
             "```",
             ""
  )
  
  # Results (if fitted)
  if (is_fitted) {
    lines <- c(lines,
               "## Model Results",
               "",
               "```{r results}",
               "summary(fit)",
               "```",
               "",
               "### Model Plots",
               "",
               "```{r plots, fig.width=8, fig.height=6}",
               "plot(fit)",
               "```",
               ""
    )
  }
  
  # Diagnostics
  if (include_diagnostics && is_fitted) {
    lines <- c(lines,
               "## Model Diagnostics",
               "",
               "```{r diagnostics}",
               "diagnostics <- diagnose_model(fit)",
               "print(diagnostics)",
               "```",
               "",
               "```{r diagnostic_plots, fig.width=8, fig.height=6}",
               "plot(diagnostics)",
               "```",
               ""
    )
  }
  
  # Session info
  lines <- c(lines,
             "## Session Information",
             "",
             "```{r session_info}",
             "sessionInfo()",
             "```",
             ""
  )
  
  return(lines)
}

#' Export as Plain Text
#' @keywords internal
.export_as_text <- function(model, include_data, include_diagnostics, is_fitted) {
  
  lines <- character(0)
  
  # Header
  lines <- c(lines,
             "=============================================================================",
             "qbrms Model Specification",
             sprintf("Generated: %s", Sys.time()),
             "=============================================================================",
             ""
  )
  
  # Extract components
  if (is_fitted) {
    formula <- model$formula
    family_name <- model$family
  } else {
    formula <- model$formula
    family_name <- model$family$name
  }
  
  # Data section
  if (include_data) {
    lines <- c(lines,
               "DATA SUMMARY",
               "------------",
               if (is_fitted) {
                 c(sprintf("Observations: %d", nrow(model$data)),
                   sprintf("Variables: %d", ncol(model$data)),
                   sprintf("Response: %s", as.character(formula[[2]])))
               } else {
                 c(sprintf("Observations: %d", nrow(model$data)),
                   sprintf("Variables: %d", ncol(model$data)),
                   sprintf("Response: %s (%s)", model$response_info$name, 
                           model$response_info$type))
               },
               ""
    )
  }
  
  # Model specification
  lines <- c(lines,
             "MODEL SPECIFICATION",
             "-------------------",
             sprintf("Formula: %s", deparse(formula)),
             sprintf("Family: %s", family_name),
             ""
  )
  
  # Model code
  lines <- c(lines,
             "MODEL CODE",
             "----------",
             if (is_fitted) {
               sprintf("fit <- qbrms(\n  formula = %s,\n  data = model_data,\n  family = %s()\n)",
                       deparse(formula), family_name)
             } else {
               model$model_code
             },
             ""
  )
  
  # Results summary (if fitted)
  if (is_fitted) {
    lines <- c(lines,
               "MODEL FITTED",
               "------------",
               "Use summary(fit) to view detailed results",
               "Use plot(fit) to create diagnostic plots",
               ""
    )
  }
  
  return(lines)
}

#' Export as JSON
#' @keywords internal
.export_as_json <- function(model, include_data, is_fitted) {
  
  # Extract components
  if (is_fitted) {
    formula_str <- deparse(model$formula)
    family_name <- model$family
    data_summary <- list(
      n_obs = nrow(model$data),
      n_vars = ncol(model$data),
      response = as.character(model$formula[[2]])
    )
  } else {
    formula_str <- deparse(model$formula)
    family_name <- model$family$name
    data_summary <- list(
      n_obs = nrow(model$data),
      n_vars = ncol(model$data),
      response = model$response_info$name,
      response_type = model$response_info$type
    )
  }
  
  # Create JSON structure
  json_obj <- list(
    metadata = list(
      generated = as.character(Sys.time()),
      package = "qbrms",
      model_type = if (is_fitted) "fitted" else "specification"
    ),
    model = list(
      formula = formula_str,
      family = family_name
    )
  )
  
  if (include_data) {
    json_obj$data <- data_summary
  }
  
  if (is_fitted) {
    # Add summary statistics if available
    json_obj$fitted <- TRUE
  }
  
  # Convert to JSON string
  json_str <- jsonlite::toJSON(json_obj, pretty = TRUE, auto_unbox = TRUE)
  
  return(as.character(json_str))
}

#' Import Model Specification from JSON
#'
#' @description
#' Import a previously exported model specification from JSON format.
#'
#' @param file Character string specifying JSON file path
#'
#' @return A list containing the model specification components
#'
#' @examples
#' \dontrun{
#' # Import model
#' spec <- import_model("my_model.json")
#'
#' # Recreate model
#' fit <- qbrms(
#'   formula = as.formula(spec$model$formula),
#'   data = my_data,
#'   family = get(spec$model$family)()
#' )
#' }
#'
#' @export
import_model <- function(file) {
  
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }
  
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for JSON import", call. = FALSE)
  }
  
  # Read JSON
  json_content <- readLines(file, warn = FALSE)
  spec <- jsonlite::fromJSON(paste(json_content, collapse = "\n"))
  
  message("Model specification imported from: ", file)
  message("Formula: ", spec$model$formula)
  message("Family: ", spec$model$family)
  
  return(spec)
}