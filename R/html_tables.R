#' Create HTML Table for qbrms Models with Enhanced Styling
#'
#' @description
#' Generate APA-style HTML tables for qbrms model outputs with customizable styling options.
#'
#' @param ... One or more qbrms_fit objects to display in the table
#' @param show.ci Logical; show credible intervals (default TRUE)
#' @param ci.lvl Credible interval level (default 0.95)
#' @param show.rope Logical; show ROPE analysis (default FALSE)
#' @param rope Numeric vector c(lower, upper) for ROPE bounds
#' @param show.p_sig Logical; show probability of practical significance (default FALSE)
#' @param show.pd Logical; show probability of direction (default FALSE)
#' @param show.bf Logical; show Bayes factors (default FALSE)
#' @param digits Number of decimal places (default 2)
#' @param title Character; table title
#' @param file Character; file path to save HTML output (optional)
#' @param CSS Character; custom CSS styling (optional)
#' @param dv.labels Character vector of dependent variable labels
#' @param pred.labels Named character vector for predictor labels
#' @param show.intercept Logical; show intercept row (default TRUE)
#' @param show.r2 Logical; show R-squared if available (default FALSE)
#' @param show.icc Logical; show ICC for mixed models (default FALSE)
#' @param show.nobs Logical; show number of observations (default TRUE)
#' @param bootstrap Logical; use Bootstrap CSS framework (default TRUE)
#' @param table.style Character; table style theme. Options: "default", "minimal", "academic", "modern"
#' @param font.family Character; CSS font family (default "system-ui")
#' @param font.size Character; base font size (default "14px")
#' @param header.bg Character; header background colour (default "#f8f9fa")
#' @param stripe.bg Character; striped row background colour (default "#f9f9f9")
#' @param verbose Logical; print progress (default FALSE)
#'
#' @return An object of class "qbrms_html_table" containing the HTML code
#'
#' @export
tab_model <- function(...,
                      show.ci = TRUE,
                      ci.lvl = 0.95,
                      show.rope = FALSE,
                      rope = c(-0.1, 0.1),
                      show.p_sig = FALSE,
                      show.pd = FALSE,
                      show.bf = FALSE,
                      digits = 2,
                      title = "Model Results",
                      file = NULL,
                      CSS = NULL,
                      dv.labels = NULL,
                      pred.labels = NULL,
                      show.intercept = TRUE,
                      show.r2 = FALSE,
                      show.icc = FALSE,
                      show.nobs = TRUE,
                      bootstrap = TRUE,
                      table.style = "default",
                      font.family = "system-ui, -apple-system, sans-serif",
                      font.size = "14px",
                      header.bg = "#f8f9fa",
                      stripe.bg = "#f9f9f9",
                      verbose = FALSE) {
  
  # Get model objects
  models <- list(...)
  n_models <- length(models)
  
  # Validate inputs
  if (n_models == 0) stop("No models provided")
  
  # Check all are qbrms_fit objects
  if (!all(sapply(models, inherits, "qbrms_fit"))) {
    stop("All objects must be qbrms_fit objects")
  }
  
  if (verbose) cat("Creating HTML table for", n_models, "model(s)\n")
  
  # Extract model information
  model_info <- lapply(models, extract_model_info, 
                       ci.lvl = ci.lvl, 
                       rope = rope,
                       show.rope = show.rope,
                       show.p_sig = show.p_sig,
                       show.pd = show.pd,
                       show.bf = show.bf,
                       verbose = verbose)
  
  # Get all unique predictors across models
  all_predictors <- unique(unlist(lapply(model_info, function(x) x$predictors)))
  
  # Apply predictor labels if provided
  if (!is.null(pred.labels)) {
    display_predictors <- ifelse(all_predictors %in% names(pred.labels),
                                 pred.labels[all_predictors],
                                 all_predictors)
  } else {
    display_predictors <- all_predictors
  }
  
  # Apply DV labels if provided
  if (!is.null(dv.labels)) {
    if (length(dv.labels) != n_models) {
      warning("Length of dv.labels doesn't match number of models")
      dv.labels <- paste("Model", 1:n_models)
    }
  } else {
    dv.labels <- paste("Model", 1:n_models)
  }
  
  # Build HTML with styling
  html <- build_html_table_styled(model_info = model_info,
                                  all_predictors = all_predictors,
                                  display_predictors = display_predictors,
                                  dv.labels = dv.labels,
                                  show.ci = show.ci,
                                  show.rope = show.rope,
                                  show.p_sig = show.p_sig,
                                  show.pd = show.pd,
                                  show.bf = show.bf,
                                  show.intercept = show.intercept,
                                  show.r2 = show.r2,
                                  show.icc = show.icc,
                                  show.nobs = show.nobs,
                                  digits = digits,
                                  title = title,
                                  bootstrap = bootstrap,
                                  table.style = table.style,
                                  font.family = font.family,
                                  font.size = font.size,
                                  header.bg = header.bg,
                                  stripe.bg = stripe.bg,
                                  CSS = CSS,
                                  verbose = verbose)
  
  # Save to file if requested
  if (!is.null(file)) {
    if (verbose) cat("Saving table to", file, "\n")
    writeLines(html, file)
  }
  
  # Create return object
  result <- list(
    html = html,
    n_models = n_models,
    file = file,
    style = table.style
  )
  
  class(result) <- "qbrms_html_table"
  return(result)
}

#' Build HTML Table with Enhanced Styling
#' @keywords internal
build_html_table_styled <- function(model_info, all_predictors, display_predictors,
                                    dv.labels, show.ci, show.rope, show.p_sig,
                                    show.pd, show.bf, show.intercept, show.r2,
                                    show.icc, show.nobs, digits, title,
                                    bootstrap, table.style, font.family, font.size,
                                    header.bg, stripe.bg, CSS, verbose = FALSE) {
  
  n_models <- length(model_info)
  
  # Generate CSS based on style
  style_css <- generate_table_css(table.style, font.family, font.size, 
                                  header.bg, stripe.bg)
  
  # Start HTML
  html <- character()
  
  # Add HTML header
  if (bootstrap) {
    html <- c(html, 
              '<!DOCTYPE html>',
              '<html>',
              '<head>',
              '<meta charset="utf-8">',
              '<meta name="viewport" content="width=device-width, initial-scale=1">',
              paste0('<title>', title, '</title>'),
              '<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">',
              '<style>',
              style_css,
              if (!is.null(CSS)) CSS else '',
              '</style>',
              '</head>',
              '<body>')
  }
  
  # Start table with style-specific classes
  table_class <- switch(table.style,
                        "minimal" = "table table-minimal",
                        "academic" = "table table-academic", 
                        "modern" = "table table-modern",
                        "table table-striped table-hover")  # default
  
  html <- c(html,
            '<div class="container">',
            paste0('<h3 class="table-title">', title, '</h3>'),
            paste0('<table class="', table_class, '">'))
  
  # Table header
  html <- c(html, '<thead>')
  
  # First header row - model names
  header1 <- '<tr><th rowspan="2">Predictors</th>'
  for (i in 1:n_models) {
    colspan <- 1  # Base estimate column
    if (show.ci) colspan <- colspan + 1
    if (show.rope) colspan <- colspan + 1
    if (show.p_sig) colspan <- colspan + 1
    if (show.pd) colspan <- colspan + 1
    if (show.bf) colspan <- colspan + 1
    
    header1 <- paste0(header1, 
                      '<th colspan="', colspan, 
                      '" class="text-center model-header">', 
                      dv.labels[i], '</th>')
  }
  header1 <- paste0(header1, '</tr>')
  html <- c(html, header1)
  
  # Second header row - column names
  header2 <- '<tr>'
  for (i in 1:n_models) {
    header2 <- paste0(header2, '<th>Estimate</th>')
    if (show.ci) header2 <- paste0(header2, '<th>95% CI</th>')
    if (show.rope) header2 <- paste0(header2, '<th>ROPE %</th>')
    if (show.p_sig) header2 <- paste0(header2, '<th>p(sig)</th>')
    if (show.pd) header2 <- paste0(header2, '<th>pd</th>')
    if (show.bf) header2 <- paste0(header2, '<th>BF</th>')
  }
  header2 <- paste0(header2, '</tr>')
  html <- c(html, header2)
  html <- c(html, '</thead>')
  
  # Table body
  html <- c(html, '<tbody>')
  
  # Add rows for each predictor
  for (j in seq_along(all_predictors)) {
    pred <- all_predictors[j]
    
    # Skip intercept if requested
    if (!show.intercept && pred == "(Intercept)") next
    
    row <- paste0('<tr><td class="predictor-name"><strong>', display_predictors[j], '</strong></td>')
    
    for (i in 1:n_models) {
      if (pred %in% model_info[[i]]$predictors) {
        idx <- which(model_info[[i]]$predictors == pred)
        
        # Estimate
        est <- format_number(model_info[[i]]$estimates[idx], digits)
        se <- format_number(model_info[[i]]$se[idx], digits)
        row <- paste0(row, '<td class="estimate">', est)
        if (!is.na(model_info[[i]]$se[idx])) {
          row <- paste0(row, ' <span class="se">(', se, ')</span>')
        }
        row <- paste0(row, '</td>')
        
        # CI
        if (show.ci) {
          ci_low <- format_number(model_info[[i]]$ci_lower[idx], digits)
          ci_high <- format_number(model_info[[i]]$ci_upper[idx], digits)
          row <- paste0(row, '<td class="ci">[', ci_low, ', ', ci_high, ']</td>')
        }
        
        # ROPE
        if (show.rope) {
          rope_val <- format_percentage(model_info[[i]]$rope_pct[idx])
          row <- paste0(row, '<td class="rope">', rope_val, '</td>')
        }
        
        # p_significance
        if (show.p_sig) {
          p_sig_val <- format_number(model_info[[i]]$p_sig[idx], digits)
          row <- paste0(row, '<td class="p-sig">', p_sig_val, '</td>')
        }
        
        # Probability of direction
        if (show.pd) {
          pd_val <- format_number(model_info[[i]]$pd[idx], digits)
          row <- paste0(row, '<td class="pd">', pd_val, '</td>')
        }
        
        # Bayes factor
        if (show.bf) {
          bf_val <- format_bf(model_info[[i]]$bf[idx])
          row <- paste0(row, '<td class="bf">', bf_val, '</td>')
        }
        
      } else {
        # Empty cells for predictors not in this model
        ncols <- 1
        if (show.ci) ncols <- ncols + 1
        if (show.rope) ncols <- ncols + 1
        if (show.p_sig) ncols <- ncols + 1
        if (show.pd) ncols <- ncols + 1
        if (show.bf) ncols <- ncols + 1
        
        for (k in 1:ncols) {
          row <- paste0(row, '<td class="empty">-</td>')
        }
      }
    }
    
    row <- paste0(row, '</tr>')
    html <- c(html, row)
  }
  
  html <- c(html, '</tbody>')
  
  # Add footer with model statistics
  if (show.nobs || show.r2 || show.icc) {
    html <- c(html, '<tfoot>')
    
    if (show.nobs) {
      row <- '<tr><td class="footer-label"><em>N</em></td>'
      for (i in 1:n_models) {
        row <- paste0(row, '<td class="footer-stat">', model_info[[i]]$n_obs, '</td>')
        # Add empty cells for extra columns
        ncols <- 0
        if (show.ci) ncols <- ncols + 1
        if (show.rope) ncols <- ncols + 1
        if (show.p_sig) ncols <- ncols + 1
        if (show.pd) ncols <- ncols + 1
        if (show.bf) ncols <- ncols + 1
        for (k in 1:ncols) {
          row <- paste0(row, '<td class="footer-empty"></td>')
        }
      }
      row <- paste0(row, '</tr>')
      html <- c(html, row)
    }
    
    html <- c(html, '</tfoot>')
  }
  
  # Close table
  html <- c(html, '</table>', '</div>')
  
  # Close HTML if bootstrap
  if (bootstrap) {
    html <- c(html, '</body>', '</html>')
  }
  
  # Convert html vector to string
  html_content <- paste(html, collapse = '\n')
  
  # Display in RStudio viewer if interactive
  if (interactive()) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      temp_file <- tempfile(fileext = ".html")
      writeLines(html_content, temp_file)
      rstudioapi::viewer(temp_file)
    }
  }
  
  # Return the HTML content (remove the duplicate return)
  return(html_content)
}

#' Generate CSS for different table styles
#' @keywords internal
generate_table_css <- function(style, font.family, font.size, header.bg, stripe.bg) {
  
base_css <- paste0('
body { 
  padding: 20px; 
  font-family: ', font.family, ';
  font-size: ', font.size, ';
}

.container {
  max-width: 1200px;
  margin: 0 auto;
}

.table-title {
  margin-bottom: 20px;
  text-align: center;
}
')
  
  style_specific <- switch(style,
                           
                           "minimal" = paste0('
.table-minimal {
  border-collapse: collapse;
  width: 100%;
  background: white;
}

.table-minimal th,
.table-minimal td {
  padding: 8px 12px;
  text-align: left;
  border-bottom: 1px solid #ddd;
}

.table-minimal th {
  background: white;
  font-weight: 600;
  border-bottom: 2px solid #333;
}

.table-minimal .model-header {
  background: white !important;
  border-bottom: 1px solid #333;
}

.predictor-name {
  font-weight: 500;
}

.se {
  color: #666;
  font-size: 0.9em;
}
'),
                           
                           "academic" = paste0('
.table-academic {
  border-collapse: collapse;
  width: 100%;
  background: white;
  font-family: "Times New Roman", serif;
}

.table-academic th,
.table-academic td {
  padding: 6px 10px;
  text-align: left;
  border: none;
}

.table-academic th {
  background: white;
  font-weight: bold;
  border-top: 2px solid black;
  border-bottom: 1px solid black;
}

.table-academic .model-header {
  background: white !important;
  border-bottom: 1px solid black;
}

.table-academic tbody tr:last-child td {
  border-bottom: 2px solid black;
}

.predictor-name {
  font-style: italic;
}

.se {
  color: #333;
  font-size: 0.9em;
}
'),
                           
                           "modern" = paste0('
.table-modern {
  border-collapse: separate;
  border-spacing: 0;
  width: 100%;
  background: white;
  border-radius: 8px;
  overflow: hidden;
  box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
}

.table-modern th,
.table-modern td {
  padding: 12px 16px;
  text-align: left;
}

.table-modern th {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color: white;
  font-weight: 600;
  text-transform: uppercase;
  font-size: 0.85em;
  letter-spacing: 0.5px;
}

.table-modern .model-header {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
}

.table-modern tbody tr:nth-child(even) {
  background: #f8f9fa;
}

.table-modern tbody tr:hover {
  background: #e8f4fd;
  transition: background-color 0.3s ease;
}

.predictor-name {
  font-weight: 600;
  color: #2c3e50;
}

.estimate {
  font-weight: 500;
}

.se {
  color: #7f8c8d;
  font-size: 0.9em;
}

.ci { color: #34495e; }
.rope { color: #27ae60; }
.p-sig { color: #3498db; }
.pd { color: #9b59b6; }
.bf { color: #f39c12; }
'),
                           
                           # Default style
                           paste0('
.table {
  font-size: 0.9em;
  border-collapse: collapse;
  width: 100%;
}

.table th {
  background-color: ', header.bg, ';
  padding: 10px;
  font-weight: 600;
}

.table td {
  padding: 8px 10px;
}

.table-striped tbody tr:nth-of-type(odd) {
  background-color: ', stripe.bg, ';
}

.table-hover tbody tr:hover {
  background-color: rgba(0, 0, 0, 0.075);
}

.model-header {
  background-color: #e9ecef !important;
  font-weight: bold;
}

.se {
  color: #6c757d;
  font-size: 0.85em;
}

.ci { color: #6c757d; font-size: 0.85em; }
.rope { color: #28a745; }
.p-sig { color: #007bff; }
.pd { color: #17a2b8; }
.bf { color: #ffc107; }
')
  )
  
  return(paste0(base_css, style_specific))
}

#' Display HTML Table in Viewer
#' @param x A `qbrms_html_table` created by `tab_model()`.
#' @export
view_table <- function(x) { 
  if (!inherits(x, "qbrms_html_table")) {
    stop("Object must be created by tab_model()")
  }
  
  # Create temporary file
  temp_file <- tempfile(fileext = ".html")
  writeLines(x$html, temp_file)
  
  # Try to open in viewer or browser
  if (requireNamespace("rstudioapi", quietly = TRUE) && 
      rstudioapi::isAvailable()) {
    rstudioapi::viewer(temp_file)
  } else {
    utils::browseURL(temp_file)
  }
  
  invisible(x)
}

# =============================================================================
# Helper Functions for html_tables.R
# =============================================================================

#' Extract Model Information for HTML Table
#' @keywords internal
extract_model_info <- function(model, ci.lvl, rope, show.rope, show.p_sig, 
                               show.pd, show.bf, verbose = FALSE) {
  
  info <- list()
  
  # Get response variable
  info$response <- all.vars(model$original_formula)[1]
  
  # Get coefficients
  if (!is.null(model$fit$summary.fixed)) {
    coef_table <- model$fit$summary.fixed
    info$predictors <- rownames(coef_table)
    info$estimates <- coef_table[, "mean"]
    info$se <- coef_table[, "sd"]
    
    # Calculate CI
    z <- qnorm((1 + ci.lvl) / 2)
    info$ci_lower <- info$estimates - z * info$se
    info$ci_upper <- info$estimates + z * info$se
  } else {
    # Fallback for other model types
    coefs <- coef(model)
    info$predictors <- names(coefs)
    info$estimates <- coefs
    info$se <- rep(NA, length(coefs))
    info$ci_lower <- rep(NA, length(coefs))
    info$ci_upper <- rep(NA, length(coefs))
  }
  
  # Add Bayesian statistics if requested
  nsim <- 1000
  
  if (show.rope) {
    if (verbose) cat("  Calculating ROPE...\n")
    rope_res <- rope_analysis(model, parameters = info$predictors, 
                              rope = rope, nsim = nsim)
    info$rope_pct <- rope_res$ROPE_percentage
  }
  
  if (show.p_sig) {
    if (verbose) cat("  Calculating p_significance...\n")
    # FIXED: Use threshold instead of rope
    p_sig_res <- p_significance(model, parameters = info$predictors,
                                threshold = rope, nsim = nsim, verbose = FALSE)
    info$p_sig <- p_sig_res$ps
  }
  
  if (show.pd) {
    if (verbose) cat("  Calculating probability of direction...\n")
    pd_res <- p_direction(model, parameters = info$predictors, nsim = nsim)
    info$pd <- pd_res$PD
  }
  
  if (show.bf) {
    if (verbose) cat("  Calculating Bayes factors...\n")
    info$bf <- sapply(info$predictors, function(p) {
      tryCatch({
        bf_res <- bayesfactor(model, 
                              hypothesis = paste(p, "= 0"),
                              nsim = nsim,
                              verbose = FALSE)
        bf_res$bayes_factor
      }, error = function(e) NA)
    })
  }
  
  # Get model fit statistics
  info$n_obs <- nrow(model$data)
  
  # Extract family name safely
  if (!is.null(model$family)) {
    if (is.character(model$family)) {
      info$family <- model$family
    } else if (is.list(model$family) && !is.null(model$family$family)) {
      info$family <- model$family$family
    } else {
      info$family <- "unknown"
    }
  } else {
    info$family <- "gaussian"  # default
  }
  
  # Add R-squared if available (would need to implement)
  info$r2 <- NA
  
  # Add ICC for mixed models if available
  if (!is.null(model$group_var)) {
    info$icc <- NA  # Would need proper implementation
  }
  
  return(info)
}

#' Format Number for Display
#' @keywords internal
format_number <- function(x, digits = 2) {
  if (is.na(x)) return("-")
  if (abs(x) < 0.001) return("<0.001")
  sprintf(paste0("%.", digits, "f"), x)
}

#' Format Percentage for Display
#' @keywords internal
format_percentage <- function(x, digits = 1) {
  if (is.na(x)) return("-")
  paste0(sprintf(paste0("%.", digits, "f"), x * 100), "%")
}

#' Format Bayes Factor for Display
#' @keywords internal
format_bf <- function(x) {
  if (is.na(x)) return("-")
  if (x > 100) return(">100")
  if (x < 0.01) return("<0.01")
  if (x > 10) return(sprintf("%.0f", x))
  if (x > 1) return(sprintf("%.1f", x))
  return(sprintf("%.2f", x))
}

#' Extract family name safely
#' @keywords internal  
extract_family_name <- function(family) {
  if (is.null(family)) return("gaussian")
  
  if (is.character(family)) {
    return(family)
  } else if (is.list(family)) {
    if (!is.null(family$family)) {
      return(family$family)
    } else if (!is.null(family$link)) {
      return(paste("unknown", family$link, sep = "_"))
    }
  }
  
  return("unknown")
}
