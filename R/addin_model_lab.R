# =============================================================================
# R/addin_model_lab.R
# =============================================================================
#' qbrms Model Lab (RStudio Add-in)
#'
#' @description
#' Compare plausible families, run prior/posterior checks, plot conditional
#' effects, compute diagnostics, and emit reproducible code. Does not load 'brms'.
#'
#' @return This function is called for its side effects (launching a Shiny
#'   gadget in RStudio). It returns \code{NULL} invisibly.
#'
#' @importFrom graphics hist barplot
#' @export
model_lab_addin <- function() {
  if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable())
    stop("This add-in requires RStudio.", call. = FALSE)
  if (!requireNamespace("shiny",  quietly = TRUE)) stop("Package 'shiny' is required.",  call. = FALSE)
  if (!requireNamespace("miniUI", quietly = TRUE)) stop("Package 'miniUI' is required.", call. = FALSE)
  
  have_gg <- requireNamespace("ggplot2", quietly = TRUE)
  
  # Canonical keys accepted by your converter/qbrm()
  supported_families <- c(
    "gaussian", "student", "skew_normal",
    "lognormal", "gamma",
    "poisson", "negbinomial",
    "binomial", "beta_binomial",
    "beta", "simplex",
    "zero_inflated_poisson", "zero_inflated_negbinomial",
    "weibull", "exponential",
    "multinomial", "gev", "circular_normal"
  )
  
  family_label <- function(f) switch(
    f,
    gaussian   = "Gaussian",
    student    = "Student t",
    skew_normal= "Skew-normal",
    lognormal  = "Log-normal",
    gamma      = "Gamma",
    poisson    = "Poisson",
    negbinomial= "Neg. binomial",
    binomial   = "Binomial",
    beta_binomial = "Beta-binomial",
    beta       = "Beta",
    simplex    = "Simplex",
    zero_inflated_poisson     = "Zero-inflated Poisson",
    zero_inflated_negbinomial = "Zero-inflated NegBin",
    weibull    = "Weibull",
    exponential= "Exponential",
    multinomial= "Multinomial",
    gev        = "GEV",
    circular_normal = "Circular normal",
    f
  )
  
  # Named choices used both at creation and on every update (prevents “disappearing” items)
  family_choices_map <- setNames(supported_families, vapply(supported_families, family_label, ""))
  
  # ------------------------------- UI -----------------------------------------
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("qbrms Model Lab"),
    miniUI::miniTabstripPanel(
      
      miniUI::miniTabPanel("Data & Spec", icon = shiny::icon("table"),
                           miniUI::miniContentPanel(
                             shiny::h4("Step 1: Select data and build a formula"),
                             
                             shiny::radioButtons("data_source", "Data source:",
                                                 choices = c("Existing dataset" = "existing", "Upload CSV" = "upload"),
                                                 selected = "existing", inline = TRUE
                             ),
                             
                             shiny::conditionalPanel("input.data_source == 'existing'",
                                                     shiny::selectInput("dataset", "Choose dataset:",
                                                                        choices = .get_dataframes_in_env(), multiple = FALSE
                                                     )
                             ),
                             
                             shiny::conditionalPanel("input.data_source == 'upload'",
                                                     shiny::fileInput("file_upload", "Upload CSV", accept = c(".csv", ".txt")),
                                                     shiny::checkboxInput("header", "Header", TRUE),
                                                     shiny::radioButtons("sep", "Separator:",
                                                                         c("Comma" = ",", "Semicolon" = ";", "Tab" = "\t"),
                                                                         selected = ",", inline = TRUE)
                             ),
                             
                             shiny::hr(),
                             shiny::uiOutput("response_selector"),
                             shiny::uiOutput("predictor_selector"),
                             shiny::div(style = "margin-top:6px;",
                                        shiny::actionLink("sel_all", "Select all"),
                                        " / ",
                                        shiny::actionLink("sel_none", "Deselect all")
                             ),
                             
                             shiny::hr(),
                             shiny::selectizeInput(
                               "families", "Candidate families",
                               choices  = family_choices_map,
                               multiple = TRUE,
                               selected = c("gaussian", "student"),
                               options = list(plugins = list("remove_button"))
                             ),
                             shiny::div(style = "margin-top:6px;",
                                        shiny::actionLink("apply_reco", "Use recommended families")
                             ),
                             shiny::verbatimTextOutput("family_reco_text"),
                             
                             shiny::hr(),
                             shiny::actionButton("fit_models", "Fit selected families", icon = shiny::icon("play")),
                             shiny::verbatimTextOutput("spec_summary")
                           )
      ),
      
      miniUI::miniTabPanel("Compare", icon = shiny::icon("balance-scale"),
                           miniUI::miniContentPanel(
                             shiny::h4("Model comparison across families"),
                             shiny::tableOutput("cmp_table"),
                             shiny::textOutput("cmp_best"),
                             shiny::verbatimTextOutput("cmp_notes")
                           )
      ),
      
      miniUI::miniTabPanel("Priors", icon = shiny::icon("sliders"),
                           miniUI::miniContentPanel(
                             shiny::h4("Prior visualisation and prior predictive checks"),
                             shiny::actionButton("run_prior_pp", "Run prior predictive check"),
                             shiny::checkboxInput("show_obs_in_prior", "Overlay observed data on prior check", FALSE),
                             shiny::plotOutput("prior_plot", height = "320px"),
                             shiny::verbatimTextOutput("prior_msg")
                           )
      ),
      
      miniUI::miniTabPanel("Posterior checks", icon = shiny::icon("chart-line"),
                           miniUI::miniContentPanel(
                             shiny::h4("Posterior predictive checks"),
                             shiny::selectInput("pp_type", "Check type",
                                                choices = c("dens_overlay", "hist", "scatter", "scatter_avg"),
                                                selected = "dens_overlay"
                             ),
                             shiny::actionButton("run_pp_check", "Create posterior check plot"),
                             shiny::plotOutput("pp_plot", height = "360px"),
                             shiny::verbatimTextOutput("pp_msg")
                           )
      ),
      
      miniUI::miniTabPanel("Conditional effects", icon = shiny::icon("project-diagram"),
                           miniUI::miniContentPanel(
                             shiny::h4("Conditional effects"),
                             shiny::uiOutput("cond_effect_var"),
                             shiny::actionButton("run_cond_eff", "Plot conditional effect"),
                             shiny::plotOutput("cond_plot", height = "360px"),
                             shiny::verbatimTextOutput("cond_msg")
                           )
      ),
      
      miniUI::miniTabPanel("Diagnostics", icon = shiny::icon("stethoscope"),
                           miniUI::miniContentPanel(
                             shiny::h4("Diagnostics"),
                             shiny::actionButton("run_diag", "Run diagnostics"),
                             shiny::verbatimTextOutput("diag_text"),
                             shiny::plotOutput("diag_plot", height = "380px")
                           )
      ),
      
      miniUI::miniTabPanel("Code", icon = shiny::icon("code"),
                           miniUI::miniContentPanel(
                             shiny::h4("Reproducible code for the selected or best model"),
                             shiny::verbatimTextOutput("model_code"),
                             shiny::actionButton("insert_code", "Insert at Cursor", icon = shiny::icon("arrow-right"))
                           )
      )
    )
  )
  
  # ----------------------------- Server ---------------------------------------
  server <- function(input, output, session) {
    rv <- shiny::reactiveValues(
      data = NULL, response = NULL, predictors = character(0),
      fits = list(), cmp = NULL, best_key = NULL, cmp_error = NULL,
      reco_families = NULL, reco_note = NULL
    )
    
    # Load data
    shiny::observe({
      if (input$data_source == "existing") {
        shiny::req(input$dataset)
        rv$data <- get(input$dataset, envir = .GlobalEnv)
      } else if (input$data_source == "upload") {
        shiny::req(input$file_upload)
        rv$data <- tryCatch({
          utils::read.csv(
            input$file_upload$datapath,
            header = isTRUE(input$header),
            sep    = input$sep,
            stringsAsFactors = TRUE
          )
        }, error = function(e) {
          shiny::showNotification(e$message, type = "error")
          NULL
        })
      }
    })
    
    # Reset selectors when data changes
    shiny::observeEvent(rv$data, {
      if (is.null(rv$data)) return(NULL)
      cols <- names(rv$data)
      if (!length(cols)) return(NULL)
      # Response
      shiny::updateSelectInput(session, "response", choices = cols, selected = cols[1])
      # Predictors
      shiny::updateSelectizeInput(session, "predictors",
                                  choices = cols, selected = character(0), server = TRUE
      )
      # Families: reset to defaults to avoid empty state
      shiny::updateSelectizeInput(session, "families",
                                  choices = family_choices_map, selected = c("gaussian", "student"), server = TRUE
      )
    }, ignoreInit = TRUE)
    
    # Variable selectors
    output$response_selector <- shiny::renderUI({
      shiny::req(rv$data)
      shiny::selectInput("response", "Response variable:", names(rv$data))
    })
    
    output$predictor_selector <- shiny::renderUI({
      shiny::req(rv$data, input$response)
      vars <- setdiff(names(rv$data), input$response)
      shiny::selectizeInput(
        "predictors", "Predictors:",
        choices = vars, selected = character(0), multiple = TRUE,
        options = list(plugins = list("remove_button"))
      )
    })
    
    # Keep predictors always valid
    shiny::observeEvent(list(rv$data, input$response), {
      shiny::req(rv$data, input$response)
      available <- setdiff(names(rv$data), input$response)
      sel <- intersect(input$predictors %||% character(0), available)
      shiny::updateSelectizeInput(session, "predictors",
                                  choices = available, selected = sel, server = TRUE
      )
    }, ignoreInit = TRUE)
    
    # Select all / none
    shiny::observeEvent(input$sel_all,  {
      shiny::req(rv$data, input$response)
      vars <- setdiff(names(rv$data), input$response)
      shiny::updateSelectizeInput(session, "predictors", selected = vars, server = TRUE)
    })
    shiny::observeEvent(input$sel_none, {
      shiny::updateSelectizeInput(session, "predictors", selected = character(0), server = TRUE)
    })
    
    # Recommend families on response change
    shiny::observe({
      shiny::req(rv$data, input$response)
      y <- rv$data[[input$response]]
      rec <- .infer_recommended_families(y)
      rec_ok <- intersect(rec$families, supported_families)
      if (!length(rec_ok)) rec_ok <- c("gaussian", "student")
      rv$reco_families <- unique(rec_ok)
      rv$reco_note <- rec$note
      # IMPORTANT: always pass 'choices' when updating (prevents disappearing selections)
      shiny::updateSelectizeInput(
        session, "families",
        choices  = family_choices_map,
        selected = rv$reco_families,
        server   = TRUE
      )
    })
    
    # If a user types an unsupported value, drop it immediately and keep choices visible
    shiny::observeEvent(input$families, {
      keep <- intersect(input$families %||% character(0), supported_families)
      if (!identical(keep, input$families)) {
        shiny::updateSelectizeInput(session, "families",
                                    choices = family_choices_map, selected = keep, server = TRUE
        )
      }
    }, ignoreInit = TRUE)
    
    # Apply recommendation explicitly
    shiny::observeEvent(input$apply_reco, {
      if (length(rv$reco_families)) {
        shiny::updateSelectizeInput(session, "families",
                                    choices = family_choices_map, selected = rv$reco_families, server = TRUE
        )
      }
    })
    
    # Recommendation text
    output$family_reco_text <- shiny::renderText({
      if (!length(rv$reco_families)) return("No recommendation available yet.")
      paste0(
        "Recommended families based on the response: ",
        paste(vapply(rv$reco_families, family_label, ""), collapse = ", "),
        if (!is.null(rv$reco_note) && nzchar(rv$reco_note)) paste0("\nReason: ", rv$reco_note) else ""
      )
    })
    
    # Spec summary
    output$spec_summary <- shiny::renderText({
      shiny::req(rv$data, input$response)
      rhs <- if (length(input$predictors)) paste(input$predictors, collapse = " + ") else "1"
      form <- paste(input$response, "~", rhs)
      fams <- input$families
      paste0("Formula: ", form, "\nFamilies: ", if (length(fams)) paste(fams, collapse = ", ") else "(none selected)")
    })
    
    # -------------------- Fit models across families and compare ----------------
    shiny::observeEvent(input$fit_models, {
      shiny::req(rv$data, input$response)
      shiny::validate(shiny::need(length(input$families) >= 1, "Select at least one family."))
      
      # Ensure predictors are valid at click time
      avail <- setdiff(names(rv$data), input$response)
      sel   <- intersect(input$predictors %||% character(0), avail)
      
      rhs  <- if (length(sel)) paste(sel, collapse = " + ") else "1"
      form <- stats::as.formula(paste(input$response, "~", rhs))
      
      # Reset
      rv$fits <- list(); rv$cmp <- NULL; rv$best_key <- NULL; rv$cmp_error <- NULL
      
      # Fit per family (normalise "t" -> "student")
      for (fam in input$families) {
        fam_low <- tolower(fam)
        fam_in  <- if (identical(fam_low, "t")) "student" else fam_low
        key <- fam_in
        
        fit <- tryCatch({
          # Merge any user-specified control.compute with the IC flags we need
          cc_user <- if (exists("control.compute", inherits = FALSE)) control.compute else NULL
          cc <- modifyList(list(waic = TRUE, dic = TRUE, cpo = TRUE), if (is.list(cc_user)) cc_user else list())
          
          qbrm(
            formula = form,
            data    = rv$data,
            family  = fam_in,
            control.compute = cc
          )
        }, error = function(e) e)
        
        if (inherits(fit, "error")) {
          shiny::showNotification(
            paste0("Fit failed for '", fam, "': ", fit$message),
            type = "error", duration = 10
          )
        } else {
          rv$fits[[key]] <- fit
        }
      }
      
      # Compare if two or more succeeded
      if (length(rv$fits) >= 2) {
        rv$cmp <- tryCatch(
          do.call(
            compare_models,
            c(rv$fits, list(criterion = "auto", compare_predictions = TRUE, weights = TRUE))
          ),
          error = function(e) { rv$cmp_error <- e$message; NULL }
        )
        
        if (!is.null(rv$cmp) && !is.null(rv$cmp$best_model)) {
          rv$best_key <- rv$cmp$best_model
        } else if (length(rv$fits) > 0) {
          rv$best_key <- names(rv$fits)[1]
        }
        
        if (!is.null(rv$cmp) && is.list(rv$cmp$criteria) && !is.null(rv$cmp$criteria$values)) {
          vals <- rv$cmp$criteria$values
          all_na <- is.data.frame(vals) && all(vapply(vals, function(col) all(is.na(col)), logical(1)))
          if (isTRUE(all_na)) {
            shiny::showNotification(
              "Information criteria are NA (likely a non-INLA fallback fit). Adjust specification and refit.",
              type = "warning", duration = 10
            )
          }
        }
      } else if (length(rv$fits) == 1) {
        rv$best_key <- names(rv$fits)[1]
      }
    })
    
    # Comparison table and notes
    output$cmp_table <- shiny::renderTable({
      shiny::validate(
        shiny::need(!is.null(rv$cmp), rv$cmp_error %||% "Fit at least two models to compare.")
      )
      as.data.frame(rv$cmp$comparison_table)
    }, striped = TRUE, hover = TRUE)
    
    output$cmp_best <- shiny::renderText({
      if (is.null(rv$cmp)) return("")
      paste("Best model by", toupper(rv$cmp$criterion_used), ":", rv$cmp$best_model)
    })
    
    output$cmp_notes <- shiny::renderText({
      if (!is.null(rv$cmp_error)) rv$cmp_error else ""
    })
    
    # Priors
    output$prior_plot <- shiny::renderPlot({
      shiny::req(rv$data, input$response)
      if (is.null(rv$best_key) || is.null(rv$fits[[rv$best_key]])) return(NULL)
      try({
        if (exists("prior_predictive_check", where = asNamespace("qbrms"), inherits = FALSE)) {
          rhs  <- if (length(input$predictors)) paste(input$predictors, collapse = " + ") else "1"
          form <- stats::as.formula(paste(input$response, "~", rhs))
          pfun <- get("prior_predictive_check", asNamespace("qbrms"))
          p <- pfun(formula = form, data = rv$data, family = rv$best_key, prior = NULL, n_samples = 1000)
          if (inherits(p, "ggplot")) print(p)
        }
      }, silent = TRUE)
    })
    
    output$prior_msg <- shiny::renderText({
      if (!have_gg) "Install 'ggplot2' to see prior plots." else ""
    })
    
    shiny::observeEvent(input$run_prior_pp, { output$prior_plot <- output$prior_plot })
    
    # Posterior predictive checks
    output$pp_plot <- shiny::renderPlot({
      shiny::req(rv$best_key, rv$fits[[rv$best_key]])
      try({
        if (exists("pp_check", where = asNamespace("qbrms"), inherits = TRUE)) {
          p <- get("pp_check", asNamespace("qbrms"))(rv$fits[[rv$best_key]], type = input$pp_type)
          if (inherits(p, "ggplot")) print(p)
        } else {
          p <- pp_check(rv$fits[[rv$best_key]], type = input$pp_type)
          if (inherits(p, "ggplot")) print(p)
        }
      }, silent = TRUE)
    })
    
    output$pp_msg <- shiny::renderText({
      if (!have_gg) "Install 'ggplot2' to see posterior plots." else ""
    })
    
    shiny::observeEvent(input$run_pp_check, { output$pp_plot <- output$pp_plot })
    
    # Conditional effects
    output$cond_effect_var <- shiny::renderUI({
      shiny::req(rv$data, input$response)
      choices <- setdiff(names(rv$data), input$response)
      shiny::selectInput("cond_var", "Effect of:",
                         choices = choices,
                         selected = if (length(choices)) choices[1] else NULL
      )
    })
    
    output$cond_plot <- shiny::renderPlot({
      shiny::req(rv$best_key, rv$fits[[rv$best_key]])
      shiny::req(!is.null(input$cond_var), nzchar(input$cond_var))
      try({
        if (exists("conditional_effects", where = asNamespace("qbrms"), inherits = FALSE)) {
          fun <- get("conditional_effects", asNamespace("qbrms"))
          p <- fun(rv$fits[[rv$best_key]], effect = input$cond_var)
          if (inherits(p, "ggplot")) print(p) else if (!is.null(p)) print(p)
        } else {
          if (!is.null(rv$data[[input$cond_var]]) && !is.null(rv$data[[input$response]])) {
            plot(rv$data[[input$cond_var]], rv$data[[input$response]], pch = 16,
                 main = paste("Effect of", input$cond_var),
                 xlab = input$cond_var, ylab = input$response)
          }
        }
      }, silent = TRUE)
    })
    
    output$cond_msg <- shiny::renderText({
      if (!have_gg) "Install 'ggplot2' for enhanced conditional plots." else ""
    })
    
    shiny::observeEvent(input$run_cond_eff, { output$cond_plot <- output$cond_plot })
    
    # Diagnostics
    output$diag_text <- shiny::renderText({ "" })
    output$diag_plot <- shiny::renderPlot({ NULL })
    
    shiny::observeEvent(input$run_diag, {
      shiny::req(rv$best_key, rv$fits[[rv$best_key]])
      dg <- tryCatch({
        if (exists("diagnose_model", where = asNamespace("qbrms"), inherits = FALSE)) {
          get("diagnose_model", asNamespace("qbrms"))(rv$fits[[rv$best_key]])
        } else {
          stop("diagnose_model() not found in qbrms.")
        }
      }, error = function(e) e)
      
      if (inherits(dg, "error")) {
        output$diag_text <- shiny::renderText(dg$message)
        output$diag_plot <- shiny::renderPlot(NULL)
      } else {
        output$diag_text <- shiny::renderText(paste(capture.output(print(dg)), collapse = "\n"))
        output$diag_plot <- shiny::renderPlot({ try(print(plot(dg)), silent = TRUE) })
      }
    })
    
    # Code export
    output$model_code <- shiny::renderText({
      shiny::req(rv$data, input$response, rv$best_key)
      fam_print <- switch(rv$best_key,
                          "student"    = "student()",
                          "gaussian"   = "gaussian()",
                          "skew_normal"= "\"skew_normal\"",
                          "gamma"      = "gamma()",
                          "lognormal"  = "qbrms::lognormal()",
                          "binomial"   = "binomial()",
                          "beta_binomial" = "qbrms::beta_binomial()",
                          "poisson"    = "poisson()",
                          "negbinomial"= "\"negbinomial\"",
                          "beta"       = "qbrms::Beta()",
                          "simplex"    = "qbrms::simplex()",
                          "zero_inflated_poisson"     = "qbrms::zero_inflated_poisson()",
                          "zero_inflated_negbinomial" = "qbrms::zero_inflated_negbinomial()",
                          "weibull"    = "qbrms::weibull()",
                          "exponential"= "qbrms::exponential()",
                          "multinomial"= "\"multinomial\"",
                          "gev"        = "qbrms::gev()",
                          "circular_normal" = "qbrms::circular_normal()",
                          paste0("\"", rv$best_key, "\"")
      )
      avail <- setdiff(names(rv$data), input$response)
      sel   <- intersect(input$predictors %||% character(0), avail)
      rhs   <- if (length(sel)) paste(sel, collapse = " + ") else "1"
      form_str <- paste(input$response, "~", rhs)
      data_sym <- if (!is.null(input$dataset) && nzchar(input$dataset)) input$dataset else "uploaded_data"
      
      paste0(
        "# qbrms model specification (generated by Model Lab)\n",
        "fit <- qbrm(\n",
        "  formula = ", form_str, ",\n",
        "  data    = ", data_sym, ",\n",
        "  family  = ", fam_print, "\n",
        ")\n\n",
        "summary(fit)\n",
        "# pp_check(fit, type = \"dens_overlay\")\n",
        "# diagnose_model(fit)\n"
      )
    })
    
    shiny::observeEvent(input$insert_code, {
      if (requireNamespace("rstudioapi", quietly = TRUE)) {
        rstudioapi::insertText(text = output$model_code() %||% "")
        shiny::showNotification("Code inserted at cursor.", type = "message")
      }
    })
    
    # Done / Cancel
    shiny::observeEvent(input$done,   { shiny::stopApp(invisible(NULL)) })
    shiny::observeEvent(input$cancel, { shiny::stopApp(invisible(NULL)) })
  }
  
  viewer <- shiny::dialogViewer("qbrms Model Lab", width = 1000, height = 820)
  shiny::runGadget(ui, server, viewer = viewer)
}

# Helper
`%||%` <- function(x, y) if (is.null(x) || (is.character(x) && !length(x))) y else x

# --------------------------- Family recommendation ----------------------------

#' Infer recommended families from a response vector
#' @keywords internal
.infer_recommended_families <- function(y, tol = 1e-8) {
  note <- ""
  if (is.factor(y) || is.character(y)) {
    f <- if (is.factor(y)) y else factor(y)
    k <- length(levels(f))
    if (k == 2L) {
      return(list(families = c("binomial", "beta_binomial"), note = "Binary outcome detected."))
    } else {
      if (is.ordered(f)) {
        return(list(
          families = c("multinomial"),
          note = "Ordered categorical outcome detected. Use qbrmO for ordinal models where available."
        ))
      } else {
        return(list(families = c("multinomial"), note = "Unordered categorical outcome detected."))
      }
    }
  }
  
  if (!is.numeric(y)) {
    return(list(families = c("gaussian", "student"), note = "Unknown type; defaulting to robust continuous families."))
  }
  
  y <- y[!is.na(y)]
  if (!length(y)) {
    return(list(families = c("gaussian", "student"), note = "No non-missing values; using defaults."))
  }
  
  is_int     <- all(abs(y - round(y)) < tol)
  all_nonneg <- all(y >= 0)
  all_pos    <- all(y > 0)
  in_open01  <- all(y > 0 & y < 1)
  in_closed01<- all(y >= 0 & y <= 1)
  zero_frac  <- if (length(y)) mean(y == 0) else 0
  nuniq      <- length(unique(y))
  
  if (is_int && all_nonneg) {
    fam <- c("poisson", "negbinomial")
    if (is.finite(zero_frac) && zero_frac > 0.10) {
      fam <- c(fam, "zero_inflated_poisson", "zero_inflated_negbinomial")
      note <- "Count data with excess zeros detected."
    } else {
      note <- "Count data detected."
    }
    return(list(families = unique(fam), note = note))
  }
  
  if (in_open01 || in_closed01) {
    note <- if (in_closed01) "Proportions in [0,1] with boundary values; consider slight adjustment to (0,1)." else "Proportions in (0,1) detected."
    return(list(families = c("beta", "simplex"), note = note))
  }
  
  if (all_pos) {
    return(list(families = c("lognormal", "gamma", "student", "gaussian"), note = "Positive continuous response detected."))
  }
  
  # Unbounded continuous
  fam <- c("gaussian", "student", "skew_normal")
  if (nuniq < 20) fam <- c("gaussian", "student")
  return(list(families = fam, note = "Continuous unbounded response detected."))
}

# Discover data frames in the global environment
#' @keywords internal
.get_dataframes_in_env <- function() {
  objs <- ls(envir = .GlobalEnv)
  if (!length(objs)) return(character(0))
  is_df <- vapply(objs, function(x) {
    obj <- get(x, envir = .GlobalEnv)
    is.data.frame(obj)
  }, logical(1))
  objs[is_df]
}