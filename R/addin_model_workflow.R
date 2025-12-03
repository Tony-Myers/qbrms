# =============================================================================
# R/addin_model_workflow.R
# =============================================================================

#' Launch Guided Bayesian Workflow (RStudio Add-in)
#'
#' @description
#' A comprehensive, step-by-step assistant for Bayesian model building with qbrms.
#'
#' @return No return value. This function launches an interactive Shiny gadget 
#'   for model building and code generation.
#'
#' @export
model_workflow_addin <- function() {
  
  # --- Dependencies ---
  if (!requireNamespace("shiny", quietly = TRUE)) stop("Package 'shiny' is required.")
  if (!requireNamespace("miniUI", quietly = TRUE)) stop("Package 'miniUI' is required.")
  if (!requireNamespace("rstudioapi", quietly = TRUE)) stop("Package 'rstudioapi' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  
  # --- UI Definition ---
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("qbrms Bayesian Coach", 
                           left = miniUI::miniTitleBarCancelButton(),
                           right = miniUI::miniTitleBarButton("done", "Insert Code", primary = TRUE)),
    
    miniUI::miniTabstripPanel(
      
      # --- Tab 1: Data & Structure ---
      miniUI::miniTabPanel(
        "1. Data & Structure",
        icon = shiny::icon("database"),
        miniUI::miniContentPanel(
          shiny::fluidRow(
            shiny::column(4, 
                          shiny::wellPanel(
                            shiny::h4("Data Selection"),
                            shiny::radioButtons("data_source", "Source:", choices = c("Environment" = "env", "Load File" = "file"), inline = TRUE),
                            shiny::conditionalPanel(condition = "input.data_source == 'env'", shiny::selectInput("dataset_env", "Choose Data Frame:", choices = NULL)),
                            shiny::conditionalPanel(condition = "input.data_source == 'file'", 
                                                    shiny::actionButton("btn_browse", "Browse...", icon = shiny::icon("folder-open"), class = "btn-info", width = "100%"),
                                                    shiny::helpText("Supports: .csv, .xlsx, .rds, .dta, .sav"),
                                                    shiny::textOutput("file_path_display")
                            ),
                            shiny::hr(),
                            shiny::h4("Variables"),
                            shiny::selectInput("response", "Response Variable (Outcome):", choices = NULL),
                            shiny::selectInput("predictors", "Fixed Effects (Predictors):", choices = NULL, multiple = TRUE),
                            shiny::selectInput("grouping", "Random Effects (Grouping):", choices = NULL, multiple = TRUE)
                          )
            ),
            shiny::column(8,
                          shiny::div(style = "padding-left: 10px;",
                                     shiny::h4("Response Distribution"),
                                     shiny::plotOutput("response_dist_plot", height = "200px"),
                                     shiny::hr(),
                                     shiny::h4("Model Advice"),
                                     shiny::uiOutput("data_advice"),
                                     shiny::br(),
                                     shiny::h4("Likelihood Family"),
                                     shiny::helpText("Check the dropdown list below for additional family options not auto-selected."),
                                     shiny::uiOutput("family_selector"),
                                     shiny::uiOutput("family_description")
                          )
            )
          )
        )
      ),
      
      # --- Tab 2: Priors ---
      miniUI::miniTabPanel(
        "2. Priors",
        icon = shiny::icon("sliders-h"),
        miniUI::miniContentPanel(
          shiny::fluidRow(
            shiny::column(4,
                          shiny::wellPanel(
                            shiny::h4("Elicit Priors"),
                            shiny::helpText("Describe expectations before seeing data."),
                            shiny::numericInput("outcome_mean", "Expected Outcome Mean (approx):", value = 0),
                            shiny::numericInput("outcome_sd", "Expected Outcome SD (uncertainty):", value = 1, min = 0.1),
                            shiny::checkboxInput("standardise", "Standardise predictors?", value = TRUE),
                            shiny::hr(),
                            shiny::actionButton("check_priors", "Simulate & Check Priors", class = "btn-primary", width = "100%")
                          )
            ),
            shiny::column(8,
                          shiny::div(style = "padding-left: 10px;",
                                     shiny::h4("Prior Predictive Check"),
                                     shiny::plotOutput("prior_plot", height = "350px"),
                                     shiny::verbatimTextOutput("prior_summary")
                          )
            )
          )
        )
      ),
      
      # --- Tab 3: Fit & Compare ---
      miniUI::miniTabPanel(
        "3. Fit & Compare",
        icon = shiny::icon("check-circle"),
        miniUI::miniContentPanel(
          shiny::fluidRow(
            shiny::column(4,
                          shiny::wellPanel(
                            shiny::h4("Model Fitting"),
                            shiny::actionButton("fit_model", "Fit Model(s)", class = "btn-success", width = "100%", icon = shiny::icon("play")),
                            shiny::hr(),
                            shiny::h4("Status"),
                            shiny::verbatimTextOutput("fit_status")
                          )
            ),
            shiny::column(8,
                          shiny::div(style = "padding-left: 10px;",
                                     shiny::h4("Model Comparison"),
                                     shiny::verbatimTextOutput("comparison_output"),
                                     shiny::uiOutput("comparison_interpretation"),
                                     shiny::hr(),
                                     shiny::h4("Posterior Diagnostics (Best Model)"),
                                     shiny::plotOutput("pp_check_plot", height = "300px")
                          )
            )
          )
        )
      ),
      
      # --- Tab 4: Interpretation ---
      miniUI::miniTabPanel(
        "4. Interpretation",
        icon = shiny::icon("lightbulb"),
        miniUI::miniContentPanel(
          shiny::fluidRow(
            # Sidebar Column
            shiny::column(4,
                          shiny::wellPanel(
                            shiny::h4("Analysis Tools"),
                            shiny::radioButtons("interp_mode", "Analysis Mode:",
                                                choices = c("Hypothesis Testing" = "hyp", "Visualisation" = "vis")),
                            shiny::hr(),
                            
                            # Hypothesis Controls
                            shiny::conditionalPanel(
                              condition = "input.interp_mode == 'hyp'",
                              shiny::selectInput("post_analysis_type", "Method:", 
                                                 choices = c("Estimated Marginal Means (EMMs)", 
                                                             "Probability of Direction (pd)",
                                                             "Practical Significance",
                                                             "ROPE Analysis")),
                              
                              # EMMs specific input
                              shiny::conditionalPanel(
                                condition = "input.post_analysis_type == 'Estimated Marginal Means (EMMs)'",
                                shiny::selectInput("emm_term", "Select Term:", choices = NULL)
                              ),
                              
                              # ROPE/Practical Sig specific input
                              shiny::conditionalPanel(
                                condition = "input.post_analysis_type == 'Practical Significance' || input.post_analysis_type == 'ROPE Analysis'",
                                shiny::numericInput("rope_range", "Threshold (+/-):", value = 0.1, step = 0.01)
                              ),
                              
                              shiny::br(),
                              shiny::actionButton("run_analysis", "Run Analysis", class = "btn-info", width = "100%")
                            ),
                            
                            # Visualisation Controls
                            shiny::conditionalPanel(
                              condition = "input.interp_mode == 'vis'",
                              shiny::checkboxInput("show_cond_eff", "Conditional Effects", value = TRUE),
                              shiny::checkboxInput("show_prior_post", "Prior vs Posterior", value = FALSE),
                              shiny::br(),
                              shiny::actionButton("refresh_plots", "Update Plots", class = "btn-default", width = "100%")
                            )
                          )
            ),
            
            # Main Content Column
            shiny::column(8,
                          shiny::div(style = "padding-left: 10px;",
                                     shiny::h4("Results"),
                                     shiny::verbatimTextOutput("analysis_text"),
                                     shiny::plotOutput("analysis_plot_render", height = "400px")
                          )
            )
          )
        )
      )
    )
  )
  
  # --- Server Logic ---
  server <- function(input, output, session) {
    
    rv <- shiny::reactiveValues(
      data = NULL, data_load_code = NULL,
      fitted_models = list(), best_model_obj = NULL,
      priors_code = NULL, formula_str = NULL, code_snippets = character(0),
      
      # Store outputs here to ensure renderPlot always has something to read
      curr_plot = NULL,
      curr_print = "Run an analysis or visualisation to see results."
    )
    
    # 1. Data Loading 
    shiny::observe({
      dfs <- ls(envir = .GlobalEnv)
      valid_dfs <- character(0)
      if (length(dfs) > 0) {
        is_df <- vapply(dfs, function(x) {
          obj <- try(get(x, envir = .GlobalEnv), silent = TRUE)
          is.data.frame(obj) && !inherits(obj, "try-error")
        }, logical(1))
        valid_dfs <- dfs[is_df]
      }
      shiny::updateSelectInput(session, "dataset_env", choices = c("", valid_dfs))
    })
    
    shiny::observeEvent(input$dataset_env, {
      shiny::req(input$data_source == "env", input$dataset_env)
      if (input$dataset_env == "") return()
      rv$data <- get(input$dataset_env, envir = .GlobalEnv)
      rv$data_load_code <- paste0("data <- ", input$dataset_env)
      update_var_choices(rv$data)
    })
    
    shiny::observeEvent(input$btn_browse, {
      path <- rstudioapi::selectFile(caption = "Select Dataset", filter = "Data Files (*.csv *.xlsx *.rds *.dta *.sav)", existing = TRUE)
      if (!is.null(path) && path != "") {
        output$file_path_display <- shiny::renderText(basename(path))
        ext <- tolower(tools::file_ext(path))
        tryCatch({
          if (ext == "csv") { rv$data <- utils::read.csv(path); rv$data_load_code <- sprintf("data <- utils::read.csv(\"%s\")", path) }
          else if (ext == "rds") { rv$data <- readRDS(path); rv$data_load_code <- sprintf("data <- readRDS(\"%s\")", path) }
          else if (ext == "xlsx") { rv$data <- readxl::read_excel(path); rv$data_load_code <- sprintf("data <- readxl::read_excel(\"%s\")", path) }
          else if (ext %in% c("dta", "sav")) { 
            if(ext=="dta") { rv$data <- haven::read_dta(path); rv$data_load_code <- sprintf("data <- haven::read_dta(\"%s\")", path) }
            else { rv$data <- haven::read_spss(path); rv$data_load_code <- sprintf("data <- haven::read_spss(\"%s\")", path) }
          }
          update_var_choices(rv$data)
        }, error = function(e) shiny::showNotification(paste("Load failed:", e$message), type = "error"))
      }
    })
    
    update_var_choices <- function(df) {
      vars <- names(df)
      is_group <- sapply(df, function(x) (is.factor(x) || is.character(x)) && length(unique(x)) < nrow(df)/2)
      shiny::updateSelectInput(session, "response", choices = vars)
      shiny::updateSelectInput(session, "predictors", choices = vars)
      shiny::updateSelectInput(session, "grouping", choices = names(df)[is_group])
    }
    
    # 2. Visualisation & Advice
    output$response_dist_plot <- shiny::renderPlot({
      shiny::req(rv$data, input$response)
      y <- rv$data[[input$response]]
      df_plot <- data.frame(Value = as.numeric(y))
      
      ggplot2::ggplot(df_plot, ggplot2::aes(x = Value)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), fill = "#e3f2fd", color = "white", bins = 30) +
        ggplot2::geom_density(alpha = 0.2, fill = "#2196F3") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Distribution of", input$response), x = input$response, y = "Density")
    })
    
    output$data_advice <- shiny::renderUI({
      shiny::req(rv$data, input$response)
      y <- rv$data[[input$response]]
      
      is_ord <- is.ordered(y); is_cont <- is.numeric(y)
      is_int <- is_cont && all(y %% 1 == 0); is_pos <- is_cont && all(y > 0)
      is_prop <- is_cont && all(y >= 0 & y <= 1); is_bin <- length(unique(y)) == 2
      unique_vals <- length(unique(y)); pot_ord <- is_int && unique_vals > 2 && unique_vals < 15
      
      advice <- list(); fams <- c()
      if (is_ord) { advice <- c(advice, "Ordered factor detected. Recommended: Cumulative."); fams <- c("cumulative") }
      else if (is_bin) { advice <- c(advice, "Binary data detected. Recommended: Binomial."); fams <- c("binomial") }
      else if (is_prop) { advice <- c(advice, "Proportion (0-1). Recommended: Beta or Simplex."); fams <- c("Beta", "simplex") }
      else if (pot_ord) { advice <- c(advice, "Few integer levels. Recommended: Cumulative or Poisson."); fams <- c("cumulative", "poisson", "negbinomial") }
      else if (is_int && is_pos) { advice <- c(advice, "Count data. Recommended: Poisson or NegBin."); fams <- c("poisson", "negbinomial") }
      else if (is_pos) { 
        advice <- c(advice, "Positive continuous. Recommended: Gamma or Lognormal.")
        advice <- c(advice, "<div style='margin-top:5px;padding:8px;background:#fff3e0;border-left:3px solid #ff9800;font-size:0.9em'><b>Warning:</b> Do not log-transform outcome manually. Use Lognormal family on raw data for valid comparison.</div>")
        fams <- c("Gamma", "lognormal", "weibull") 
      } else { advice <- c(advice, "Continuous/Mixed. Recommended: Gaussian or Student-t."); fams <- c("gaussian", "student_t", "skew_normal") }
      
      shiny::updateSelectInput(session, "family_choice", choices = fams, selected = fams[1])
      if (length(input$grouping) > 0) advice <- c(advice, "Grouping variables selected. Fitting Mixed Effects model.")
      
      shiny::HTML(paste(advice, collapse = "<br/>"))
    })
    
    output$family_selector <- shiny::renderUI({ shiny::selectInput("family_choice", "Choose Family:", choices = c("gaussian"), multiple = TRUE) })
    output$family_description <- shiny::renderUI({
      shiny::req(input$family_choice)
      desc <- lapply(input$family_choice, function(f) { paste0("<b>", f, ":</b> ", tryCatch(family_info(f), error=function(e)"")) })
      shiny::HTML(paste(desc, collapse = "<br/>"))
    })
    
    # 3. Priors
    shiny::observeEvent(input$check_priors, {
      shiny::req(input$response, input$family_choice)
      f_fixed <- paste(input$response, "~", paste(c(1, input$predictors), collapse = " + "))
      rv$formula_str <- if (length(input$grouping)>0) paste(f_fixed, "+", paste(paste0("(1|", input$grouping, ")"), collapse=" + ")) else f_fixed
      
      tryCatch({
        pb <- prior_build_from_beliefs(stats::as.formula(rv$formula_str), rv$data, input$family_choice[1], 
                                       input$outcome_mean, input$outcome_sd, input$standardise)
        output$prior_plot <- shiny::renderPlot({
          invisible(utils::capture.output({
            fit <- qbrms(stats::as.formula(rv$formula_str), rv$data, family = input$family_choice[1], 
                         prior = pb$priors, sample_prior = "only", verbose = FALSE)
          }))
          pp_check(fit, type = "dens_overlay")
        })
        output$prior_summary <- shiny::renderPrint({ print(pb) })
        rv$priors_code <- prior_code(pb, object_name = "my_priors")
      }, error = function(e) shiny::showNotification(paste("Prior check failed:", e$message), type = "error"))
    })
    
    # 4. Fitting
    shiny::observeEvent(input$fit_model, {
      shiny::req(rv$formula_str, input$family_choice)
      rv$fitted_models <- list()
      
      withProgress(message = 'Fitting Models...', value = 0, {
        for (i in seq_along(input$family_choice)) {
          fam <- input$family_choice[i]; incProgress(1/length(input$family_choice), detail = paste("Fitting", fam))
          
          try({
            invisible(utils::capture.output({
              if (fam == "cumulative") {
                fit <- qbrmO(stats::as.formula(rv$formula_str), rv$data, family = cumulative(), verbose = FALSE)
              } else {
                fam_fun <- if(exists(fam, where=asNamespace("qbrms"))) get(fam, envir=asNamespace("qbrms")) else match.fun(fam)
                fit <- qbrms(stats::as.formula(rv$formula_str), rv$data, family = fam_fun(), prior = NULL, verbose = FALSE)
              }
            }))
            rv$fitted_models[[fam]] <- fit
          })
        }
      })
      
      output$fit_status <- shiny::renderPrint({ cat("Successfully fitted:", paste(names(rv$fitted_models), collapse = ", ")) })
      
      if (length(rv$fitted_models) > 1) {
        cmp <- loo_compare(rv$fitted_models)
        output$comparison_output <- shiny::renderPrint({ print(cmp) })
        output$comparison_interpretation <- shiny::renderUI({
          best_model <- rownames(cmp)[1]
          shiny::HTML(paste0(
            "<b>Best Model: ", best_model, "</b><br/>",
            "ELPD (Expected Log Predictive Density): <b>Higher is better</b>.<br/>",
            "LOOIC (Information Criteria): <b>Lower is better</b>.<br/>",
            "The <code>elpd_diff</code> shows the difference from the best model."
          ))
        })
        rv$best_model_obj <- rv$fitted_models[[1]] 
      } else {
        rv$best_model_obj <- rv$fitted_models[[1]]
        output$comparison_output <- shiny::renderPrint({ print(loo(rv$fitted_models[[1]])) })
        output$comparison_interpretation <- shiny::renderUI({ shiny::HTML("Only one model fitted. Check LOOIC (lower is better).") })
      }
      output$pp_check_plot <- shiny::renderPlot({ pp_check(rv$best_model_obj) })
      shiny::updateSelectInput(session, "emm_term", choices = input$predictors)
    })
    
    # 5. Interpretation Logic
    output$analysis_text <- shiny::renderPrint({ cat(rv$curr_print) })
    output$analysis_plot_render <- shiny::renderPlot({ 
      if(!is.null(rv$curr_plot)) print(rv$curr_plot) 
    })
    
    shiny::observeEvent(input$run_analysis, {
      shiny::req(rv$best_model_obj)
      type <- input$post_analysis_type
      
      tryCatch({
        if (type == "Estimated Marginal Means (EMMs)") {
          shiny::req(input$emm_term)
          em <- qbrms_emmeans(rv$best_model_obj, specs = input$emm_term)
          rv$curr_print <- paste(capture.output(print(em)), collapse="\n")
          # EMMs don't have a native plot, show instructions
          rv$curr_plot <- NULL 
          rv$code_snippets <- c(rv$code_snippets, paste0("qbrms_emmeans(fit, specs = '", input$emm_term, "')"))
          
        } else if (type == "Probability of Direction (pd)") {
          pd <- p_direction(rv$best_model_obj)
          rv$curr_print <- paste(capture.output(print(pd)), collapse="\n")
          rv$curr_plot <- plot_parameters(rv$best_model_obj) # Visual proxy
          rv$code_snippets <- c(rv$code_snippets, "p_direction(fit)")
          
        } else if (type == "Practical Significance") {
          ps <- p_significance(rv$best_model_obj, threshold = input$rope_range)
          rv$curr_print <- paste(capture.output(print(ps)), collapse="\n")
          rv$curr_plot <- plot(ps)
          rv$code_snippets <- c(rv$code_snippets, paste0("p_significance(fit, threshold = ", input$rope_range, ")"))
          
        } else if (type == "ROPE Analysis") {
          ro <- rope_analysis(rv$best_model_obj, rope = c(-input$rope_range, input$rope_range))
          rv$curr_print <- paste(capture.output(print(ro)), collapse="\n")
          rv$curr_plot <- plot_parameters(rv$best_model_obj, show_prior = FALSE) # Visual proxy
          rv$code_snippets <- c(rv$code_snippets, paste0("rope_analysis(fit, rope = c(-", input$rope_range, ", ", input$rope_range, "))"))
        }
      }, error = function(e) { rv$curr_print <- paste("Analysis failed:", e$message) })
    })
    
    shiny::observeEvent(input$refresh_plots, {
      shiny::req(rv$best_model_obj)
      rv$curr_print <- "Generating plots..."
      
      plots <- list()
      if (input$show_cond_eff) { 
        try({ 
          ce <- conditional_effects(rv$best_model_obj)
          plots[["cond"]] <- plot(ce)
          rv$code_snippets <- c(rv$code_snippets, "conditional_effects(fit)") 
        }) 
      }
      if (input$show_prior_post) { 
        try({ 
          pp <- plot_parameters(rv$best_model_obj, show_prior = TRUE)
          plots[["pp"]] <- pp
          rv$code_snippets <- c(rv$code_snippets, "plot_parameters(fit, show_prior = TRUE)") 
        }) 
      }
      
      if (length(plots) == 1) {
        rv$curr_plot <- plots[[1]]
      } else if (length(plots) > 1 && requireNamespace("gridExtra", quietly=TRUE)) {
        rv$curr_plot <- gridExtra::arrangeGrob(grobs = plots, ncol = 1)
      } else if (length(plots) > 0) {
        rv$curr_plot <- plots[[1]] # Fallback
      }
      rv$curr_print <- "Plots updated."
    })
    
    # 6. Code Insertion
    shiny::observeEvent(input$done, {
      code <- c("# --- qbrms Bayesian Workflow ---", "library(qbrms)", "", "# 1. Load Data", rv$data_load_code %||% "# No data selected", "", "# 2. Priors", rv$priors_code %||% "# Default priors", "", "# 3. Fit", paste0("fit <- qbrms(formula = ", rv$formula_str, ", data = data, family = ", input$family_choice[1], "(),", if (!is.null(rv$priors_code)) " prior = my_priors" else " prior = NULL", ")"), "", "# 4. Diagnostics", "pp_check(fit)", "summary(fit)")
      if (length(rv$code_snippets) > 0) code <- c(code, "", "# 5. Post-Hoc", unique(rv$code_snippets))
      rstudioapi::insertText(paste(code, collapse = "\n"))
      shiny::stopApp()
    })
  }
  viewer <- shiny::dialogViewer("Bayesian Workflow Coach", width = 1000, height = 850)
  shiny::runGadget(ui, server, viewer = viewer)
}