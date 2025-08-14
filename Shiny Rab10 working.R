library(shiny)
library(drc)
library(dplyr)
library(colourpicker)

### --- Fit Multiple Models and Select Best by RMSE ---
fit_best_models <- function(data) {
  # Remove vehicle control data (0 conc)
  data2 <- filter(data, Conc > 0)

# Define six model specifications using different LL.* functions

  models_specs <- list(
    list(name = "model1", formula = Rab10 ~ Log, fct = LL.2(fixed = c(NA, NA), names = c("Top", "IC50")), logDose = 10),
    list(name = "model2", formula = Rab10 ~ Log, fct = LL.3(fixed = c(NA, NA, NA), names = c("Bottom", "Top", "IC50")), logDose = 10),
    list(name = "model3", formula = Rab10 ~ Log, fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hillslope", "Bottom", "Top", "IC50")), logDose = 10),
    list(name = "model4", formula = Rab10 ~ Log, fct = LL.3(fixed = c(0, NA, NA), names = c("Bottom", "Top", "IC50")), logDose = 10),
    list(name = "model5", formula = Rab10 ~ Log, fct = LL.4(fixed = c(NA, 0, NA, NA), names = c("Hillslope", "Bottom", "Top", "IC50")), logDose = 10),
    list(name = "model6", formula = Rab10 ~ Log, fct = LL.4(fixed = c(1, 0, NA, NA), names = c("Hillslope", "Bottom", "Top", "IC50")), logDose = 10)
  )

  # Loop through each compound to fit models
  results <- list()
  for (cmpd in unique(data2$Compound)) {
    data_sub <- filter(data2, Compound == cmpd) # Subset data for the compound
    data_sub$Log <- log10(data_sub$Conc) # Add log-transformed concentration column

    # Loop through model specs and try fitting each model
    for (spec in models_specs) {
      model_name <- paste0(spec$name, "_", cmpd) # Create unique name for this model-compound combo
      fit <- try(drm(
        formula = spec$formula,
        data = data_sub,
        fct = spec$fct,
        logDose = spec$logDose
      ), silent = TRUE)

      if (inherits(fit, "try-error")) next # Skip if model fitting failed
      coefs <- coef(fit)
      ic50_val <- if ("IC50" %in% names(coefs)) coefs["IC50"] else coefs[grep("IC50", names(coefs))[1]] # Get IC50 value
      AIC_val <- AIC(fit)    # Calculate AIC for model
      predicted <- predict(fit) # Get predicted values
      observed <- data_sub$Rab10  # Actual observed Rab10 values
      residuals <- observed - predicted  # Calculate residuals
      rmse_val <- sqrt(mean(residuals^2, na.rm = TRUE)) # Root Mean Square Error as model fit metric

      # Store model object and metrics in the results list
      results[[model_name]] <- list(model = fit, Compound = cmpd, Model = spec$name, IC50 = ic50_val, AIC = AIC_val, RMSE = rmse_val)
    }
  }

  # Create a summary table of all models for all compounds
  summary_table <- bind_rows(
    lapply(results, function(res) {
      data.frame(Model = res$Model, Compound = res$Compound, IC50 = res$IC50, AIC = res$AIC, RMSE = res$RMSE, stringsAsFactors = FALSE)
    })
  )

  # For each compound, select the best model (lowest RMSE)
  best_models <- summary_table %>%
    group_by(Compound) %>%
    slice_min(order_by = RMSE, n = 1, with_ties = FALSE) %>%
    ungroup()

  best_fitted_models <- list()

  # Retrieve best model for each compound
  for (i in seq_len(nrow(best_models))) {
    cmpd <- best_models$Compound[i]
    model_name <- paste0(best_models$Model[i], "_", cmpd)
    if (model_name %in% names(results)) {
      best_fitted_models[[as.character(cmpd)]] <- results[[model_name]]$model
    }
  }

  # Return a list containing the full summary and best fitted models
  return(list(summary_table = summary_table, best_models = best_models, best_fitted_models = best_fitted_models))
}
##########################################################################################################################
### --- UI Definition ---
ui <- navbarPage("Rab10 dose-response ",

### --- Tab 1: Instructions for users ---
tabPanel("Instructions",
                          fluidPage(
                            h3("Instructions"),
                            p("1. Upload your CSV data file with columns: Compound, Conc (concentration in nM), Rab10 (response)."),
                            p("2. The app fits multiple dose-response models per compound and selects the best based on RMSE."),
                            p("3. Customize the plot colours, line widths, and toggle lines & labels on/off."),
                            p("4. View the dose-response curves and model parameter table."),
                            p("5. Use the download buttons to export plots or tables.")
                          )
                 ),

# --- Tab 2: Upload CSV Data File ---
tabPanel("Data Upload",
                          fluidPage(
                            fileInput("file1", "Upload CSV file",
                                      accept = c(".csv")), # File input control
                            checkboxInput("header", "Header", TRUE),  # Whether CSV has a header row
                            tags$hr(), # Horizontal rule for visual separation
                            h4("Preview of uploaded data"),# Header above preview table
                            tableOutput("contents")   # Displays first 20 rows of uploaded file

                          )
                 ),

# --- Tab 3: Plot and Customize Panel ---
 tabPanel("Plot & Customize",
                          sidebarLayout(
                            sidebarPanel(
                              checkboxInput("show_abline", "Show IC50 and Dmax lines", value = TRUE), # Toggle lines
                              hr(),
                              h4("Colours"),
                              colourInput("col_curve", "Curve colour", value = "blue"), # Colour of fitted curve
                              colourInput("col_ic50_v", "IC50 vertical line colour", value = "blue"),
                              colourInput("col_ic50_h", "IC50 horizontal line colour", value = "black"),
                              colourInput("col_dmax_obs", "Observed Dmax line colour", value = "red"),
                              colourInput("col_dmax_pred", "Predicted Dmax line colour", value = "darkgreen"),
                              hr(),
                              h4("Line widths"),
                              sliderInput("lw_curve", "Curve line width", min = 0.5, max = 5, value = 2, step = 0.1),
                              sliderInput("lw_lines", "Lines (IC50, Dmax) line width", min = 0.5, max = 5, value = 1.2, step = 0.1),
                              hr(),
                              h4("Download Options"),
                              numericInput("exportFont", "Font size (cex)", value = 0.9, min = 0.5, max = 5, step = 0.1),
                              numericInput("exportWidth", "Plot width (inches)", value = 5, min = 3, max = 10),
                              numericInput("exportHeight", "Plot height (inches)", value = 5, min = 3, max = 10),
                              downloadButton("downloadPlot", "Download Plot"), # Button to download PNG plot
                              br(), br(),
                              downloadButton("downloadTable", "Download Table")  # Button to download CSV table
                            ),

# --- Main Panel: Outputs ---
mainPanel(plotOutput("dosePlot", height = "600px"),# Dose-response curve plot
          hr(),
          h4("Model Parameters and Outputs"),
          dataTableOutput("model_table")  # Table showing IC50, AIC, Dmax etc.
                            )
                          )
                 )
)

### --- Server Logic ---

server <- function(input, output, session) {

# --- Reactive expression: read uploaded CSV file ---
  dataInput <- reactive({
    req(input$file1) # Ensure file is uploaded before proceeding
    df <- read.csv(input$file1$datapath, header = input$header)  # Read CSV with optional header
    # Validate minimal columns:
    validate(
      need(all(c("Compound", "Conc", "Rab10") %in% colnames(df)), "Data must have Compound, Conc, and Rab10 columns")
    )
    df$Compound <- as.factor(df$Compound)
    df
  })

  # --- Show a preview of uploaded data ---
  output$contents <- renderTable({
    req(dataInput())
    head(dataInput(), 20) # Display first 20 rows
  })

  # --- Fit models once file is uploaded ---
  fitResults <- reactive({
    req(dataInput())
    fit_best_models(dataInput())
  })

  # --- Plot the dose-response curves ---
  output$dosePlot <- renderPlot({
    req(fitResults())
    df <- dataInput()
    best_models <- fitResults()$best_models
    best_fitted_models <- fitResults()$best_fitted_models

    # Prepare plot limits and ticks
    df2 <- filter(df, Conc > 0)
    x_min <- min(df2$Conc)
    x_max <- max(df2$Conc)
    xlim_extended <- c(x_min / 10, x_max * 10)

    compounds <- names(best_fitted_models)  # ✅ Moved here BEFORE the loop
    pch_vals <- rep(19:25, length.out = length(compounds))

    for (i in seq_along(compounds)) {
      cmpd <- compounds[i]
      model <- best_fitted_models[[cmpd]]
      cmpd_data <- filter(df2, Compound == cmpd)

      # Prediction curve
      curve_range <- exp(seq(log(min(cmpd_data$Conc)), log(max(cmpd_data$Conc)), length.out = 200))
      pred_data <- data.frame(Log = log10(curve_range))
      curve_vals <- predict(model, newdata = pred_data)

      # Plot raw data and fitted model
      plot(NA, NA, log = "x",
           xlim = xlim_extended,
           ylim = c(0, 1.2),
           xaxs = "i",
           yaxs = "i",
           xaxt = "n",
           xlab = "Concentration (nM)",
           ylab = "pT73 Rab10:Rab10 (normalised to DMSO)",
           main = cmpd,
           cex.main = 1.2, font.main = 2, bty = "l", las = 1)

      log_range <- floor(log10(xlim_extended[1])):ceiling(log10(xlim_extended[2]))
      x_ticks <- 10^log_range
      axis(1, at = x_ticks, labels = sub("\\.0+$", "", format(x_ticks, trim = TRUE, scientific = FALSE)))

      compounds <- names(best_fitted_models)
      pch_vals <- rep(19:25, length.out = length(compounds))

      lines(curve_range, curve_vals, col = input$col_curve, lwd = input$lw_curve)

      # Mean & SEM points
      cmpd_summary <- cmpd_data %>%
        group_by(Conc) %>%
        summarise(
          Mean = mean(Rab10, na.rm = TRUE),
          SEM = sd(Rab10, na.rm = TRUE) / sqrt(n()),
          .groups = "drop"
        )

      points(cmpd_data$Conc, cmpd_data$Rab10, col = input$col_curve, pch = pch_vals[i], cex = 1.3)


      coefs <- coef(model)
      top <- if ("Top:(Intercept)" %in% names(coefs)) coefs["Top:(Intercept)"] else max(curve_vals, na.rm = TRUE)
      bottom <- if ("Bottom:(Intercept)" %in% names(coefs)) coefs["Bottom:(Intercept)"] else min(curve_vals, na.rm = TRUE)

      ic50_log <- if ("IC50:(Intercept)" %in% names(coefs)) log10(coefs["IC50:(Intercept)"]) else coefs[grep("IC50", names(coefs))[1]]
      ic50_conc <- 10^ic50_log
      ic50_response <- (top + bottom) / 2

      # Draw vertical and horizontal lines for IC50 if toggled on
      if (input$show_abline) {
        # IC50 lines
        abline(v = ic50_conc, col = input$col_ic50_v, lty = "dashed", lwd = input$lw_lines)
        abline(h = ic50_response, col = input$col_ic50_h, lty = "dashed", lwd = input$lw_lines)

        # IC50 label near vertical dashed line, just above x-axis
        text(x = ic50_conc * 1.1, y = par("usr")[4] * 0.95,
             labels = bquote(IC[50] ~ "=" ~ .(signif(ic50_conc, 3)) ~ "nM"),
             col = input$col_curve, cex = 1, pos = 4)

        # 50% max inhibition label below horizontal line
        text(x = min(xlim_extended), y = ic50_response - 0.05,
             labels = "50% of maximum inhibition", pos = 4, cex = 0.7, col = input$col_ic50_h)
      }

      max_conc <- max(cmpd_summary$Conc)
      # Mark observed Dmax
      dmax_obs <- cmpd_summary$Mean[cmpd_summary$Conc == max_conc]
      if (input$show_abline) {
        abline(h = dmax_obs, col = input$col_dmax_obs, lty = "dashed", lwd = input$lw_lines)
      }

      # Mark predicted Dmax (bottom parameter)
      bottom_threshold <- 0.02
      if (abs(dmax_obs - bottom) > bottom_threshold && bottom <= dmax_obs) {
        if (input$show_abline) {
          abline(h = bottom, col = input$col_dmax_pred, lty = "dashed", lwd = input$lw_lines)
        }
      }
    }

    # Legend for Dmax observed and predicted
    legend_labels2 <- c()
    legend_colors2 <- c()

    for (cmpd in compounds) {
      model <- best_fitted_models[[cmpd]]
      cmpd_data <- filter(df2, Compound == cmpd)
      cmpd_summary <- cmpd_data %>%
        group_by(Conc) %>%
        summarise(Mean = mean(Rab10, na.rm = TRUE), .groups = "drop")
      max_conc <- max(cmpd_summary$Conc)
      dmax_obs_val <- cmpd_summary$Mean[cmpd_summary$Conc == max_conc]

      coefs <- coef(model)
      bottom <- if ("Bottom:(Intercept)" %in% names(coefs)) coefs["Bottom:(Intercept)"] else min(predict(model), na.rm = TRUE)

      perc_deg_obs <- (1 - dmax_obs_val) * 100
      pred_conc_100 <- if (perc_deg_obs > 0) round((100 * max_conc) / perc_deg_obs) else NA

      legend_labels2 <- c(
        legend_labels2,
        paste0("Observed Dmax = ", round(dmax_obs_val, 2),
               " (", round(perc_deg_obs), "% at ", max_conc, " nM)"),
        if (!is.na(pred_conc_100) && abs(dmax_obs_val - bottom) > 0.02 && bottom <= dmax_obs_val)
          {paste0("Predicted Dmax = ", round(bottom, 2),
                 " (100% at ", pred_conc_100, " nM)")
        } else {
          NULL
        }
      )

      legend_colors2 <- c(legend_colors2,
                          input$col_dmax_obs,
                          if (!is.na(pred_conc_100) && abs(dmax_obs_val - bottom) > 0.02) input$col_dmax_pred else NULL)
    }

    par(xpd = TRUE)
    legend("bottomleft",
           legend = legend_labels2,
           text.col = legend_colors2,
           bty = "n",
           cex = 0.7,
           x.intersp = 0.5,
           y.intersp = 1.1)
    par(xpd = FALSE)
  })

  # Table output: Model parameters summary + dmax info
  output$model_table <- renderDataTable({
    req(fitResults())
    summary_table <- fitResults()$summary_table
    best_models <- fitResults()$best_models
    best_fitted_models <- fitResults()$best_fitted_models
    df <- dataInput()
    df2 <- filter(df, Conc > 0)

    # For each compound get observed and predicted Dmax:
    dmax_info <- lapply(names(best_fitted_models), function(cmpd) {
      model <- best_fitted_models[[cmpd]]
      cmpd_data <- filter(df2, Compound == cmpd)
      cmpd_summary <- cmpd_data %>%
        group_by(Conc) %>%
        summarise(Mean = mean(Rab10, na.rm = TRUE), .groups = "drop")
      max_conc <- max(cmpd_summary$Conc)
      dmax_obs_val <- cmpd_summary$Mean[cmpd_summary$Conc == max_conc]

      coefs <- coef(model)
      bottom <- if ("Bottom:(Intercept)" %in% names(coefs)) coefs["Bottom:(Intercept)"] else min(predict(model), na.rm = TRUE)

      perc_deg_obs <- (1 - dmax_obs_val) * 100
      pred_conc_100 <- if (perc_deg_obs > 0) round((100 * max_conc) / perc_deg_obs) else NA

      data.frame(
        Compound = cmpd,
        Observed_Dmax = round(dmax_obs_val, 3),
        Observed_Dmax_pct = round(perc_deg_obs, 1),
        Observed_Conc = max_conc,
        Predicted_Dmax = if (bottom <= dmax_obs_val) round(bottom, 3) else NA,
        Predicted_Conc_100pct = if (bottom <= dmax_obs_val) pred_conc_100 else NA
      )
    }) %>% bind_rows()

    # Join with best_models summary
    final_tbl <- best_models %>%
      left_join(dmax_info, by = "Compound") %>%
      select(Compound, Model, IC50, AIC, RMSE,
             Observed_Dmax, Observed_Dmax_pct, Observed_Conc,
             Predicted_Dmax, Predicted_Conc_100pct)

    final_tbl
  }, options = list(pageLength = 10, scrollX = TRUE))

  # --- Download Plot as PNG ---
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("dose_response_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      png(file,
          width = input$exportWidth,
          height = input$exportHeight,
          units = "in",
          res = 300)
      par(cex = input$exportFont)  # This applies the larger font scaling      # Recreate the plot
      df <- dataInput()
      best_models <- fitResults()$best_models
      best_fitted_models <- fitResults()$best_fitted_models

      df2 <- filter(df, Conc > 0)
      x_min <- min(df2$Conc)
      x_max <- max(df2$Conc)
      xlim_extended <- c(x_min / 10, x_max * 10)

      plot(NA, NA, log = "x",
           xlim = xlim_extended,
           ylim = c(0, 1.2),
           xaxs = "i",
           yaxs = "i",
           xaxt = "n",
           xlab = "Concentration (nM)",
           ylab = "pT73 Rab10:Rab10 (normalised to DMSO)",
           main = cmpd,
           cex.main = 1.2, font.main = 2, bty = "l", las = 1)

      log_range <- floor(log10(xlim_extended[1])):ceiling(log10(xlim_extended[2]))
      x_ticks <- 10^log_range
      axis(1, at = x_ticks, labels = sub("\\.0+$", "", format(x_ticks, trim = TRUE, scientific = FALSE)))

      compounds <- names(best_fitted_models)
      pch_vals <- rep(19:25, length.out = length(compounds))

      for (i in seq_along(compounds)) {
        cmpd <- compounds[i]
        model <- best_fitted_models[[cmpd]]
        cmpd_data <- filter(df2, Compound == cmpd)

        curve_range <- exp(seq(log(min(cmpd_data$Conc)), log(max(cmpd_data$Conc)), length.out = 200))
        pred_data <- data.frame(Log = log10(curve_range))
        curve_vals <- predict(model, newdata = pred_data)

        lines(curve_range, curve_vals, col = input$col_curve, lwd = input$lw_curve)

        cmpd_summary <- cmpd_data %>%
          group_by(Conc) %>%
          summarise(
            Mean = mean(Rab10, na.rm = TRUE),
            SEM = sd(Rab10, na.rm = TRUE) / sqrt(n()),
            .groups = "drop"
          )

        points(cmpd_data$Conc, cmpd_data$Rab10, col = input$col_curve, pch = pch_vals[i], cex = 1.3)

        coefs <- coef(model)
        top <- if ("Top:(Intercept)" %in% names(coefs)) coefs["Top:(Intercept)"] else max(curve_vals, na.rm = TRUE)
        bottom <- if ("Bottom:(Intercept)" %in% names(coefs)) coefs["Bottom:(Intercept)"] else min(curve_vals, na.rm = TRUE)

        ic50_log <- if ("IC50:(Intercept)" %in% names(coefs)) log10(coefs["IC50:(Intercept)"]) else coefs[grep("IC50", names(coefs))[1]]
        ic50_conc <- 10^ic50_log
        ic50_response <- (top + bottom) / 2

        if (input$show_abline) {
          abline(v = ic50_conc, col = input$col_ic50_v, lty = "dashed", lwd = input$lw_lines)
          abline(h = ic50_response, col = input$col_ic50_h, lty = "dashed", lwd = input$lw_lines)

          text(x = ic50_conc * 1.1, y = par("usr")[4] * 0.95,
               labels = bquote(IC[50] ~ "=" ~ .(signif(ic50_conc, 3)) ~ "nM"),
               col = input$col_curve, cex = 1, pos = 4)

          text(x = min(xlim_extended), y = ic50_response - 0.05,
               labels = "50% of maximum inhibition", pos = 4, cex = 0.7, col = input$col_ic50_h)
        }

        max_conc <- max(cmpd_summary$Conc)
        dmax_obs <- cmpd_summary$Mean[cmpd_summary$Conc == max_conc]
        if (input$show_abline) {
          abline(h = dmax_obs, col = input$col_dmax_obs, lty = "dashed", lwd = input$lw_lines)
        }

        bottom_threshold <- 0.02
        if (abs(dmax_obs - bottom) > bottom_threshold) {
          if (input$show_abline) {
            abline(h = bottom, col = input$col_dmax_pred, lty = "dashed", lwd = input$lw_lines)
          }
        }
      }

      legend_labels2 <- c()
      legend_colors2 <- c()

      for (cmpd in compounds) {
        model <- best_fitted_models[[cmpd]]
        cmpd_data <- filter(df2, Compound == cmpd)
        cmpd_summary <- cmpd_data %>%
          group_by(Conc) %>%
          summarise(Mean = mean(Rab10, na.rm = TRUE), .groups = "drop")
        max_conc <- max(cmpd_summary$Conc)
        dmax_obs_val <- cmpd_summary$Mean[cmpd_summary$Conc == max_conc]

        coefs <- coef(model)
        bottom <- if ("Bottom:(Intercept)" %in% names(coefs)) coefs["Bottom:(Intercept)"] else min(predict(model), na.rm = TRUE)

        perc_deg_obs <- (1 - dmax_obs_val) * 100
        pred_conc_100 <- if (perc_deg_obs > 0) round((100 * max_conc) / perc_deg_obs) else NA

        legend_labels2 <- c(
          legend_labels2,
          paste0("Observed Dmax = ", round(dmax_obs_val, 2),
                 " (", round(perc_deg_obs), "% at ", max_conc, " nM)"),
          if (!is.na(pred_conc_100) && abs(dmax_obs_val - bottom) > 0.02) {
            paste0("Predicted Dmax = ", round(bottom, 2),
                   " (100% at ", pred_conc_100, " nM)")
          } else {
            NULL
          }
        )

        legend_colors2 <- c(legend_colors2,
                            input$col_dmax_obs,
                            if (!is.na(pred_conc_100) && abs(dmax_obs_val - bottom) > 0.02) input$col_dmax_pred else NULL)
      }

      par(xpd = TRUE)
      legend("bottomleft",
             legend = legend_labels2,
             text.col = legend_colors2,
             bty = "n",
             cex = 0.7,
             x.intersp = 0.5,
             y.intersp = 1.1)
      par(xpd = FALSE)
      dev.off()
    }
  )

  # Download table handler
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("dose_response_table_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(fitResults())
      summary_table <- fitResults()$summary_table
      best_models <- fitResults()$best_models
      best_fitted_models <- fitResults()$best_fitted_models
      df <- dataInput()
      df2 <- filter(df, Conc > 0)

      dmax_info <- lapply(names(best_fitted_models), function(cmpd) {
        model <- best_fitted_models[[cmpd]]
        cmpd_data <- filter(df2, Compound == cmpd)
        cmpd_summary <- cmpd_data %>%
          group_by(Conc) %>%
          summarise(Mean = mean(Rab10, na.rm = TRUE), .groups = "drop")
        max_conc <- max(cmpd_summary$Conc)
        dmax_obs_val <- cmpd_summary$Mean[cmpd_summary$Conc == max_conc]

        coefs <- coef(model)
        bottom <- if ("Bottom:(Intercept)" %in% names(coefs)) coefs["Bottom:(Intercept)"] else min(predict(model), na.rm = TRUE)

        perc_deg_obs <- (1 - dmax_obs_val) * 100
        pred_conc_100 <- if (perc_deg_obs > 0) round((100 * max_conc) / perc_deg_obs) else NA

        data.frame(
          Compound = cmpd,
          Observed_Dmax = round(dmax_obs_val, 3),
          Observed_Dmax_pct = round(perc_deg_obs, 1),
          Observed_Conc = max_conc,
          Predicted_Dmax = round(bottom, 3),
          Predicted_Conc_100pct = pred_conc_100
        )
      }) %>% bind_rows()

      final_tbl <- best_models %>%
        left_join(dmax_info, by = "Compound") %>%
        select(Compound, Model, IC50, AIC, RMSE,
               Observed_Dmax, Observed_Dmax_pct, Observed_Conc,
               Predicted_Dmax, Predicted_Conc_100pct)

      write.csv(final_tbl, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
