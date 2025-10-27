library(shiny)
library(bslib)
library(dplyr)
library(reactable)
library(bsicons)
library(sparkline)
library(writexl) # 

source('global.R')
source('R/chromatogram.R')
source('R/sequence-viewer.R')

sidebar <- sidebar(
  tags$div(
    style = "padding: 10px;",
    tags$h4("Upload Sanger AB1 Files", style = "margin-bottom: 15px;"),
    fileInput(
      "ab1",
      label = NULL,
      accept = ".ab1",
      multiple = TRUE,
      buttonLabel = "Browse...",
      placeholder = "No files selected"
    ),
    #tags$hr(),
    tags$p(
      "Select one or more .ab1 files to view summary statistics and traces.",
      style = "color: #888; font-size: 13px;"
    ),
    #checkboxInput('shorten_names', 'Shorten names', value = F),
    tags$br(),
    actionButton('settings', 'QC settings', style = "margin-bottom: 15px;"),
    tags$br(),
    tags$hr(),
    tags$span(
      bsicons::bs_icon("info-circle"),
      HTML(" Supported:<br>Sanger AB1 files"),
      style = "color: #888; font-size: 13px;"
    )
  )
)

ui <- page_navbar(
  sidebar = sidebar,
  title = "",
  header = tags$head(
    tags$style(
      # CSS for chromatogram plot container and quality colors
      HTML("
        .scrollable-plot-container {
          overflow-x: auto;
          width: 100%;
          white-space: nowrap;
        }
        /* Quality score colors */
        .qv-high { color: #4CAF50; }
        .qv-medium { color: #FFC107; }
        .qv-low { color: #F44336; }
        /* Chromatogram container styles */
        .chromatogram-output {
          border: 1px solid #eee;
          border-radius: 4px;
          background: white;
          margin: 8px 0;
        }
      ")
    ),
    # Add D3.js dependency
    tags$script(src = "https://d3js.org/d3.v7.min.js"),
    # Add our visualization modules
    tags$script(src = "js/chromatogram.js"),
    tags$script(src = "js/sequence-viewer.js")
  ),
  nav_panel('Basecall',
            reactableOutput('table1'),
            # Use uiOutput to conditionally render the download button and footer
            uiOutput('qc_controls_and_footer') 
  ),
  nav_panel('Raw signal',
            reactableOutput('table2')
  ),
  nav_panel('Sequence', 
            reactableOutput('table3')
  ),
  nav_panel('CRL plots',
            reactableOutput('table4')        
  ),
  #nav_panel('QC Summary'),
  nav_panel('Run info', 
            reactableOutput('table5')
  ),
  theme = bs_theme(bootswatch = "simplex")
)

##############################
# Store default QC thresholds 
default_qc <- data.frame(
  crl_window_size = 20,
  crl_qv_threshold = 20,
  crl20_fail = 500,
  crl20_suspect = 800,
  basesQ20_fail = 500,
  basesQ20_suspect = 800,
  sn_fail = 60,
  sn_suspect = 180
)
##############################

server <- function(input, output, session) {
  
  # Show modal when QC settings button is clicked
  observeEvent(input$settings, {
    showModal(
      modalDialog(
        size = 'l',
        title = "QC Settings",
        tags$p('CRL settings', style = "font-weight: bold; font-size: 15px;"),
        
        tags$div(
          style = "display: flex; gap: 16px; width: 100%; font-size: 12px; color: #1C398E;",
          sliderInput(
            "crl_window_size",
            "CRL window size",
            min = 5,
            max = 100,
            value = max(5, min(100, qc_thresholds$crl_window_size)),
            width = "100%"
          ),
          sliderInput(
            "crl_qv_threshold",
            "CRL QV threshold",
            min = 5,
            max = 50,
            value = max(5, min(50, qc_thresholds$crl_qv_threshold)),
            width = "100%"
          )
        ),
        tags$hr(style = "border-top: 1px solid #D6D3D1; margin: 16px 0;"),
        
        tags$p("Thresholds for sample QC flagging:", style = "font-weight: bold; font-size: 15px;"),
        tags$div(
          style = "width: 100%; font-size: 12px; color: #1C398E;",
          sliderInput(
            "qc_crl20",
            "CRL thresholds (fail, suspect)",
            min = 0,
            max = 1500,
            value = c(qc_thresholds$crl20_fail, qc_thresholds$crl20_suspect),
            step = 10,
            width = "100%"
          )
        ),
        tags$div(
          style = "width: 100%; font-size: 12px; color: #1C398E;",
          sliderInput(
            "qc_basesQ20",
            "Q20+ thresholds (fail, suspect)",
            min = 0,
            max = 1500,
            value = c(
              qc_thresholds$basesQ20_fail,
              qc_thresholds$basesQ20_suspect
            ),
            step = 10,
            width = "100%"
          )
        ),
        tags$div(
          style = "width: 100%; font-size: 12px; color: #1C398E;",
          sliderInput(
            "qc_sn",
            "Signal/Noice (fail, suspect)",
            min = 0,
            max = 1000,
            value = c(
              qc_thresholds$sn_fail,
              qc_thresholds$sn_suspect
            ),
            step = 1,
            width = "100%"
          )
        ),
        footer = tagList(
          actionButton('restore_qc', 'Restore defaults'),
          modalButton("Cancel"),
          actionButton("apply_qc_settings", "Apply")
        ),
        easyClose = TRUE
      )
    )
  })
  
  # Store QC thresholds in reactive values
  qc_thresholds <- reactiveValues(
    crl_window_size = 20,
    crl_qv_threshold = 20,
    crl20_fail = 500,
    crl20_suspect = 800,
    basesQ20_fail = 500,
    basesQ20_suspect = 800,
    sn_fail = 60,
    sn_suspect = 180
  )
  
  # Update thresholds when Apply is clicked
  observeEvent(input$apply_qc_settings, {
    qc_thresholds$crl_window_size <- max(5, input$crl_window_size)
    qc_thresholds$crl_qv_threshold <- max(5, input$crl_qv_threshold)
    qc_thresholds$crl20_fail <- input$qc_crl20[1]
    qc_thresholds$crl20_suspect <- input$qc_crl20[2]
    qc_thresholds$basesQ20_fail <- input$qc_basesQ20[1]
    qc_thresholds$basesQ20_suspect <- input$qc_basesQ20[2]
    qc_thresholds$sn_fail <- input$qc_sn[1]
    qc_thresholds$sn_suspect <- input$qc_sn[2]
    removeModal()
  })
  
  # Restore thresholds when Restore is clicked
  observeEvent(input$restore_qc, {
    slider_updates <- list(
      crl_window_size = default_qc$crl_window_size,
      crl_qv_threshold = default_qc$crl_qv_threshold,
      qc_crl20 = c(default_qc$crl20_fail, default_qc$crl20_suspect),
      qc_basesQ20 = c(default_qc$basesQ20_fail, default_qc$basesQ20_suspect),
      qc_sn = c(default_qc$sn_fail, default_qc$sn_suspect)
    )
    
    # Use lapply to iterate through the list and update each slider input
    lapply(names(slider_updates), function(input_id) {
      updateSliderInput(
        inputId = input_id,
        value = slider_updates[[input_id]],
        session = session
      )
    })
  })
  
  
  # Create a reusable list of QC col defs with conditional formatting, it is for df2() that is used 3 times
  create_qc_coldefs <- function(nsamples, crl_fail, crl_suspect, qv_fail, qv_suspect, sn_fail, sn_suspect) {
    list(
      # Sample name
      sample = colDef(
        minWidth = 200, 
        footer = paste0("Total ", nsamples, " samples")),
        #footer = paste0("Total ", "samples")),
      # Other
      rawSeqLen = colDef(minWidth = 50),
      well = colDef(minWidth = 30, show = FALSE),
      data = colDef(show = FALSE),
      crl_start = colDef(show = FALSE),
      crl_end = colDef(show = FALSE),
      # CRL
      crl20 = colDef(
        name = paste0("CRL",  qc_thresholds$crl_qv_threshold), 
        minWidth = 60,
        html = T,
        footer = function(values) {
          paste0("Min: ", round(min(values), 0), "<br>", "Max: ", round(max(values), 0), "<br>", "Mean: ", round(mean(values), 0))
        },
        style = function(value) {
          if (is.na(value)) return(list(color = "#eee"))
          if (value < crl_fail) return(list(color = "#F44336", fontWeight = "normal"))
          if (value < crl_suspect) return(list(color = "#FFC107", fontWeight = "normal"))
          list(color = "#4CAF50", fontWeight = "normal")
        }
      ),
      basesQ20 = colDef(
        name = "QV20+",
        minWidth = 60,
        html = T,
        footer = function(values) {
          paste0("Min: ", round(min(values), 0), "<br>", "Max: ", round(max(values), 0), "<br>", "Mean: ", round(mean(values), 0))
        },
        style = function(value) {
          if (is.na(value)) return(list(color = "#eee"))
          if (value < qv_fail) return(list(color = "#F44336", fontWeight = "normal"))
          if (value < qv_suspect) return(list(color = "#FFC107", fontWeight = "normal"))
          list(color = "#4CAF50", fontWeight = "normal")
        }
      ),
      meanSNratio = colDef(
        name = "Signal/Noice", 
        format = colFormat(digits = 0),
        minWidth = 60,
        html = T,
        footer = function(values) {
          paste0(
            "Min: ", round(min(values), 0), "<br>",
            "Max: ", round(max(values), 0), "<br>",
            "Mean: ", round(mean(values), 0))
        },
        style = function(value) {
          if (is.na(value)) return(list(color = "#eee"))
          if (value < sn_fail) return(list(color = "#F44336", fontWeight = "normal"))
          if (value < sn_suspect) return(list(color = "#FFC107", fontWeight = "normal"))
          list(color = "#4CAF50", fontWeight = "normal")
        }
      ),
      QC_flag = colDef(
        name = "QC flag",
        html = T, 
        align = 'right',
        minWidth = 60,
        footer = function(value){
          paste0(
            '<a style="color:#4CAF50;">Pass: </a><b>', str_count(str_flatten(value), 'pass'), 
            '</b><br><a style="color:#FFC107;">Suspect: </a><b>', str_count(str_flatten(value), 'suspect'), 
            '</b><br><a style="color:#F44336;">Fail: </a><b>', str_count(str_flatten(value), 'fail')
          )
        },
        cell = function(value) {
          color <- switch(
            value,
            "fail" = "#F44336",
            "suspect" = "#FFC107",
            "pass" = "#4CAF50",
            "#eee"
          )
          htmltools::tagList(
            value,
            htmltools::tags$span(
              style = paste0(
                "display: inline-block; width: 14px; height: 14px; border-radius: 50%; background:", 
                color, "; margin-left: 8px; margin-right: 2px; vertical-align: middle;"
              ),
              ""
            )
          )
        }
      )
    )
  }
  
  # Create column defs list
  # coldefs <- reactive({
  #   create_qc_coldefs(
  #     crl_fail = qc_thresholds$crl20_fail, crl_suspect = qc_thresholds$crl20_suspect, 
  #     qv_fail = qc_thresholds$basesQ20_fail, qv_suspect = qc_thresholds$basesQ20_suspect,
  #     sn_fail = qc_thresholds$sn_fail, sn_suspect = qc_thresholds$sn_suspect
  #   )
  # })
  # 
  df <- reactive({
    ab1list <- input$ab1$datapath
    req(ab1list)
    withProgress(message = "Processing AB1 files...", value = 0, {
      n <- length(ab1list)
      results <- vector("list", n)
      for (i in seq_along(ab1list)) {
        results[[i]] <- get_ab1(ab1list[i])
        incProgress(1/n, detail = paste("File", i, "of", n))
        #print(results[[i]]$sample)
      }
      bind_rows(results)
    })
  })
  
  df2 <- reactive({
    req(df()) 
    df() %>%
      rowwise() %>%
      mutate(
        crl_results = list(crl(
          qscores = unlist(qscores), window_size = qc_thresholds$crl_window_size, qval = qc_thresholds$crl_qv_threshold
        )),
        crl20 = crl_results$crl,
        crl_start = crl_results$crl_start,
        crl_end = crl_results$crl_end
      ) %>%
      select(-crl_results)
    # window size and QV threshold
  })
  
  # Reactive data frame for table1, which is also used for download
  df_table1_data <- reactive({
    req(df2()) # Ensure data exists
    df2() %>%
      select('sample', 'well', 'rawSeqLen', 'crl20', 'basesQ20', 'data', 'meanSNratio','crl_start', 'crl_end') %>%
      rowwise() %>%
      mutate(
        #sample = ifelse(input$shorten_names, str_trunc(sample, width = 42), sample),
        #rawMeanQscore = round(rawMeanQscore),
        QC_flag = case_when(
          crl20 < qc_thresholds$crl20_fail |
            meanSNratio < qc_thresholds$sn_fail |
              basesQ20 < qc_thresholds$basesQ20_fail ~ "fail",
          crl20 < qc_thresholds$crl20_suspect |
            meanSNratio < qc_thresholds$sn_suspect |
              basesQ20 < qc_thresholds$basesQ20_suspect ~ "suspect",
          TRUE ~ "pass"
        )
      ) %>%
      ungroup() # Convert back to a standard data frame
  })
  
  
  # Basecall chromatogram table 
  output$table1 <- renderReactable({
    req(df_table1_data()) 
    reactable(
      df_table1_data(),
      pagination = FALSE, searchable = TRUE, highlight = TRUE, bordered = TRUE, striped = FALSE, compact = TRUE, resizable = TRUE,
      style = list(fontSize = "14px"),
      defaultColDef = colDef(footerStyle = list(color='grey', fontWeight = 'normal')),
      onClick = "expand", # Expand row details on click
      columns = create_qc_coldefs(
        nsamples = nrow(df_table1_data()),
        crl_fail = qc_thresholds$crl20_fail, crl_suspect = qc_thresholds$crl20_suspect,
        qv_fail = qc_thresholds$basesQ20_fail, qv_suspect = qc_thresholds$basesQ20_suspect,
        sn_fail = qc_thresholds$sn_fail, sn_suspect = qc_thresholds$sn_suspect
      ),
      #columns = coldefs(),
      # ADD THIS rowStyle TO DRAW A LINE ABOVE THE FOOTER
      rowStyle = function(index) {
        if (index == nrow(df_table1_data())) {
          list(borderBottom = "1px solid grey") #
        }
      },
      #################
      details = function(index) {
        # Define unique IDs for the chromatogram
        chrom_output_id <- paste0('chrom_', index)
        detail_container_id <- paste0('detail_chrom', index)
        
        # Render the chromatogram when the row is expanded
        local({
          output[[chrom_output_id]] <- renderChromatogram({
            raw_abif_data <- df2()[index, ]$data
            
            
            if (length(raw_abif_data) > 0) {
              # Format data for D3
              list(
                traces = list(
                  DATA1 = raw_abif_data$data$DATA.9,
                  DATA2 = raw_abif_data$data$DATA.10,
                  DATA3 = raw_abif_data$data$DATA.11,
                  DATA4 = raw_abif_data$data$DATA.12
                ),
                bases = strsplit(raw_abif_data$data$PBAS.1, "")[[1]],
                # Convert PLOC.1 to zero-based indices for JS
                peakLocations = as.numeric(raw_abif_data$data$PLOC.1) - 1,
                qualityScores = raw_abif_data$data$PCON.1,
                plotWidth = 11000  # Set a wide plot width like the original ggplot
              )
            }
          })
        })
        
        # Return the container with our custom chromatogram output
        htmltools::div(
          id = detail_container_id,
          style = "padding: 10px; padding-top: 0px; padding-bottom: 10px;",
          htmltools::div(
            class = "scrollable-plot-container",
            style = "overflow-x: auto; width: 100%; white-space: nowrap;",
            chromatogramOutput(chrom_output_id, width = "11000px", height = "250px")
          )
        )
      }
    )
  })
  
  # Conditional rendering of the download button and the footer
  output$qc_controls_and_footer <- renderUI({
    # Ensure data is present before rendering UI elements that rely on it
    req(df_table1_data()) 
    
    # Use a main div with display: flex to put the two items side-by-side
    tags$div(
      style = "display: flex; justify-content: space-between; align-items: flex-start; margin-top: 10px;",
      
      # LEFT COLUMN:
     
      tags$div(
        style = "flex: 1 1 50%; padding-top: 5px;", 
        tags$details(
          #open = NA, # This makes the details expanded by default
          style = "margin-top: 0px;", # Removed margin-top since it's now handled by the parent div
          tags$summary(
            style = "font-size: 13px; color: #37474f; cursor: pointer;",
            "QC settings and thresholds used"
          ),
          tags$div(
            style = "border-radius: 6px; padding: 10px 16px; font-size: 12px; color: #37474f;",
            tags$ul(
              style = "margin: 8px 0 0 18px; padding: 0;",
              tags$li(
                paste0("CRL window size: ", qc_thresholds$crl_window_size)
              ),
              tags$li(
                paste0("CRL QV threshold: ", qc_thresholds$crl_qv_threshold)
              ),
              tags$li(
                paste0("CRL20: fail < ", qc_thresholds$crl20_fail, ", suspect < ", qc_thresholds$crl20_suspect)
              ),
              tags$li(
                paste0("Q20+: fail < ", qc_thresholds$basesQ20_fail, ", suspect < ", qc_thresholds$basesQ20_suspect)
              ),
              tags$li(
                paste0("Signal/Noice: fail < ", qc_thresholds$sn_fail, ", suspect < ", qc_thresholds$sn_suspect)
              )
            )
          )
        )
      ),
      # RIGHT COLUMN:
      tags$div(
        style = "display: flex; justify-content: flex-end;", 
        downloadButton("download_table1", "Download QC Data (Excel)", class = "btn-info btn-sm")
      )
    )
  })
  
  # DOWNLOAD HANDLER FOR TABLE1
  output$download_table1 <- downloadHandler(
    filename = function() {
      paste0("sanger_qc_data-", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      # Prepare data for download
      data_to_export <- df_table1_data() %>%
        # Rename columns 
        rename(
          !!paste0("CRL", qc_thresholds$crl_qv_threshold) := crl20,
          QV20_plus = basesQ20,
          Mean_Signal_Noise = meanSNratio,
          Raw_Seq_Length = rawSeqLen
        ) %>%
        # Select and order the final columns
        select(
          Sample = sample, 
          Raw_Seq_Length, 
          !!paste0("CRL", qc_thresholds$crl_qv_threshold), 
          QV20_plus, 
          Mean_Signal_Noise, 
          QC_flag
        )
      
      # Prepare metadata/settings to include in a separate sheet
      qc_settings_df <- data.frame(
        Setting = c(
          "CRL window size",
          "CRL QV threshold",
          paste0("CRL", qc_thresholds$crl_qv_threshold, " Fail Threshold"),
          paste0("CRL", qc_thresholds$crl_qv_threshold, " Suspect Threshold"),
          "QV20+ Fail Threshold",
          "QV20+ Suspect Threshold",
          "SN Ratio Fail Threshold",
          "SN Ratio Suspect Threshold"
        ),
        Value = c(
          qc_thresholds$crl_window_size,
          qc_thresholds$crl_qv_threshold,
          qc_thresholds$crl20_fail,
          qc_thresholds$crl20_suspect,
          qc_thresholds$basesQ20_fail,
          qc_thresholds$basesQ20_suspect,
          qc_thresholds$sn_fail,
          qc_thresholds$sn_suspect
        )
      )
      
      # Write both data and settings to a multi-sheet Excel file
      writexl::write_xlsx(
        list("QC_Summary" = data_to_export, "QC_Settings" = qc_settings_df), 
        path = file
      )
    }
  )
  
  # Raw signal
  output$table2 <- renderReactable({

    reactable(
      df_table1_data(),
      pagination = FALSE, searchable = TRUE, highlight = TRUE, bordered = TRUE, striped = FALSE, compact = TRUE, resizable = TRUE,
      style = list(fontSize = "14px"),
      defaultColDef = colDef(footerStyle = list(color='grey', fontWeight = 'normal')),
      onClick = "expand", # Expand row details on click
      columns = create_qc_coldefs(
        nsamples = nrow(df_table1_data()),
        crl_fail = qc_thresholds$crl20_fail, crl_suspect = qc_thresholds$crl20_suspect,
        qv_fail = qc_thresholds$basesQ20_fail, qv_suspect = qc_thresholds$basesQ20_suspect,
        sn_fail = qc_thresholds$sn_fail, sn_suspect = qc_thresholds$sn_suspect
      ),
      #columns = coldefs(),
      # ADD THIS rowStyle TO DRAW A LINE ABOVE THE FOOTER
      rowStyle = function(index) {
        if (index == nrow(df_table1_data())) {
          list(borderBottom = "1px solid grey") #
        }
      },
      #################
      details = function(index) {
        # Define unique IDs for the chromatogram
        chrom_output_id <- paste0('raw_', index)
        detail_container_id <- paste0('detail_raw_', index)
        
        # Render the chromatogram when the row is expanded
        local({
          output[[chrom_output_id]] <- renderChromatogram({
            raw_abif_data <- df2()[index, ]$data
            # subsampling
            MAX_POINTS_FOR_PLOT <- 2000
            step_size <- ceiling(length(raw_abif_data$data$DATA.1) / MAX_POINTS_FOR_PLOT)
            subsample_index <- seq(1, length(raw_abif_data$data$DATA.1), by = step_size)
            
            if (length(raw_abif_data) > 0) {
              # Format data for D3 - raw signal only
              list(
                traces = list(
                  DATA1 = raw_abif_data$data$DATA.1[subsample_index],
                  DATA2 = raw_abif_data$data$DATA.2[subsample_index],
                  DATA3 = raw_abif_data$data$DATA.3[subsample_index],
                  DATA4 = raw_abif_data$data$DATA.4[subsample_index]
                ),
                # No bases or quality scores for raw signal view
                plotWidth = 1100  # Smaller width than basecall view
              )
            }
          })
        })
        
        # Return the container with our custom chromatogram output
        htmltools::div(
          id = detail_container_id,
          style = "padding: 10px; padding-top: 0px; padding-bottom: 10px;",
          htmltools::div(
            class = "scrollable-plot-container",
            style = "overflow-x: auto; width: 100%; white-space: nowrap;",
            chromatogramOutput(chrom_output_id, width = "5000px", height = "250px")
          )
        )
      }
    )
  })
  
  # Sequence viewer output
  output$table3 <- renderReactable({
    reactable(
      df_table1_data(),
      pagination = FALSE, searchable = TRUE, highlight = TRUE, bordered = TRUE, striped = FALSE, compact = TRUE, resizable = TRUE,
      style = list(fontSize = "14px"),
      defaultColDef = colDef(footerStyle = list(color='grey', fontWeight = 'normal')),
      onClick = "expand", # Expand row details on click
      columns = create_qc_coldefs(
        nsamples = nrow(df_table1_data()),
        crl_fail = qc_thresholds$crl20_fail, crl_suspect = qc_thresholds$crl20_suspect,
        qv_fail = qc_thresholds$basesQ20_fail, qv_suspect = qc_thresholds$basesQ20_suspect,
        sn_fail = qc_thresholds$sn_fail, sn_suspect = qc_thresholds$sn_suspect
      ),
      rowStyle = function(index) {
        if (index == nrow(df_table1_data())) {
          list(borderBottom = "1px solid grey")
        }
      },
      details = function(index) {
        # Define unique IDs for sequence viewer
        seq_output_id <- paste0('seq_', index)
        detail_container_id <- paste0('detail_seq_', index)
        
        # Render sequence when the row is expanded
        local({
          output[[seq_output_id]] <- renderSequence({
            data_row <- df2()[index, ]
            raw_abif_data <- data_row$data
            
            if (length(raw_abif_data) == 0 || is.null(raw_abif_data$data$PBAS.1)) {
              warning("No sequence data available")
              return(NULL)
            }
            
            # Process and format data for sequence viewer
            bases <- tryCatch({
              strsplit(raw_abif_data$data$PBAS.1, "")[[1]]
            }, error = function(e) {
              warning("Error processing bases: ", e$message)
              NULL
            })
            
            if (is.null(bases) || length(bases) == 0) {
              warning("No base sequence found")
              return(NULL)
            }
            
            # Ensure bases and qscores have the same length by truncating to the shorter
            qscores_vec <- as.numeric(raw_abif_data$data$PCON.1)
            min_len <- min(length(bases), length(qscores_vec))
            if (min_len < 1) {
              warning("Empty sequence after processing")
              return(NULL)
            }
            if (min_len < length(bases) || min_len < length(qscores_vec)) {
              warning(sprintf("Truncating sequence to min length %d (bases:%d, qscores:%d)",
                              min_len, length(bases), length(qscores_vec)))
            }
            # Clamp CRL bounds to the new length
            crl_start_val <- as.numeric(data_row$crl_start) + 1
            crl_end_val <- as.numeric(data_row$crl_end) + 1
            crl_start_val <- min(crl_start_val, min_len)
            crl_end_val <- min(crl_end_val, min_len)
            if (crl_start_val > crl_end_val) crl_start_val <- crl_end_val

            list(
              bases = bases[seq_len(min_len)],
              qscores = qscores_vec[seq_len(min_len)],
              crlStart = crl_start_val,
              crlEnd = crl_end_val
            )
          })
        })
        
        # Return the container with sequence viewer
        htmltools::div(
          id = detail_container_id,
          style = "padding: 10px; padding-top: 0px; padding-bottom: 10px;",
          sequenceOutput(seq_output_id)
        )
      }
    )
  })
  
  #CRL Traces
  output$table4 <- renderReactable({
    data <- df2() %>%
      select('sample', 'seq', 'qscores', 'crl20', 'crl_start', 'crl_end', 'rawSeqLen')
    #qc_thresholds$crl_window_size, qval = qc_thresholds$crl_qv_threshold
    reactable(
      data,
      pagination = FALSE, searchable = TRUE, highlight = TRUE, bordered = TRUE, 
      striped = FALSE, compact = TRUE, resizable = TRUE, wrap = FALSE,
      style = list(fontSize = "14px"),
      columns = list(
        sample = colDef(minWidth = 200, vAlign = 'bottom'),
        seq = colDef(show = FALSE),
        crl_start = colDef(show = F),
        crl_end = colDef(show = F),
        rawSeqLen = colDef(show = F),
        crl20 = colDef(name = paste0("CRL",  qc_thresholds$crl_qv_threshold), minWidth = 70),
        qscores = colDef(
          vAlign = 'bottom',
          name = paste0("QV roll mean (CRL", qc_thresholds$crl_qv_threshold, ")"),
          cell = function(value, index) {
            sample_name <- data[index,]$sample
            crlstart <- data[index,]$crl_start
            crlend <- data[index,]$crl_end
            
            sp1 <- sparkline(
              #value,
              round(RcppRoll::roll_mean(value, n = qc_thresholds$crl_window_size, by = 1, na.rm = T, fill = c(0,0,0)), 0), # fill to avoid crashes with rawSeqLen < window
              type = 'line', 
              lineColor = "#e6550d",       # 
              normalRangeColor = "lightgrey", #
              #normalRangeMin = 0,
              #normalRangeMax = qc_thresholds$crl_qv_threshold,
              width = 800, height = 45,
              chartRangeMin = 1,
              chartRangeMax = 70,
              fillColor = NA, 
              lineWidth = 3,
              #tooltipFormatter = js_formatter
              tooltipFormat = paste0(
                '<a style = "font-size:12px; text-align: right;">',
                '<i>',sample_name, '</i><br>',
                'CRL start-end: <b>', crlstart, '-', crlend,'</b><br>',
                'Pos: <b>{{x}}</b><br>Roll mean QV: <b>{{y}}</b></a>'
                )
            )
            # find out where to place the bars
            a <- rep(0, data[index, ]$crl_start)
            b <- rep(70, data[index, ]$crl_end - data[index, ]$crl_start)
            c <- rep(0, data[index, ]$rawSeqLen - data[index, ]$crl_end)
            #
            sp2 <- sparkline(
              c(a, 70, b, 70, c),
              type = 'bar', barColor = 'rgb(161,217,155, 0.25)', zeroColor = 'grey', disableTooltips = FALSE, width = 800, height = 45
            )
            spk_composite(sp2,sp1)
          },
          minWidth = 700
        )
      )
      #style = list(fontSize = "14px")
    )
  })
  
  # Run info
  output$table5 <- renderReactable({
    data <- df() %>%
      select('sample', 'containerID','rundate', 'instrument', 'machine', 'well', 'capillary', 'gel_type', 'analysis_prot', 'data_coll_modfile', 'dyeset_name')
    reactable(
      data, 
      pagination = FALSE, searchable = TRUE, highlight = TRUE, 
      bordered = TRUE, striped = FALSE, compact = TRUE,
      wrap = FALSE, resizable = TRUE,
      style = list(fontSize = "14px"),
      columns = list(
        sample = colDef(minWidth = 150),
        containerID = colDef(minWidth = 70),
        rundate = colDef(minWidth = 70),
        instrument = colDef(minWidth = 50),
        well = colDef(minWidth = 50),
        capillary = colDef(minWidth = 50)
      )
    )
    
  })
  
}

shinyApp(ui, server)