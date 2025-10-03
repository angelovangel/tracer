library(shiny)
library(bslib)
library(dplyr)
library(reactable)
library(bsicons)
library(purrr)
library(sparkline)
library(writexl) # 
#library(ggrastr)
#library(sangerseqR)

source('global.R')
source('bases_as_html.R')

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
      # CSS to make the chromatogram plot container horizontally scrollable
      HTML("
        .scrollable-plot-container {
          overflow-x: auto;
          width: 100%; /* Ensures the container takes up the full width of the detail row */
          white-space: nowrap; /* Important to keep the content on one line for horizontal scrolling */
        }
        /* START OF NEW GLOBAL CSS FOR TOOLTIP WIDTH */
        .jqstooltip {
          /* Use !important to override any conflicting styles from Shiny/Bootstrap/Theme */
          width: 100px !important; 
          min-width: 100px !important;
          /* You might also need to ensure it doesn't try to inherit 100% width */
          box-sizing: content-box !important;
        }
        /* END OF NEW GLOBAL CSS */
      ")
    )
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
            "qc_trimMeanQscore",
            "Trim Qscore thresholds (fail, suspect)",
            min = 0,
            max = 60,
            value = c(
              qc_thresholds$trimMeanQscore_fail,
              qc_thresholds$trimMeanQscore_suspect
            ),
            step = 1,
            width = "100%"
          )
        ),
        footer = tagList(
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
    trimMeanQscore_fail = 30,
    trimMeanQscore_suspect = 40
  )
  
  # Update thresholds when Apply is clicked
  observeEvent(input$apply_qc_settings, {
    qc_thresholds$crl_window_size <- max(5, input$crl_window_size)
    qc_thresholds$crl_qv_threshold <- max(5, input$crl_qv_threshold)
    qc_thresholds$crl20_fail <- input$qc_crl20[1]
    qc_thresholds$crl20_suspect <- input$qc_crl20[2]
    qc_thresholds$basesQ20_fail <- input$qc_basesQ20[1]
    qc_thresholds$basesQ20_suspect <- input$qc_basesQ20[2]
    qc_thresholds$trimMeanQscore_fail <- input$qc_trimMeanQscore[1]
    qc_thresholds$trimMeanQscore_suspect <- input$qc_trimMeanQscore[2]
    removeModal()
  })
  
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
        crl20 = crl(qscores = unlist(qscores), window_size = qc_thresholds$crl_window_size, qval = qc_thresholds$crl_qv_threshold)$crl,
        crl_start = crl(qscores = unlist(qscores), window_size = qc_thresholds$crl_window_size, qval = qc_thresholds$crl_qv_threshold)$crl_start,
        crl_end = crl(qscores = unlist(qscores), window_size = qc_thresholds$crl_window_size, qval = qc_thresholds$crl_qv_threshold)$crl_end
      )
    # window size and QV threshold
  })
  
  # Reactive data frame for table1, which is also used for download
  df_table1_data <- reactive({
    req(df2()) # Ensure data exists
    df2() %>%
      select('sample', 'well', 'rawSeqLen', 'crl20', 'basesQ20', 'trimMeanQscore', 'data') %>%
      rowwise() %>%
      mutate(
        #sample = ifelse(input$shorten_names, str_trunc(sample, width = 42), sample),
        trimMeanQscore = round(trimMeanQscore),
        QC_flag = case_when(
          crl20 < qc_thresholds$crl20_fail |
            basesQ20 < qc_thresholds$basesQ20_fail |
            trimMeanQscore < qc_thresholds$trimMeanQscore_fail ~ "fail",
          crl20 < qc_thresholds$crl20_suspect |
            basesQ20 < qc_thresholds$basesQ20_suspect |
            trimMeanQscore < qc_thresholds$trimMeanQscore_suspect ~ "suspect",
          TRUE ~ "pass"
        )
      ) %>%
      ungroup() # Convert back to a standard data frame
  })
  
  # Create a reusable list of QC col defs with conditional formatting, it is for df2() that is used 3 times
  create_qc_coldefs <- function(crl_fail, crl_suspect, qv_fail, qv_suspect, trim_fail, trim_suspect) {
    list(
      # Sample name
      sample = colDef(
        minWidth = 200, 
        #footer = paste0("Total ", nrow(df2()), " samples")),
        footer = paste0("Total ", "samples")),
      # Other
      rawSeqLen = colDef(minWidth = 50),
      well = colDef(minWidth = 30, show = FALSE),
      data = colDef(show = FALSE),
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
      trimMeanQscore = colDef(
        name = "Qscore",
        minWidth = 50,
        html = T,
        footer = function(values) {
          paste0(
            "Min: ", round(min(values), 0), "<br>",
            "Max: ", round(max(values), 0), "<br>", 
            "Mean: ", round(qscore_mean(values), 0))
        },
        style = function(value) {
          if (is.na(value)) return(list(color = "#eee"))
          if (value < trim_fail) return(list(color = "#F44336", fontWeight = "normal"))
          if (value < trim_suspect) return(list(color = "#FFC107", fontWeight = "normal"))
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
  coldefs <- reactive({
    create_qc_coldefs(
      crl_fail = qc_thresholds$crl20_fail, crl_suspect = qc_thresholds$crl20_suspect, 
      qv_fail = qc_thresholds$basesQ20_fail, qv_suspect = qc_thresholds$basesQ20_suspect, 
      trim_fail = qc_thresholds$trimMeanQscore_fail, trim_suspect = qc_thresholds$trimMeanQscore_suspect
    )
  })
  
  # Basecall chromatogram table 
  output$table1 <- renderReactable({
    
    reactable(
      df_table1_data(),
      pagination = FALSE, searchable = TRUE, highlight = TRUE, bordered = TRUE, striped = FALSE, compact = TRUE, resizable = TRUE,
      style = list(fontSize = "14px"),
      defaultColDef = colDef(footerStyle = list(color='grey', fontWeight = 'normal')),
      onClick = "expand", # Expand row details on click
      columns = coldefs(),
      # ADD THIS rowStyle TO DRAW A LINE ABOVE THE FOOTER
      rowStyle = function(index) {
        if (index == nrow(df_table1_data())) {
          list(borderBottom = "1px solid grey") #
        }
      },
      #################
      details = function(index) {
        # 1. Define a unique ID for the plotOutput for this row and the detail container 
        plot_output_id <- paste0('chrom_', index)
        detail_container_id <- paste0('detail_chrom', index)

        # 2. Render the plot inside the details row when it's opened
        local({
          output[[plot_output_id]] <- renderPlot({
          raw_abif_data <- df2()[index, ]$data
          
          if (length(raw_abif_data) > 0) {
            plot_abif_chromatogram(raw_abif_data, type = 'basecall')
          } else {
            ggplot() + labs(title = "No ABIF raw data found for this sample.")
          }
          }
          #width = 12000,
          #height = 250
          )
        })

        # 3. Return the HTML structure containing the plotOutput, scrolling container
        htmltools::div(
          id = detail_container_id, # Assign a unique ID to scroll to
          # Keep padding-bottom to prevent table footer overlap
          style = "padding: 10px; padding-top: 0px; padding-bottom: 10px;",
          htmltools::div(
          class = "scrollable-plot-container",
            plotOutput(plot_output_id, width = "12000px", height = "250px")
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
                paste0("Qscore: fail < ", qc_thresholds$trimMeanQscore_fail, ", suspect < ", qc_thresholds$trimMeanQscore_suspect)
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
        # Rename columns to be descriptive in the Excel file
        rename(
          !!paste0("CRL", qc_thresholds$crl_qv_threshold) := crl20,
          QV20_plus = basesQ20,
          Trimmed_Mean_Qscore = trimMeanQscore,
          Raw_Seq_Length = rawSeqLen
        ) %>%
        # Select and order the final columns
        select(
          Sample = sample, 
          Raw_Seq_Length, 
          !!paste0("CRL", qc_thresholds$crl_qv_threshold), 
          QV20_plus, 
          Trimmed_Mean_Qscore, 
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
          "Qscore Fail Threshold",
          "Qscore Suspect Threshold"
        ),
        Value = c(
          qc_thresholds$crl_window_size,
          qc_thresholds$crl_qv_threshold,
          qc_thresholds$crl20_fail,
          qc_thresholds$crl20_suspect,
          qc_thresholds$basesQ20_fail,
          qc_thresholds$basesQ20_suspect,
          qc_thresholds$trimMeanQscore_fail,
          qc_thresholds$trimMeanQscore_suspect
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
      columns = coldefs(),
      # ADD THIS rowStyle TO DRAW A LINE ABOVE THE FOOTER
      rowStyle = function(index) {
        if (index == nrow(df_table1_data())) {
          list(borderBottom = "1px solid grey") #
        }
      },
      #################
      details = function(index) {
        # 1. Define a unique ID for the plotOutput for this row and the detail container 
        plot_output_id <- paste0('raw_', index)
        detail_container_id <- paste0('detail_raw_', index)
        
        # 2. Render the plot inside the details row when it's opened
        local({
          output[[plot_output_id]] <- renderPlot({
            raw_abif_data <- df2()[index, ]$data
            
            if (length(raw_abif_data) > 0) {
              plot_abif_chromatogram(raw_abif_data, type = 'rawsignal')
            } else {
              ggplot() + labs(title = "No ABIF raw data found for this sample.")
            }
          }
          #width = 12000,
          #height = 250
          )
        })
        
        # 3. Return the HTML structure containing the plotOutput, scrolling container
        htmltools::div(
          id = detail_container_id, # Assign a unique ID to scroll to
          # Keep padding-bottom to prevent table footer overlap
          style = "padding: 10px; padding-top: 0px; padding-bottom: 10px;",
          htmltools::div(
            class = "scrollable-plot-container",
            plotOutput(plot_output_id, width = "1200px", height = "250px")
          )
        )
      }
    )
    #make_qc_table(fdata = df_table1_data(), ftype = 'rawsignal', fwidth = 1200, fheight = 200)
  })
  
  # HTML sequence
  output$table3 <- renderReactable({
    data <-  df_table1_data()
    
    reactable(
      data,
      pagination = FALSE, searchable = TRUE, highlight = TRUE, bordered = TRUE, striped = FALSE, compact = TRUE, resizable = TRUE,
      style = list(fontSize = "14px"),
      defaultColDef = colDef(footerStyle = list(color='grey', fontWeight = 'normal')),
      onClick = "expand", # Expand row details on click
      columns = coldefs(),
      # ADD THIS rowStyle TO DRAW A LINE ABOVE THE FOOTER
      rowStyle = function(index) {
        if (index == nrow(df_table1_data())) {
          list(borderBottom = "1px solid grey") #
        }
      },
      
      details = colDef(
        name = "",
        details = function(index){
          bases <- data[index,]$data$data$PBAS.1 %>% str_split('') %>% unlist
          qscores <- data[index,]$data$data$PCON.1
          content <-format_bases_as_html(bases, qscores)
          content
        },
        html = TRUE
      )
    )
      
  })
  
  #CRL Traces
  output$table4 <- renderReactable({
    data <- df2() %>%
      select('sample', 'seq', 'qscores', 'crl20', 'crl_start', 'crl_end', 'rawSeqLen')
    #qc_thresholds$crl_window_size, qval = qc_thresholds$crl_qv_threshold
    reactable(
      data,
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
            sp1 <- sparkline(
              round(RcppRoll::roll_mean(value, n = qc_thresholds$crl_window_size, by = 1), 0),
              type = 'line', 
              lineColor = "darkred",       # Color for values outside normalRange (i.e., > crl_qv_threshold)
              normalRangeColor = "lightgrey", # Color for values inside normalRange (i.e., <= crl_qv_threshold)
              normalRangeMin = 0,
              normalRangeMax = qc_thresholds$crl_qv_threshold,
              width = 800, height = 50,
              chartRangeMin = 1,
              chartRangeMax = 70,
              fillColor = NA, 
              lineWidth = 3,
              #tooltipFormatter = js_formatter
              tooltipFormat = 'Position: <b>{{x}}</b><br>Roll mean QV: <b>{{y}}</b>'
            )
            # find out where to place the bars
            a <- rep(0, data[index, ]$crl_start)
            b <- rep(0, data[index, ]$crl_end - data[index, ]$crl_start)
            c <- rep(0, data[index, ]$rawSeqLen - data[index, ]$crl_end)
            #
            sp2 <- sparkline(
              c(a, 1, b, 1, c),
              type = 'bar', barColor = 'black', zeroColor = 'grey', disableTooltips = TRUE
            )
            spk_composite(sp1, sp2, options = list(width = 800, height = 50))
          },
          minWidth = 800
        )
      ),
      pagination = FALSE, searchable = TRUE, highlight = FALSE,
      bordered = F, striped = FALSE, compact = F,
      wrap = FALSE, resizable = TRUE
      #style = list(fontSize = "14px")
    )
  })
  
  # Run info
  output$table5 <- renderReactable({
    data <- df() %>%
      select('sample', 'rundate', 'instrument', 'machine', 'well', 'capillary', 'gel_type', 'analysis_prot', 'data_coll_modfile', 'dyeset_name')
    reactable(
      data, 
      pagination = FALSE, searchable = TRUE, highlight = TRUE, 
      bordered = TRUE, striped = FALSE, compact = TRUE,
      wrap = FALSE, resizable = TRUE,
      style = list(fontSize = "14px"),
      columns = list(
        sample = colDef(minWidth = 200),
        rundate = colDef(minWidth = 70),
        instrument = colDef(minWidth = 70),
        well = colDef(minWidth = 50),
        capillary = colDef(minWidth = 50)
      )
    )
    
  })
  
}

shinyApp(ui, server)