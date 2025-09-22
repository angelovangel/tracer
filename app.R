library(shiny)
library(bslib)
library(dplyr)
library(reactable)
library(bsicons)
library(purrr)
library(sparkline)
#library(sangerseqR)

source('global.R')

sidebar <- sidebar(
  tags$div(
    style = "padding: 20px;",
    tags$h4("Upload Sanger AB1 Files", style = "margin-bottom: 15px;"),
    fileInput(
      "ab1",
      label = NULL,
      accept = ".ab1",
      multiple = TRUE,
      buttonLabel = "Browse...",
      placeholder = "No files selected"
    ),
    tags$hr(),
    tags$p(
      "Select one or more .ab1 files to view summary statistics and traces.",
      style = "font-size: 14px; color: #555;"
    ),
    checkboxInput('shorten_names', 'Shorten names', value = F),
    tags$br(),
    actionButton('settings', 'QC settings'),
    tags$br(),
    tags$span(
      bsicons::bs_icon("info-circle"),
      " Supported: Sanger AB1 files",
      style = "color: #888; font-size: 13px;"
    )
  )
)

ui <- page_navbar(
  sidebar = sidebar,
  title = "",
  nav_panel('QC flags',
    reactableOutput('table1'),
    htmlOutput('qc_footer')
  ),
  nav_panel('Run info', 
    reactableOutput('table2')
  ),
  nav_panel('QC Summary', 'Under construction'),
  theme = bs_theme(bootswatch = "simplex")
  #
)

server <- function(input, output, session) {

  # Show modal when QC settings button is clicked
  observeEvent(input$settings, {
    showModal(modalDialog(size = 'l',
      title = "QC Settings",
      tags$p('CRL settings', style = "font-weight: bold; font-size: 15px;"),
      tags$div(
        style = "display: flex; gap: 16px; align-items: center; font-size: 13px; color: #1C398E;",
        numericInput(
          'crl_window_size',
          label = 'CRL window size',
          min = 5,
          max = 100,
          value = max(5, min(100, qc_thresholds$crl_window_size)),
          width = "110px"
        ),
        numericInput(
          'crl_qv_threshold',
          label = 'CRL QV threshold',
          min = 5,
          max = 50,
          value = max(5, min(50, qc_thresholds$crl_qv_threshold)),
          width = "110px"
        ),
      ),
      tags$hr(style = "border-top: 1px solid #D6D3D1; margin: 16px 0;"),
      
      tags$p("Thresholds for sample QC flagging:", style = "font-weight: bold; font-size: 15px;"),
      tags$div(
        style = "width: 100%; font-size: 12px; color: #1C398E;",
        sliderInput("qc_crl20", "CRL thresholds (fail, suspect)", min = 0, max = 1500, value = c(qc_thresholds$crl20_fail, qc_thresholds$crl20_suspect), step = 10, width = "100%")
      ),
      tags$div(
        style = "width: 100%; font-size: 12px; color: #1C398E;",
        sliderInput("qc_basesQ20", "Q20+ thresholds (fail, suspect)", min = 0, max = 1500, value = c(qc_thresholds$basesQ20_fail, qc_thresholds$basesQ20_suspect), step = 10, width = "100%")
      ),
      tags$div(
        style = "width: 100%; font-size: 12px; color: #1C398E;",
        sliderInput("qc_trimMeanQscore", "Trim Qscore thresholds (fail, suspect)", min = 0, max = 60, value = c(qc_thresholds$trimMeanQscore_fail, qc_thresholds$trimMeanQscore_suspect), step = 1, width = "100%")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("apply_qc_settings", "Apply")
      ),
      easyClose = TRUE
    ))
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
    df() %>%
      rowwise() %>%
      mutate(crl20 = crl(qscores = unlist(qscores), window_size = qc_thresholds$crl_window_size, qval = qc_thresholds$crl_qv_threshold))
     # window size and QV threshold
  })
  
  output$table1 <- renderReactable({
    data <- df2() %>%
      select('sample', 'well', 'rawSeqLen', 'crl20', 'basesQ20', 'trimMeanQscore') %>%
      rowwise() %>%
      mutate(
        sample = ifelse(input$shorten_names, str_trunc(sample, width = 42), sample),
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
      )
    reactable(
      data, pagination = FALSE, searchable = TRUE, highlight = TRUE, bordered = TRUE, striped = FALSE, compact = TRUE,
      columns = list(
        sample = colDef(minWidth = 200),
        rawSeqLen = colDef(minWidth = 50),
        well = colDef(minWidth = 30),
        crl20 = colDef(
          name = "CRL",
          minWidth = 50,
          style = function(value) {
            if (is.na(value)) return(list(color = "#eee"))
            if (value < qc_thresholds$crl20_fail) return(list(color = "#F44336", fontWeight = "bold"))
            if (value < qc_thresholds$crl20_suspect) return(list(color = "#FFC107", fontWeight = "bold"))
            list(color = "#4CAF50", fontWeight = "bold")
          }
        ),
        basesQ20 = colDef(
          name = "QV20+",
          minWidth = 50,
          style = function(value) {
            if (is.na(value)) return(list(color = "#eee"))
            if (value < qc_thresholds$basesQ20_fail) return(list(color = "#F44336", fontWeight = "bold"))
            if (value < qc_thresholds$basesQ20_suspect) return(list(color = "#FFC107", fontWeight = "bold"))
            list(color = "#4CAF50", fontWeight = "bold")
          }
        ),
        trimMeanQscore = colDef(
          name = "Qscore",
          minWidth = 50,
          style = function(value) {
            if (is.na(value)) return(list(color = "#eee"))
            if (value < qc_thresholds$trimMeanQscore_fail) return(list(color = "#F44336", fontWeight = "bold"))
            if (value < qc_thresholds$trimMeanQscore_suspect) return(list(color = "#FFC107", fontWeight = "bold"))
            list(color = "#4CAF50", fontWeight = "bold")
          }
        ),
        QC_flag = colDef(
          name = "QC flag",
          minWidth = 80,
          cell = function(value) {
            color <- switch(
              value,
              "fail" = "#F44336",
              "suspect" = "#FFC107",
              "pass" = "#4CAF50",
              "#eee"
            )
            htmltools::tagList(
              htmltools::tags$span(
                style = paste0(
                  "display: inline-block; width: 14px; height: 14px; border-radius: 50%; background:", color, "; margin-right: 8px; vertical-align: middle;"
                ),
                ""
              ),
              value
            )
          }
        )
      )
    ) 
  })

  output$qc_footer <- renderUI({
  req(df())
  tags$details(
    open = NA, # This makes the details expanded by default
    style = "margin-top: 10px;",
    tags$summary(
      style = "font-size: 13px; color: #37474f; cursor: pointer;",
      "QC settings and thresholds used"
    ),
    tags$div(
      style = "background: #f5f5f5; border-radius: 6px; padding: 10px 16px; font-size: 12px; color: #37474f;",
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
})
  
  output$table2 <- renderReactable({
    data <- df() %>%
      select('sample', 'rundate', 'instrument', 'machine', 'capillary', 'analysis_prot', 'data_coll_modfile', 'dyeset_name', 'polymer_expdate')
    reactable(
      data, 
      pagination = FALSE, searchable = TRUE, highlight = TRUE, 
      bordered = TRUE, striped = FALSE, compact = TRUE,
      style = list(fontSize = "14px") # Decrease font size
    )
      
  })

}

shinyApp(ui, server)